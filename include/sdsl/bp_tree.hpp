/* sdsl - succinct data structures library
    Copyright (C) 20014 Simon Gog

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see http://www.gnu.org/licenses/ .
*/
/*! \file bp_tree.hpp
    \brief bp_tree.hpp contains an implementation of a balanced parentheses tree.
    \author Simon Gog, Roberto Konow
*/
#ifndef INCLUDED_SDSL_BP_TREE
#define INCLUDED_SDSL_BP_TREE
//#pragma once
#include "bp_support.hpp"

namespace sdsl {


// TODO: Discarded the 'node' class as in louds, since it would require some extra-processing (pre_order_ranks...)
// This could be added and maybe using  template we could decide if we use it or not...

//template<typename bp_support_t=bp_support_sada<> >
//class bp_node {
//    typedef bp_support_t::size_type     size_type;
//private:
//    size_type m_nr;  // node number;
//    size_type m_pos; // position in the Bitmap;
//public:
//    const size_type& nr;
//    const size_type& pos;
//
//    bp_node(size_type f_nr=0, size_type f_pos=0):m_nr(f_nr), m_pos(f_pos),nr(m_nr), pos(m_pos)
//    {}
//    bool operator==(const bp_node& v)const {
//        return m_nr == v.m_nr and m_pos ==v.m_pos;
//    }
//
//    bool operator!=(const  bp_node& v)const {
//        return !(v==*this);
//    }
//};


template<typename t_bp_support=bp_support_sada<>,
        typename t_bit_vec = bit_vector>
class bp_tree {

    public:
        typedef bit_vector::size_type    size_type;
        typedef t_bp_support             bp_support_type;
        typedef t_bit_vec                bit_vector_type;

    private:
            bp_support_type              m_bp;
            bit_vector_type              m_bv;
    public:
        const bit_vector_type&           bv;

    public:

        bp_tree() : m_bp(), m_bv(), bv(m_bv) { }

        bp_tree(const bit_vector_type& bp) : m_bp(), bv(m_bv) {
            m_bv = bit_vector_type(std::move(bp));
            util::init_support(m_bp, &m_bv);
        }

        bp_tree(const bp_tree &lt) : bv(m_bv) {
            *this = lt;
        }

        bp_tree(bp_tree &&lt) : bv(m_bv) {
            *this = std::move(lt);
        }

        bp_tree &operator=(const bp_tree& lt)  {
            if (this != &lt) {
                m_bp = lt.m_bp;
                m_bv = lt.m_bv;
            }
            return *this;
        }

        bp_tree operator=(bp_tree &&lt) {
            if (this != &lt) {
                m_bv = std::move(lt.m_bv);
                m_bp = std::move(lt.m_bp);
                m_bp.set_vector(&m_bv);
            }
            return *this;
        }

        void swap(bp_tree& tree) {
            m_bv.swap(tree.m_bv);
            util::swap_support(m_bp, tree.m_bp, &m_bv, &(tree.m_bv));
        }

        size_type root() const {
            return 1;
        }

        size_type length() const {
            return m_bv.size();
        }

        size_type nodes() const {
            return (m_bp.size()+1)/2;
        }

        size_type first_child(const size_type i) const {
            if (!isleaf(i)) {
                return i+1;
            } else {
                return std::numeric_limits<size_type>::max();
            }
        }

        size_type next_sibling(const size_type i) const {
            size_type aux = m_bp.find_close(i)+1;
            if (m_bv[aux] == 1)
                return aux;
            else
                return std::numeric_limits<size_type>::max();
        }

        size_type parent(const size_type v) const {
            return m_bp.enclose(v);
        }

        size_type preorder_select(const size_type v) const {
            return m_bp.select(v);
        }

        size_type preorder_rank(const size_type v) const {
            return m_bp.rank(v);
        }

        size_type node_depth(const size_type v) const {
            return 2*m_bp.rank(v) - v;
        }

        size_type close(const size_type v) const {
            return m_bp.find_close(v);
        }

        size_type enclose(const size_type v) const {
            return m_bp.enclose(v);
        }

        bool isleaf(const size_type v) const {
            return m_bv[v + 1] == 0;
        }

        size_type subtree_size(const size_type v) const {
            return (m_bp.find_close(v) -  v) / 2;
        }

        size_type serialize(std::ostream &out, structure_tree_node *v = nullptr, std::string name = "") const {
            structure_tree_node *child = structure_tree::add_child(v, name, util::class_name(*this));
            size_type written_bytes = 0;
            m_bv.serialize(out, child, "bit_vector");
            m_bp.serialize(out,child, "bp_support");
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        void load(std::istream &in) {
            m_bv.load(in);
            m_bp.load(in,&m_bv);
            m_bp.set_vector(&m_bv);
        }
    };
}

#endif