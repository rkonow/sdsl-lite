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
//#ifndef INCLUDED_SDSL_BP_TREE
//#define INCLUDED_SDSL_BP_TREE
#pragma once
#include "bp_support.hpp"

namespace sdsl {

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


template<typename bp_support_t=bp_support_sada<>,
         typename rank_support_t=rank_support_v5<0>>
class bp_tree {

    public:
        typedef bp_support_t             bp_support_type;
        typedef rank_support_t           rank_support_type;
        typedef size_t                   size_type;

    private:
            bp_support_type                        m_bp;

    public:
        bp_tree(const bit_vector *bp = nullptr) {
            m_bp = bp_support_type(bp);
        }

        bp_tree(const bp_tree &lt) {
            *this = lt;
        }
        // default?
        bp_tree(bp_tree &&lt) {
            *this = std::move(lt);
        }

        bp_tree &operator=(const bp_tree &lt) {
            if (this != &lt) {
                m_bp = lt.m_bp;
            }
            return *this;
        }

        bp_tree operator=(bp_tree &&lt) {
            if (this != &lt) {
                m_bp = std::move(lt.m_bp);
            }
            return *this;

        }

        size_type root() const {
            return 0;
        }

        size_type length() const {
            return m_bp.size();
        }

        size_type nodes() const {
            return (m_bp.size()+1)/2;
        }

    //    size_type degree(const node_type &v) const {
    //
    //    }

        size_type first_child(const size_type i) const {
            if (!isleaf(i)) {
                return i+1;
            } else {
                return std::numeric_limits<size_type>::max();
            }
        }

        size_type next_sibling(const size_type i) const {
            size_type aux = m_bp.find_close(i)+1;
            if (m_bp.access(aux) == 1)
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

        size_type close(const size_type v) const {
            return m_bp.find_close(v);
        }

        bool isleaf(const size_type v) const {
            return m_bp.access(v + 1) == 0;
        }

        size_type subtree_size(const size_type v) const {
            return (m_bp.find_close(v) -  v) / 2;
        }

        size_type serialize(std::ostream &out, structure_tree_node *v = nullptr, std::string name = "") const {
            structure_tree_node *child = structure_tree::add_child(v, name, util::class_name(*this));
            size_type written_bytes = 0;
            m_bp.serialize(out, child, "bp_tree");
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        void load(std::istream &in) {
            m_bp.load(in);
        }
    };
}

//#endif