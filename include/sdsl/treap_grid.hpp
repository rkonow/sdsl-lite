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
/*! \file treap_grid.hpp
    \brief treap_grid.hpp contains an implementation of a grid representation using a compact treap
    \author Simon Gog, Roberto Konow
*/
//#ifndef INCLUDED_SDSL_TREAP_GRID
//#define INCLUDED_SDSL_TREAP_GRID
#pragma  once

#include "sdsl/bp_support_sada.hpp"
#include "sdsl/dac_vector.hpp"
#include "sdsl/int_vector.hpp"
#include "sdsl/int_vector_buffer.hpp"
#include "sdsl/bp_tree.hpp"
#include "bp_tree.hpp"
//#include "sdsl/wt_algorithm.hpp"

#include <iostream>
#include <memory>
#include <complex>
#include <algorithm>
#include <queue>
#include <array>

using std::cout;
using std::endl;

namespace sdsl
{
template<size_t  t_levels=0,
        typename t_bp=bp_support_sada<>,
        typename t_rmq=rmq_succinct_sct<false>,
        typename t_weight_vec=dac_vector<>,
        typename t_y_vec=dac_vector<>
        >

class treap_grid {
    struct treap_node;

    public:
        enum { permuted_x = false};
        typedef int_vector<>::size_type         size_type;
        typedef std::shared_ptr<treap_node>     node_ptr;
        typedef std::complex<uint64_t>          point_type;
        typedef std::pair<size_type,size_type>  range_type;



        class top_k_iterator
        {
            public:
                typedef void(*t_mfptr)();
                typedef std::pair<point_type, uint64_t> t_point_val;
                typedef std::array<uint64_t, 4> t_state;

            private:
                const treap_grid*            m_topk;
                std::priority_queue<t_state> m_pq;
                t_point_val                  m_point_val;
                bool                         m_valid = false;
            public:
                top_k_iterator() = default;
                top_k_iterator(const top_k_iterator&) = default;
                top_k_iterator(top_k_iterator&&) = default;
                top_k_iterator& operator=(const top_k_iterator&) = default;
                top_k_iterator& operator=(top_k_iterator&&) = default;

                top_k_iterator(const treap_grid& topk, point_type p1, point_type p2) :
                        m_topk(&topk), m_valid(topk.size()>0)
                {
                    if (m_topk->size() > 0) {
                        auto iv_it = map_to_sorted_sequence(*m_topk,
                                {real(p1), real(p2)}, {imag(p1), imag(p2)});
                        while (iv_it) {
                            auto x = std::get<0>(*iv_it);
                            auto y = std::get<1>(*iv_it);
                            auto w = std::get<2>(*iv_it);
                            auto max_w = std::get<3>(*iv_it);
                            m_pq.push({max_w, w, x, y});
                            ++iv_it;
                        }
                        ++(*this);
                    }
                }


                //! Prefix increment of the iterator
                top_k_iterator& operator++()
                {
                    m_valid = false;
                    if (!m_pq.empty()) {
                        auto s = m_pq.top();
                        m_point_val = { {s[2],s[3]}, s[1] };
                        m_pq.pop();
                        m_valid = true;
                    }
                    return *this;
                };

                //! Postfix increment of the iterator
                top_k_iterator operator++(int)
                {
                    top_k_iterator it = *this;
                    ++(*this);
                    return it;
                }

                t_point_val operator*() const
                {
                    return m_point_val;
                }

                operator t_mfptr() const
                {
                    return (t_mfptr)(m_valid);
                }
        };

    private:

        struct treap_node {
            public:
                size_type    m_x_value;
                size_type    m_y_value;
                size_type    m_y_dest_value;

                node_ptr     m_left;
                node_ptr     m_right;

                treap_node() = default;
                treap_node(const treap_node& tr) = default;
                treap_node& operator=(const treap_node& tr) = default;
                treap_node& operator=(treap_node&& tr) = default;
                treap_node(treap_node&& tr) = default;

                treap_node(size_t x_value, size_t y_value, size_t y_dest=0) : m_x_value(x_value),
                                                             m_y_value(y_value),
                                                             m_y_dest_value(y_dest),
                                                             m_left(nullptr),
                                                             m_right(nullptr)
                { }
        };

    public:
        // start_pos, end_pos, parent_node, left(0),right(1)
        typedef std::tuple<size_type, size_type, node_ptr, bool>  t_state;
        typedef std::array<size_type,6> node_data;
        typedef std::tuple<size_type,size_type, size_type, size_type> node_data_type;

        node_ptr                    m_root;
        bp_tree<t_bp>               m_tree;
        t_rmq                       m_rmq;
        t_weight_vec                m_weights;
        t_y_vec                     m_y_values;
        size_type                   m_size;
        size_type                   m_levels = 0;
        bit_vector                  m_bp;
        node_data_type              m_root_data; // fix this for serialize

        treap_grid() = default;
        treap_grid(const treap_grid& tr) = default;
        treap_grid& operator=(const treap_grid& tr) = default;
        treap_grid& operator=(treap_grid&& tr) = default;
        treap_grid(treap_grid&& tr) = default;

        size_type size() const
        {
            return m_size;
        }



        void swap(treap_grid& tr)
        {
            if (this != &tr) {
                m_bp.swap(tr.m_bp);
                m_rmq.swap(tr.m_rmq);
                m_weights.swap(tr.m_weights);
                m_y_values.swap(tr.m_y_values);
                m_tree.swap(tr.m_tree);
                tr.m_root_data = m_root_data;
                tr.m_size = m_size;
                // TODO: fix this swap thing...
//                tr.m_tree = bp_tree<t_bp>(m_bp);
            }
        }

        treap_grid(int_vector_buffer<>&,
                   int_vector_buffer<>& buf_y,
                   int_vector_buffer<>& buf_w)
        {
            m_levels = t_levels;
            std::string temp_grid_w_filename = "treap_grid_w" + std::to_string(util::pid())
                    + "_" + std::to_string(util::id());
            std::string temp_grid_y_filename =  "treap_grid_y_" + std::to_string(util::pid())
                    + "_" + std::to_string(util::id());

            // (1) we create the treap with pointers
            {
                std::stack<t_state> st;
                int_vector<> y_vec;
                load_from_file(y_vec, buf_y.filename());
                m_size = y_vec.size();
                size_t min_pos = findMin(0, m_size - 1, y_vec);
                m_root = node_ptr(new treap_node(min_pos, y_vec[min_pos], y_vec[min_pos]));
                st.emplace(0, min_pos - 1, m_root, 0);
                st.emplace(min_pos + 1, m_size - 1, m_root, 1);

                while (!st.empty()) {
                    t_state values = st.top();
                    st.pop();
                    node_ptr node = std::get<2>(values);
                    size_type start = std::get<0>(values);
                    size_type end = std::get<1>(values);
                    bool left = std::get<3>(values);
                    // TODO: clean up these IFs...
                    if (start > end)
                        continue;

                    if (start == std::numeric_limits<size_type>::max())
                        continue;

                    if (end == std::numeric_limits<size_type>::max())
                        continue;

                    if (start == end)
                        min_pos = start;
                    else
                        min_pos = findMin(start, end, y_vec);
                    node_ptr new_node = node_ptr(new treap_node(min_pos, y_vec[min_pos]-node->m_y_dest_value, y_vec[min_pos]));

                    if (left == 0) {
                        node->m_left = new_node;
                        node = node->m_left;
                    } else {
                        node->m_right = new_node;
                        node = node->m_right;
                    }
                    st.emplace(start, min_pos - 1, new_node, 0);
                    st.emplace(min_pos + 1, end, new_node, 1);
                }
            }
            // (2) traverse the trep with pointer in pre-order and save the y-dest values and weights
            {
                int_vector<> weights;
                load_from_file(weights, buf_w.filename());

                int_vector_buffer<> permuted_weights(temp_grid_w_filename, std::ios::out);
                int_vector_buffer<> permuted_y(temp_grid_y_filename, std::ios::out);
                std::stack<node_ptr> st;

                st.push(m_root);
                while(!st.empty()) {
                    node_ptr tmp = st.top();
                    permuted_weights.push_back(weights[tmp->m_x_value]);
                    permuted_y.push_back(tmp->m_y_value);
                    st.pop();
                    if (tmp->m_right != nullptr)
                        st.push(tmp->m_right);
                    if (tmp->m_left != nullptr)
                        st.push(tmp->m_left);
                }
                permuted_y.close(false);
                permuted_weights.close(false);
            }
            // (3) load the weights and y-dest values into the class.
            {
                int_vector<> weights;
                load_from_file(weights, temp_grid_w_filename);
                m_weights = t_weight_vec(weights);
                m_rmq = t_rmq(&m_weights);

                int_vector<> y_values;
                load_from_file(y_values, temp_grid_y_filename);
                m_y_values = t_y_vec(y_values);

            }
            // (4) we convert the treap into a balanced parenthesis representation and then delete the pointer version
            // we also keep track of the root values for traversal
            {
                convert_bp();
                delete_pointers();
                m_root_data = get_data(1, m_y_values[0]);

            }
            sdsl::remove(temp_grid_w_filename);
            sdsl::remove(temp_grid_y_filename);

        }

         range_type get_range(size_type node, size_type pos) const {
            size_t left = pos;
            size_t rank_param = m_tree.close(m_tree.enclose(node));
            size_t right = m_tree.preorder_rank(rank_param)-2;
            return {left,right};
        }


        range_type get_weight_data(size_type node, size_type pos) const {
            range_type rng = get_range(node, pos);
            size_type max_pos = m_rmq(rng.first,rng.second);
            size_type w_value = m_weights[pos];
            size_type w_max_value = m_weights[max_pos];
            return {w_max_value, w_value};
        }

        node_data_type get_data(size_type node, size_type y_old) const {
            size_type pos = m_tree.preorder_rank(node)-2;
            size_type y_value = y_old + m_y_values[pos];
            size_type bp_close =  m_tree.close(node);
            size_type x_value =  bp_close - m_tree.preorder_rank(bp_close);
            return {x_value, y_value,  node, pos};
        }

        node_data_type move_left(size_type node, size_type y_old) const {
            size_type left_node = m_tree.first_child(node);
            if (left_node != std::numeric_limits<size_t>::max()) {
                return get_data(left_node,y_old);
            } else {
                size_type max = std::numeric_limits<size_t>::max();
                return {max,max,max};
            }
        }

        node_data_type move_right(size_type node, size_type y_old) const {
            size_type left_node = m_tree.next_sibling(node);
            if (left_node != std::numeric_limits<size_t>::max()) {
                return get_data(left_node,y_old);
            } else {
                size_type max = std::numeric_limits<size_t>::max();
                return {max,max,max};
            }
        }


    void convert_bp() {
            // TODO: change vector<bool> to something of sdsl
            std::vector<bool> bp;
            _convert_bp(m_root, bp);
            m_bp = bit_vector(bp.size(), 0);
            // TODO: efficient way to do this?
            for (size_type i = 0; i < bp.size(); i++) {
                m_bp[i] = bp[i];
            }
            m_tree = bp_tree<t_bp>(&m_bp);
        }

        void delete_pointers() {
            if (m_root != nullptr) {
                _delete_pointers(m_root);
            }
        }

        size_type getLeft(size_type v) const {
            return m_tree.first_child(v);
        }

        size_type getRight(size_type v) const {
            return m_tree.next_sibling(v);
        }

        // TODO: auxiliary functions, we're going to delete these
        void print_inorder(std::ostream& os) const {
            if (m_root != nullptr) {
                _print_inorder(m_root, os);
            }
        }

        void print_preorder(std::ostream& os)  {
            if (m_root != nullptr) {
                _print_preorder(m_root,os);
            }
        }

        // TODO: change this to print bp representation, not the pointer representation
        friend std::ostream& operator<< (std::ostream& os, treap_grid& tg) {
            tg.print_preorder(os);
            return os;
        }

        top_k_iterator topk(point_type p1, point_type p2) const  {
            return top_k_iterator(*this, p1, p2);
        }

        //! Serializes the data structure into the given ostream
        size_type serialize(std::ostream& out, structure_tree_node* v=nullptr,
                std::string name="")const
        {
            structure_tree_node* child = structure_tree::add_child(
                    v, name, util::class_name(*this));
            size_type written_bytes = 0;
            written_bytes += write_member(m_size, out, child, "n");
            written_bytes += write_member(m_levels, out, child, "levels");
            written_bytes += m_tree.serialize(out, child, "bp");
            written_bytes += m_weights.serialize(out, child, "weights");
            written_bytes += m_y_values.serialize(out,child, "y_values");
            written_bytes += m_rmq.serialieze(out,child, "rmq");
            structure_tree::add_size(child, written_bytes);
            return written_bytes;
        }

        //! Loads the data structure from the given istream.
        void load(std::istream& in) {
            read_member(m_size, in);
            read_member(m_levels, in);
            m_tree.load(in);
            m_weights.load(in);
            m_y_values.load(in);
            m_rmq.load(in);
        }

    private:
        void _delete_pointers(node_ptr node) {
            if (node != nullptr) {
                _delete_pointers(node->m_left);
                _delete_pointers(node->m_right);
                node.reset();
            }
        }

        void _convert_bp(node_ptr node, std::vector<bool>& bitmap) {
            bitmap.push_back(true);
            while(node != nullptr) {
                _convert_bp(node->m_left, bitmap);
                node = node->m_right;
            }
            bitmap.push_back(false);
        }

        // TODO: this is probably super-inefficient
        // Change this to RMQ with sub-intervals...
        size_type findMin(size_type start, size_type end, int_vector<>& y_vec)
        {
            size_type mid = (start+end)/2;
            size_type min_distance = start-end;
            size_type min_pos = 0;
            size_type min_value = std::numeric_limits<size_type>::max();
            for (size_type i = start ; i <= end; i++) {
                if (y_vec[i] < min_value) {
                    min_pos = i;
                    min_value = y_vec[i];
                }
                if (y_vec[i] == min_value) {
                    size_type distance = size_type_abs(i,mid);
                    if (distance < min_distance) {
                        min_distance = distance;
                        min_pos = i;
                    }
                }
            }
            return min_pos;
        }

        void _print_inorder(node_ptr& node, std::ostream& os) const {
            if (node != nullptr) {
                _print_inorder(node->m_left, os);
                os << "(" << node->m_x_value << "," << node->m_y_value << ") ";
                _print_inorder(node->m_right, os);
            }
        }

        void _print_preorder(node_ptr& node, std::ostream& os) const {
            if (node != nullptr) {
                os << "(" << node->m_x_value << "," << node->m_y_value << ") ";
                _print_preorder(node->m_left, os);
                _print_preorder(node->m_right, os);
            }
        }

        size_type size_type_abs(size_type a, size_type b) {
            return a < b ? b - a : a - b;
        }
};
}
//#endif