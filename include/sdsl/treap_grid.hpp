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
#include <iostream>
#include <memory>

using std::cout;
using std::endl;
namespace sdsl
{
template<int levels=4,
        typename t_bp=bp_support_sada<>,
        typename t_rmq=rmq_succinct_sct<false>,
        typename t_weight_vec=dac_vector<>,
        typename t_y_vec=dac_vector<>
        >

class treap_grid {
    struct treap_node;

    public:
        typedef int_vector<>::size_type size_type;
        enum { permuted_x = false};
        typedef std::shared_ptr<treap_node> node_ptr;


    private:

        struct cursor {
            size_type    m_x_value;
            size_type    m_y_value;
            size_type    m_w_value;
            size_type    m_node;
            bool         disabled;

            cursor() = default;

            cursor(size_type m_x_value, size_type m_y_value, size_type m_w_value)
                    : m_x_value(m_x_value),
                      m_y_value(m_y_value),
                      m_w_value(m_w_value),
                      disabled(false)
            { }
        };
        struct treap_node {
            public:
                size_type    m_x_value;
                size_type    m_y_value;
                size_type    m_preorder;

                node_ptr     m_left;
                node_ptr     m_right;

                treap_node() = default;
                treap_node(const treap_node& tr) = default;
                treap_node& operator=(const treap_node& tr) = default;
                treap_node& operator=(treap_node&& tr) = default;
                treap_node(treap_node&& tr) = default;

                treap_node(size_t x_value, size_t y_value) : m_x_value(x_value),
                                                             m_y_value(y_value),
                                                             m_left(nullptr),
                                                             m_right(nullptr)
                { }
        };
    public:
        // start_pos, end_pos, parent_node, left(0),right(1)
        typedef std::tuple<size_type, size_type, node_ptr, bool>  t_state;
        typedef std::pair<size_type,size_type>  t_range;
        typedef std::array<size_type,4> t_data;

        node_ptr                    m_root;
        bp_tree<t_bp>               m_tree;
        t_rmq                       m_rmq;
        t_weight_vec                m_weights;
        t_y_vec                     m_y_values;
        size_type                   m_size;
        bit_vector                  m_bp;
        cursor                      m_cursor;

        int_vector<>                m_x_preorder;

        size_type size() const
        {
            return m_size;
        }

        treap_grid() = default;
        treap_grid(const treap_grid& tr) = default;
        treap_grid& operator=(const treap_grid& tr) = default;
        treap_grid& operator=(treap_grid&& tr) = default;
        treap_grid(treap_grid&& tr) = default;

        treap_grid(int_vector_buffer<>&,
                   int_vector_buffer<>& buf_y,
                   int_vector_buffer<>& buf_w)
        {

            // (1) we create the treap with pointers
            {
                std::stack<t_state> st;
                int_vector<> y_vec;
                load_from_file(y_vec, buf_y.filename());

                m_size = y_vec.size();
                size_t min_pos = findMin(0, m_size - 1, y_vec);
                m_root = node_ptr(new treap_node(min_pos, y_vec[min_pos]));
                st.emplace(0, min_pos - 1, m_root, 0);
                st.emplace(min_pos + 1, m_size - 1, m_root, 1);

                while (!st.empty()) {
                    t_state values = st.top();
                    st.pop();
                    node_ptr node = std::get<2>(values);
                    size_type start = std::get<0>(values);
                    size_type end = std::get<1>(values);
                    bool left = std::get<3>(values);

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

                    node_ptr new_node = node_ptr(new treap_node(min_pos, y_vec[min_pos]));

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
            // (2) we create the data structures for weights and y-values
            {
                int_vector<> weights;
                load_from_file(weights, buf_w.filename());
                m_rmq = t_rmq(&weights);
                m_weights = t_weight_vec(weights);

                int_vector<> y_values;
                load_from_file(y_values, buf_y.filename());
                m_y_values = t_y_vec(y_values);

            }
            convert_bp();

            size_type x_root = m_root->m_x_value;
            m_cursor = cursor(m_root->m_x_value, m_y_values[x_root], m_weights[x_root]);

        }

        void _convert_bp(node_ptr node, std::vector<bool>& bitmap) {
            bitmap.push_back(true);

            while(node != nullptr) {
                cout << node->m_x_value << " , " << node->m_y_value << endl;
                _convert_bp(node->m_left, bitmap);
                node = node->m_right;
            }

            bitmap.push_back(false);
        }

        void convert_bp() {
            std::vector<bool> bp;
            _convert_bp(m_root, bp);
            cout << endl;
            m_bp = bit_vector(bp.size(), 0);

            // efficient way to do this?
            for (size_type i = 0; i < bp.size(); i++) {
                m_bp[i] = bp[i];
            }
            m_tree = bp_tree<t_bp>(&m_bp);
        }

        size_type getLeft(size_type v) {
            return m_tree.first_child(v);
        }

        size_type getRight(size_type v) {
            return m_tree.next_sibling(v);
        }

//        void getLeft(size_type v) {
//            cout << "left..." << endl;
//            size_type node = m_tree.first_child(v);
//            cout << "node = " << node << endl;
//            m_cursor.m_node = node;
//            if (node == std::numeric_limits<size_type>::max()) {
//                m_cursor.disabled = true;
//                return;
//            }
//            size_type x_value = m_tree.preorder_rank(m_tree.close(v))-2;
//            m_cursor.m_x_value = x_value;
//            m_cursor.m_y_value = m_y_values[x_value];
//            m_cursor.m_w_value = m_weights[x_value];
//        }
//
//        void getRight(size_type v) {
//            cout << "right..." << endl;
//            size_type node =  m_tree.next_sibling(v);
//            cout << "node = " << node << endl;
//            m_cursor.m_node = node;
//
//            if (node == std::numeric_limits<size_type>::max()) {
//                m_cursor.disabled = true;
//                return;
//            }
//            size_type x_value = m_tree.preorder_rank(m_tree.close(v))-2;
//            m_cursor.m_x_value = x_value;
//            m_cursor.m_y_value = m_y_values[x_value];
//            m_cursor.m_w_value = m_weights[x_value];
//        }

//        size_type getX(size_type v) {
//            size_type pos = m_tree.preorden_rank(v)-1;
//            size_type x = m_tree.preorden_rank(m_tree.close(v));
//            cout << "x = " << x << endl;
//            cout << "pos = " << pos << endl;
//            return x;
//        }



        void print_inorder() {
            if (m_root != nullptr) {
                _print_inorder(m_root);
            }
        }

        void print_preorder() {
            if (m_root != nullptr) {
                _print_preorder(m_root);
            }
        }

        bp_tree<> getTree() {
            return m_tree;
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

    private:
        void _print_inorder(node_ptr& node) {
            if (node != nullptr) {
                _print_inorder(node->m_left);
                std::cout << "(" << node->m_x_value << "," << node->m_y_value << ")" << std::endl;
                _print_inorder(node->m_right);
            }
        }

        void _print_preorder(node_ptr& node) {
            if (node != nullptr) {
                std::cout << "(" << node->m_x_value << "," << node->m_y_value << ")" << std::endl;
                _print_preorder(node->m_left);
                _print_preorder(node->m_right);
            }
        }

        size_type size_type_abs(size_type a, size_type b) {
            return a < b ? b - a : a - b;
        }
};

}
//#endif