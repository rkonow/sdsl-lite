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
#ifndef INCLUDED_SDSL_TREAP_GRID
#define INCLUDED_SDSL_TREAP_GRID

#include "sdsl/bp_support_sada.hpp"
#include "sdsl/dac_vector.hpp"
#include "sdsl/int_vector.hpp"
#include "sdsl/int_vector_buffer.hpp"
#include "sdsl/bp_tree.hpp"
#include "sdsl/wavelet_trees.hpp"
#include "sdsl/wt_algorithm.hpp"
#include "bp_tree.hpp"

#include <iostream>
#include <memory>
#include <complex>
#include <algorithm>
#include <queue>
#include <array>

using std::cout;
using std::endl;

namespace sdsl {
template<size_t t_levels = 2,
        typename t_bp=bp_support_sada<>,
        typename t_rmq=rmq_succinct_sct<false>,
        typename t_weight_vec=dac_vector<>,
        typename t_y_vec=dac_vector<>
>

class treap_grid {

public:
    struct treap_node;

    enum {
        permuted_x = false
    };
    typedef int_vector<>::size_type size_type;
    typedef treap_node *node_ptr;
    typedef std::pair<size_type, size_type> point_type;
    typedef std::pair<size_type, size_type> range_type;


    class top_k_iterator {
    public:
        typedef void(*t_mfptr)();

        typedef std::pair<point_type, uint64_t> t_point_val;

    private:
        const treap_grid *m_topk;
        t_point_val m_point_val;
        bool m_valid = false;

    public:
        struct comp_cand {
            bool operator()(treap_node &a, treap_node &b) {
                return a.m_max_weight < b.m_max_weight;
            }
        };

        struct comp_res {
            bool operator()(treap_node &a, treap_node &b) {
                return a.m_weight < b.m_weight;
            }
        };

        std::priority_queue<treap_node, std::vector<treap_node>, comp_cand> cand;
        std::priority_queue<treap_node, std::vector<treap_node>, comp_res> res;
        point_type m_p1;
        point_type m_p2;

        top_k_iterator() = default;

        top_k_iterator(const top_k_iterator &) = default;

        top_k_iterator(top_k_iterator &&) = default;

        top_k_iterator &operator=(const top_k_iterator &) = default;

        top_k_iterator &operator=(top_k_iterator &&) = default;

        top_k_iterator(const treap_grid &topk, point_type p1, point_type p2) :
                m_topk(&topk), m_valid(topk.size() > 0), m_p1(p1), m_p2(p2) {
            cand.push(m_topk->m_root_data);
            ++(*this);
        }


        //! Prefix increment of the iterator
        top_k_iterator &operator++() {
            m_valid = false;
            while (!cand.empty()) {
                treap_node node = cand.top();
                cand.pop();
                if (node.m_y_value <= std::get<1>(m_p2)) {
                    if (node.m_x_value >= std::get<0>(m_p1)) {
                        treap_node left = m_topk->move_left(node);
                        if (left.m_x_value != std::numeric_limits<size_type>::max()) {
                            cand.push(left);
                        }
                    }
                    if (node.m_x_value <= std::get<0>(m_p2)) {
                        treap_node right = m_topk->move_right(node);
                        if (right.m_x_value != std::numeric_limits<size_type>::max()) {
                            cand.push(right);
                        }
                    }
                    if (node.m_x_value >= std::get<0>(m_p1) and node.m_x_value <= std::get<0>(m_p2)) {
                        res.push(node);
                        if (cand.top().m_max_weight <= res.top().m_weight) {
                            treap_node node = res.top();
                            m_point_val = {{node.m_x_value, node.m_y_value}, node.m_weight};
                            res.pop();
                            m_valid = true;
                            return *this;
                        }
                    }
                }
            }
            if (res.size() > 0) {
                treap_node node = res.top();
                m_point_val = {{node.m_x_value, node.m_y_value}, node.m_weight};
                res.pop();
                m_valid = true;
            }
            return *this;
        }


        //! Postfix increment of the iterator
        top_k_iterator operator++(int) {
            top_k_iterator it = *this;
            ++(*this);
            return it;
        }

        t_point_val operator*() const {
            return m_point_val;
        }

        operator t_mfptr() const {
            return (t_mfptr) (m_valid);
        }
    };

public:
    struct treap_node {
        public:
            size_type m_x_value;
            size_type m_y_value;
            size_type m_y_dest_value;
            size_type m_node;
            size_type m_pos;
            size_type m_weight;
            size_type m_max_weight;

            node_ptr m_left;
            node_ptr m_right;

            treap_node() = default;

            treap_node(const treap_node &tr) = default;

            treap_node &operator=(const treap_node &tr) = default;

            treap_node &operator=(treap_node &&tr) = default;

            treap_node(treap_node &&tr) = default;

            treap_node(size_t x_value, size_t y_value, size_t y_dest = 0) : m_x_value(x_value),
                                                                            m_y_value(y_value),
                                                                            m_y_dest_value(y_dest),
                                                                            m_left(nullptr),
                                                                            m_right(nullptr)
            { }

            treap_node(size_type m_x_value, size_type m_y_value, size_type m_y_dest_value, size_type m_node,
                    size_type m_pos, size_type m_weight, size_type m_max_weight) :
                    m_x_value(m_x_value),
                    m_y_value(m_y_value),
                    m_y_dest_value(m_y_dest_value),
                    m_node(m_node),
                    m_pos(m_pos),
                    m_weight(m_weight),
                    m_max_weight(m_max_weight),
                    m_left(nullptr),
                    m_right(nullptr) { }
        };

    // start_pos, end_pos, parent_node, left(0),right(1)
    typedef std::tuple<size_type, size_type, treap_node *, bool> t_state;

    node_ptr        m_root;
    bp_tree<t_bp>   m_tree;
    t_rmq           m_rmq;
    t_weight_vec    m_weights;
    t_y_vec         m_y_values;
    size_type       m_size;
    treap_node      m_root_data;

    treap_grid() = default;

    treap_grid(const treap_grid &tr) = default;

    treap_grid &operator=(const treap_grid &tr) = default;

    treap_grid &operator=(treap_grid &&tr) = default;

    treap_grid(treap_grid &&tr) = default;

    size_type size() const {
        return m_size;
    }

    void swap(treap_grid &tr) {
        if (this != &tr) {
            m_rmq.swap(tr.m_rmq);
            m_weights.swap(tr.m_weights);
            m_y_values.swap(tr.m_y_values);
            m_tree.swap(tr.m_tree);
            tr.m_root_data = m_root_data;
            tr.m_size = m_size;
        }
    }

    treap_grid(int_vector_buffer<> &,
            int_vector_buffer<> &buf_y,
            int_vector_buffer<> &buf_w) {
        std::string temp_grid_w_filename = "treap_grid_w" + std::to_string(util::pid())
                + "_" + std::to_string(util::id());
        std::string temp_grid_y_filename = "treap_grid_y_" + std::to_string(util::pid())
                + "_" + std::to_string(util::id());

        // (1) we create the treap with pointers
        {
            std::stack<t_state> st;
            int_vector<> y_vec;
            load_from_file(y_vec, buf_y.filename());
            wt_int<> wt_y = wt_int<>(buf_y, buf_y.size());
            m_size = y_vec.size();
            size_t min_pos = find_min(0, m_size - 1, wt_y);
            m_root = new treap_node(min_pos, y_vec[min_pos], y_vec[min_pos]);
            st.emplace(0, min_pos - 1, m_root, 0);
            st.emplace(min_pos + 1, m_size - 1, m_root, 1);
            treap_node *node;
            while (!st.empty()) {
                t_state values = st.top();
                st.pop();
                node = std::get<2>(values);
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
                    min_pos = find_min(start, end, wt_y, node->m_y_dest_value);

                if (node == nullptr) {
                    cout << "range fucked up= " << start << " , " << end << endl;
                }
                node_ptr new_node = new treap_node(min_pos, y_vec[min_pos] - node->m_y_dest_value, y_vec[min_pos]);

                if (left == 0) {
                    node->m_left = new_node;
                    node = node->m_left;
                } else {
                    node->m_right = new_node;
                    node = node->m_right;
                }
                st.emplace(start, min_pos - 1, node, 0);
                st.emplace(min_pos + 1, end, node, 1);

            }

        }
        // (2) traverse the trep with pointer in pre-order and save the y-dest values and weights
        {
            std::stack<treap_node *> st;
            int_vector<> weights;
            size_type node_id = 0;

            load_from_file(weights, buf_w.filename());
            int_vector_buffer<> permuted_weights(temp_grid_w_filename, std::ios::out);
            int_vector_buffer<> permuted_y(temp_grid_y_filename, std::ios::out);

            st.push(m_root);
            while (!st.empty()) {
                treap_node *tmp = st.top();
                tmp->m_pos = node_id;
                tmp->m_weight = weights[tmp->m_x_value];
                node_id++;
                permuted_weights.push_back(tmp->m_weight);
                permuted_y.push_back(tmp->m_y_value);
                st.pop();
                if (tmp->m_right != nullptr)
                    st.push(tmp->m_right);
                if (tmp->m_left != nullptr)
                    st.push(tmp->m_left);
            }
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
        // (4) we convert the treap into a balanced parenthesis representation and
        // we also keep track of the root values for traversal
        {
            bit_vector bv = bit_vector(2 * (m_size + 1), 0);
            size_type count = 0;
            _convert_bp(m_root, bv, count);
            m_tree = bp_tree<t_bp>(bv);
            m_root_data = *m_root;
        }
        // (5) delete all pointers up to level t_levels.
        {
            typedef std::pair<treap_node *, size_type> t_traverse;
            std::queue<t_traverse> traverse;
            size_type level = 0;
            traverse.push({&m_root_data, level});
            while (!traverse.empty()) {
                t_traverse data = traverse.front();
                treap_node *node = std::get<0>(data);
                size_type lev = std::get<1>(data);
                traverse.pop();
                if (lev >= t_levels) {
//                        node_ptr p = node_ptr(node);
                    _delete_pointers(node->m_left);
                    _delete_pointers(node->m_right);
                    node->m_left = nullptr;
                    node->m_right = nullptr;
                    set_weight_data(node);
                } else {
                    set_weight_data(node);
                    if (node->m_left != nullptr)
                        traverse.push({node->m_left, lev + 1});
                    if (node->m_right != nullptr)
                        traverse.push({node->m_right, lev + 1});
                }
            }
            delete m_root;

        }
        sdsl::remove(temp_grid_w_filename);
        sdsl::remove(temp_grid_y_filename);
    }

    ~treap_grid() {
        m_root = &m_root_data;
        _delete_pointers(m_root);
    }

    bp_tree<t_bp> getTree() {
        return m_tree;
    }

    void set_weight_data(treap_node *node) const {
        size_type left = node->m_pos;
        size_type rank_param = m_tree.close(m_tree.enclose(node->m_node));
        size_type right = m_tree.preorder_rank(rank_param) - 2;
        size_type max_pos = m_rmq(left, right);
        node->m_max_weight = m_weights[max_pos];
    }

    treap_node get_data(size_type node, size_type y_old) const {
        // weights
        size_type pos = m_tree.preorder_rank(node) - 2;
        size_type left = pos;
        size_type rank_param = m_tree.close(m_tree.enclose(node));
        size_type right = m_tree.preorder_rank(rank_param) - 2;
        size_type max_pos = m_rmq(left, right);
        size_type w_value = m_weights[pos];
        size_type w_max_value = m_weights[max_pos];
        // coordinates
        size_type y_value = y_old + m_y_values[pos];
        size_type bp_close = m_tree.close(node);
        size_type x_value = bp_close - m_tree.preorder_rank(bp_close);

        return {x_value, y_value, y_value, node, pos, w_value, w_max_value};
    }

    treap_node move_left(treap_node &node) const {
        if (node.m_left != nullptr) {
            return *node.m_left;
        }
        size_type left_node = m_tree.first_child(node.m_node);
        if (left_node != std::numeric_limits<size_t>::max()) {
            return get_data(left_node, node.m_y_value);
        } else {
            size_type max = std::numeric_limits<size_t>::max();
            return {max, max, max};
        }
    }

    treap_node move_right(treap_node &node) const {
        if (node.m_right != nullptr) {
            return *node.m_right;
        }
        size_type right_node = m_tree.next_sibling(node.m_node);
        if (right_node != std::numeric_limits<size_t>::max()) {
            return get_data(right_node, node.m_y_value);
        } else {
            size_type max = std::numeric_limits<size_t>::max();
            return {max, max, max};
        }
    }

    size_type getLeft(size_type v) const {
        return m_tree.first_child(v);
    }

    size_type getRight(size_type v) const {
        return m_tree.next_sibling(v);
    }

    // TODO: auxiliary functions, we're going to delete these
    void print_inorder(std::ostream &os) {
        m_root = &m_root_data;
        if (m_root != nullptr) {
            _print_inorder(m_root, os);
        }
    }

    void print_preorder(std::ostream &os) {
        if (m_root != nullptr) {
            _print_preorder(std::move(m_root), os);
        }
    }

    // TODO: change this to print bp representation, not the pointer representation
    friend std::ostream &operator<<(std::ostream &os, treap_grid &tg) {
        tg.print_preorder(os);
        return os;
    }

    top_k_iterator topk(point_type p1, point_type p2) const {
        return top_k_iterator(*this, p1, p2);
    }

    //! Serializes the data structure into the given ostream
    size_type serialize(std::ostream &out, structure_tree_node *v = nullptr,
            std::string name = "") const {
        structure_tree_node *child = structure_tree::add_child(
                v, name, util::class_name(*this));
        size_type written_bytes = 0;
        written_bytes += write_member(m_size, out, child, "n");
        written_bytes += m_tree.serialize(out, child, "bp_tree");
        written_bytes += m_weights.serialize(out, child, "weights");
        written_bytes += m_y_values.serialize(out, child, "y_values");
        written_bytes += m_rmq.serialize(out, child, "rmq");
        structure_tree::add_size(child, written_bytes);
        return written_bytes;
    }

    //! Loads the data structure from the given istream.
    void load(std::istream &in) {
        read_member(m_size, in);
        m_tree.load(in);
        m_weights.load(in);
        m_y_values.load(in);
        m_rmq.load(in);

        std::stack<std::pair<treap_node *, size_type> > st;
        m_root_data = get_data(1, 0);

        size_t levels = 0;
        st.push({&m_root_data, levels});
        while (!st.empty()) {
            std::pair<treap_node *, size_type> data = st.top();
            st.pop();
            size_type depth = std::get<1>(data);
            treap_node *node = std::get<0>(data);

            if (depth >= t_levels) {
                continue;
            } else {
                treap_node next_node = move_left(*node);
                if (next_node.m_x_value != std::numeric_limits<size_t>::max()) {

                    node->m_left = new treap_node(next_node.m_x_value, next_node.m_y_value, next_node.m_y_value,
                            next_node.m_node, next_node.m_pos, next_node.m_weight, next_node.m_max_weight);
                    st.push({node->m_left, depth + 1});
                }
                next_node = move_right(*node);

                if (next_node.m_x_value != std::numeric_limits<size_t>::max()) {
                    node->m_right = new treap_node(next_node.m_x_value, next_node.m_y_value, next_node.m_y_value,
                            next_node.m_node, next_node.m_pos, next_node.m_weight, next_node.m_max_weight);
                    st.push({node->m_right, depth + 1});
                }
            }
        }
    }

private:
    void _delete_pointers(treap_node *node) {
        if (node != nullptr) {
            _delete_pointers(node->m_left);
            _delete_pointers(node->m_right);
            delete node;
        }
    }

    void _convert_bp(treap_node *node, bit_vector &bitmap, size_type &count) {
        bitmap[count++] = true;
        while (node != nullptr) {
            node->m_node = count;
            _convert_bp(node->m_left, bitmap, count);
            node = node->m_right;
        }
        bitmap[count++] = false;
    }

    size_type find_min(size_type start, size_type end, wt_int<> &wt_y, size_t symbol = 0) {
        size_t mid_range = (start + end) / 2;
        size_t distance = 0;
        size_t min_distance = start - end;
        size_t min_pos = 0;
        size_t pos = 0;
        while (symbol < wt_y.sigma) {
            size_t left_rank = wt_y.rank(start, symbol);
            size_t right_rank = wt_y.rank(end + 1, symbol);
            size_t sym_count = right_rank - left_rank;
            if (sym_count == 0) {
                symbol++;
                continue;
            }
            size_t select_start = left_rank + 1;
            while (select_start != right_rank + 1) {
                pos = wt_y.select(select_start, symbol);
                distance = size_type_abs(pos, mid_range);
                if (distance > min_distance)
                    break;
                if (distance < min_distance) {
                    min_distance = distance;
                    min_pos = pos;
                }
                select_start++;
            }
            break;
        }
        return min_pos;
    }

    void _print_inorder(node_ptr &node, std::ostream &os) {
        if (node != nullptr) {
            _print_inorder(node->m_left, os);
            os << "(" << node->m_x_value << "," << node->m_y_value << ") ";
            _print_inorder(node->m_right, os);
        }
    }
//
//        void _print_preorder(node_ptr node, std::ostream& os) const {
//            if (node != nullptr) {
//                os << "(" << node->m_x_value << "," << node->m_y_value << ") ";
//                _print_preorder(std::move(node->m_left), os);
//                _print_preorder(std::move(node->m_right), os);
//            }
//        }

    size_type size_type_abs(size_type a, size_type b) {
        return a < b ? b - a : a - b;
    }
};
}
#endif