/* sdsl - succinct data structures library
    Copyright (C) 2014 Simon Gog

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
/*! \file treap_grid_algorithm.hpp
    \brief treap_grid_algorithm.hpp contains an implementation of several algorithms for treap_grid
    \author Simon Gog, Roberto Konow
*/

#ifndef INCLUDED_SDSL_TREAP_GRID_ALGORITHM
#define INCLUDED_SDSL_TREAP_GRID_ALGORITHM

#include <limits>
#include <utility>
#include <vector>
#include <stack>
#include <queue>
#include "sdsl/int_vector_buffer.hpp"
#include "sdsl/int_vector.hpp"


namespace sdsl {

// forward declaration
template<size_t levels,
        typename t_bp,
        typename t_rmq,
        typename t_weight_vec,
        typename t_y_vec
>
class treap_grid;

template<size_t t_levels,
        typename t_bp,
        typename t_rmq,
        typename t_weight_vec,
        typename t_y_vec>

typename treap_grid<t_levels, t_bp, t_rmq, t_weight_vec, t_y_vec>::top_k_iterator
top_k(const treap_grid<t_levels, t_bp, t_rmq, t_weight_vec, t_y_vec> &tg,
        typename treap_grid<t_levels, t_bp, t_rmq, t_weight_vec, t_y_vec>::point_type p1,
        typename treap_grid<t_levels, t_bp, t_rmq, t_weight_vec, t_y_vec>::point_type p2) {
    return tg.topk(p1, p2);
}

template<size_t t_levels,
        typename t_bp,
        typename t_rmq,
        typename t_weight_vec,
        typename t_y_vec>
void
construct(treap_grid<t_levels, t_bp, t_rmq, t_weight_vec, t_y_vec> &idx, std::string file) {
    int_vector_buffer<> buf_x(file + ".x", std::ios::in);
    int_vector_buffer<> buf_y(file + ".y", std::ios::in);
    int_vector_buffer<> buf_w(file + ".w", std::ios::in);
    treap_grid<t_levels, t_bp, t_rmq, t_weight_vec, t_y_vec> tmp(buf_x, buf_y, buf_w);
    tmp.swap(idx);
}

//! Specialized version of method ,,construct_im'' for wt_topk
template<size_t t_levels,
        typename t_bp,
        typename t_rmq,
        typename t_weight_vec,
        typename t_y_vec>
void
construct_im(treap_grid<t_levels, t_bp, t_rmq, t_weight_vec, t_y_vec> &idx, std::vector <std::array<uint64_t, 3>> data) {
    std::string y_file = ram_file_name(std::to_string(util::pid()) + "_" + std::to_string(util::id()));
    std::string w_file = ram_file_name(std::to_string(util::pid()) + "_" + std::to_string(util::id()));
    {
        int_vector<> y(data.size());
        int_vector<> w(data.size());
        for (size_t i = 0; i < data.size(); ++i) {
            if (data[i][0] != i) {
                throw std::logic_error("Construct treap_grid: the first elements\
                      of the tuples should form the identity permutation\
                      [0..data.size()-1]");
            }
            y[i] = data[i][1];
            w[i] = data[i][2];
        }
        util::bit_compress(y);
        store_to_file(y, y_file);
        util::bit_compress(w);
        store_to_file(w, w_file);
    }
    int_vector_buffer<> y_buf(y_file);
    int_vector_buffer<> w_buf(w_file);
    treap_grid <t_levels, t_bp, t_rmq, t_weight_vec, t_y_vec> tmp(y_buf, y_buf, w_buf);
    tmp.swap(idx);
}

// TODO: maybe change the name?
template<typename t_treap_grid>
class treap_grid_iterator
{
    public:
        typedef void(*t_mfptr)();
        typedef typename t_treap_grid::size_type                            value_type;
        typedef typename t_treap_grid::size_type                            size_type;
        typedef std::pair<int_vector<>::size_type, int_vector<>::size_type> range_type;

        typedef typename t_treap_grid::treap_node                            treap_node;
        struct comp_cand {
            bool operator()(treap_node& a, treap_node& b ) { return a.m_max_weight < b.m_max_weight; }
        };
        struct comp_res {
            bool operator()(treap_node& a, treap_node& b ) { return  a.m_weight < b.m_weight; }
        };


    // x_value, y_value, w_value, w_max //
         typedef std::array<size_type,4>                                   ret_type;

private:
        const t_treap_grid*                     m_treap_grid = nullptr;
        range_type                              m_y_range;
        range_type                              m_x_range;
        ret_type                                m_ret;
        bool                                    m_valid = false;
        std::stack<treap_node>                  m_stack;
        std::priority_queue< treap_node , std::vector<treap_node>, comp_cand > cand;
        std::priority_queue< treap_node , std::vector<treap_node>, comp_res > res;


    public:
        treap_grid_iterator() = default;
        treap_grid_iterator(const treap_grid_iterator&) = default;
        treap_grid_iterator(treap_grid_iterator&&) = default;
        treap_grid_iterator& operator=(const treap_grid_iterator&) = default;
        treap_grid_iterator& operator=(treap_grid_iterator&&) = default;



    treap_grid_iterator(const t_treap_grid* tg, const range_type& x_range,
                const range_type& y_range) : m_treap_grid(tg),
                                             m_y_range(y_range),
                                             m_x_range(x_range),
                                             m_stack(),
                                             m_ret(),
                                             m_valid(false)
        {
            if (tg != nullptr) {
                cand.push(tg->m_root_data);
                ++(*this);
            }
        }

        treap_grid_iterator& operator++()
        {
            m_valid = false;
            while(!cand.empty()) {
                treap_node node = cand.top();
                cand.pop();
                if (node.m_y_value <= std::get<1>(m_y_range)) {
                    if (node.m_x_value >= std::get<0>(m_x_range) and node.m_x_value <= std::get<1>(m_x_range)) {
                        res.push(node);
                        if (cand.top().m_max_weight <= res.top().m_weight) {
                            treap_node node = res.top();
                            m_ret = {node.m_x_value, node.m_y_value, node.m_weight, node.m_max_weight};
                            res.pop();
                            m_valid = true;
                            return *this;
                        }
                    }
                    if (node.m_x_value >= std::get<0>(m_x_range)) {
                        treap_node left = m_treap_grid->move_left(node);
                        if (left.m_x_value != std::numeric_limits<size_type>::max()) {
                            cand.push(left);
                        }
                    }
                    if (node.m_x_value <= std::get<1>(m_x_range)) {
                        treap_node right = m_treap_grid->move_right(node);
                        if (right.m_x_value != std::numeric_limits<size_type>::max()) {
                            cand.push(right);
                        }
                    }
                }
            }
            if (res.size() > 0) {
                treap_node node = res.top();
                m_ret = {node.m_x_value, node.m_y_value, node.m_weight, node.m_max_weight};
                res.pop();
                m_valid = true;
            }
            return *this;
        }

        //! Postfix increment of the iterator
        treap_grid_iterator operator++(int)
        {
            treap_grid_iterator it = *this;
            ++(*this);
            return it;
        }

        ret_type operator*() const
        {
            return m_ret;
        }

        operator t_mfptr() const
        {
            return (t_mfptr)(m_valid);
        }
};

template<typename t_treap_grid>
treap_grid_iterator<t_treap_grid>
map_to_sorted_sequence(const t_treap_grid& tg, const std::pair<size_t,size_t>& x_range, const std::pair<size_t,size_t>& y_range)
{
    return treap_grid_iterator<t_treap_grid>(&tg, x_range, y_range);
}
}
#endif