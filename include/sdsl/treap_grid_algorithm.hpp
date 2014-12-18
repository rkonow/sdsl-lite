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
        typedef typename t_treap_grid::size_type value_type;
        typedef typename t_treap_grid::size_type  size_type;
        typedef std::pair<int_vector<>::size_type, int_vector<>::size_type> range_type;
        // x,y,node,pos
        typedef std::tuple<size_type, size_type, size_type, size_type> t_data;
        // x_value, y_value, w_value, w_max //
        typedef std::array<size_type,4> ret_type;

private:
        const t_treap_grid*     m_treap_grid = nullptr;
        range_type              m_y_range;
        range_type              m_x_range;
        std::stack<t_data>    m_stack;
        ret_type                 m_ret;
        bool                    m_valid = false;

        void cond_push(t_data& data) {
            if (std::get<0>(data) != std::numeric_limits<size_type>::max()) {
                // fix this
                size_type x = std::get<0>(data);
                size_type y = std::get<1>(data);
                size_type n = std::get<2>(data);
                size_type p = std::get<3>(data);
                m_stack.emplace(x, y, n, p);
            }
        }

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
                // fix this...
                size_type x = std::get<0>(tg->m_root_data);
                size_type y = std::get<1>(tg->m_root_data);
                size_type n = std::get<2>(tg->m_root_data);
                size_type p = std::get<3>(tg->m_root_data);

                m_stack.emplace(x,y,n,p);
                ++(*this);
            }
        }

        treap_grid_iterator& operator++()
        {
            m_valid = false;
            while (!m_stack.empty()) {
                auto x_value = std::get<0>(m_stack.top());
                auto y_value = std::get<1>(m_stack.top());
                auto node    = std::get<2>(m_stack.top());
                auto pos     = std::get<3>(m_stack.top());

                m_stack.pop();
                if (std::get<0>(m_x_range) <= x_value and std::get<1>(m_x_range) >= x_value
                        and std::get<0>(m_y_range) <= y_value and std::get<1>(m_y_range) >= y_value)
                {
                    std::pair<size_type, size_type> weight_data = m_treap_grid->get_weight_data(node,pos);
                    m_ret = {x_value,y_value, weight_data.first, weight_data.second};
                    t_data data = m_treap_grid->move_left(node, y_value);
                    cond_push(data);
                    data = m_treap_grid->move_right(node, y_value);
                    cond_push(data);
                    m_valid  = true;
                    break;
                } else {
                    t_data data;
                    if (std::get<0>(m_x_range) > x_value) {
                        data = m_treap_grid->move_left(node,y_value);
                    } else if(std::get<0>(m_x_range) < x_value) {
                        data = m_treap_grid->move_right(node,y_value);
                    }
                    cond_push(data);
                }
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