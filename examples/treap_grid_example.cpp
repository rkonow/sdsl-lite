#include "sdsl/bp_tree.hpp"
#include "sdsl/treap_grid.hpp"
#include "sdsl/treap_grid_algorithm.hpp"
using namespace sdsl;
using namespace std;

void getValues(treap_grid<>& tg, size_t node) {
    size_t pos = tg.getTree().preorder_rank(node)-2;
    cout << " weight = " << tg.m_weights[pos] << " ";
    cout << " y-dest = " << tg.m_y_values[pos] << " " ;
}

void getMaxInRange(treap_grid<>& tg, pair<size_t,size_t> rng) {
    size_t max_position = tg.m_rmq(rng.first,rng.second);
    cout << " [ max position: " <<  max_position << "] ";
    cout << " [ max value:"    << tg.m_weights[max_position] <<"] ";
}

std::pair<size_t,size_t> getXRange(bp_tree<>& bp, size_t node) {
    size_t left = bp.preorder_rank(node)-2;
    size_t rank_param = bp.close(bp.enclose(node));
    size_t right = bp.preorder_rank(rank_param)-2;
    cout << " [" << left << " , " << right << "]";
    return {left,right};
}

void getX(bp_tree<>& bp, size_t node) {
    size_t bp_close =  bp.close(node);
    cout <<  "x = " << bp_close - bp.preorder_rank(bp_close);
}

void traverse(treap_grid<>& tg, size_t node) {
    if (node != std::numeric_limits<size_t>::max()) {
        bp_tree<> bp = tg.getTree();
        getX(bp,node);
        getValues(tg,node);
        pair<size_t,size_t> values = getXRange(bp,node);
        getMaxInRange(tg,values);
        cout << endl;
        traverse(tg, tg.getLeft(node));
        traverse(tg, tg.getRight(node));
    }
}

int main(int argc, char** argv) {

    {
        int_vector<> x = {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
        int_vector<> y = {0, 0, 0, 1, 2, 3, 4, 0, 0, 0, 1};
        int_vector<> w = {1, 1, 1, 1, 1, 1, 1, 10, 1, 1, 9};

        int_vector_buffer<> buf_w("test.w", std::ios::out);
        int_vector_buffer<> buf_y("test.y", std::ios::out);
        int_vector_buffer<> buf_x("test.x", std::ios::out);

        for (size_t i = 0; i < x.size(); i++) {
            buf_w.push_back(w[i]);
            buf_y.push_back(y[i]);
            buf_x.push_back(x[i]);
        }
    }
//    {
//        int_vector_buffer<> buf_w("test.w");
//        int_vector_buffer<> buf_y("test.y");
//        int_vector_buffer<> buf_x("test.x");
//        treap_grid<> tg(buf_x, buf_y, buf_w);
//        cout << "pointer nodes: ";
//        cout << tg << endl;
//        traverse(tg,1);
//        auto mts_it = map_to_sorted_sequence(tg, {0, 10}, {0,4});
//        while (mts_it) {
//            cout << get<0>(*mts_it) << " , "  << get<1>(*mts_it)  << " , " << get<2>(*mts_it) << " , " << get<3>(*mts_it) << endl;
//            ++mts_it;
//        }
//    }
    {
        treap_grid<> tg2;
        construct_im(tg2, {{0,0,2},{1,2,3},{2,1,2},{3,0,2},{4,0,1},{5,1,4},{6,0,1},{7,1,1},{8,0,8},{9,2,5}});
        traverse(tg2, 1);
        auto mts_it = map_to_sorted_sequence(tg2, {2, 7}, {0,1});
        while (mts_it) {
            cout << get<0>(*mts_it) << " , "  << get<1>(*mts_it)  << " , " << get<2>(*mts_it) << " , " << get<3>(*mts_it) << endl;
            ++mts_it;
        }
        auto topk_it2 = top_k(tg2, {2,0}, {7,1});
        while (topk_it2) {
            auto point_weight = *topk_it2;
            cout << point_weight.first <<" weight: "<<point_weight.second << endl;
            ++topk_it2;
        }
    }
}



//}


