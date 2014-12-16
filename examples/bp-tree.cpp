#include "sdsl/bp_tree.hpp"
#include "sdsl/treap_grid.hpp"
#include "../include/sdsl/treap_grid.hpp"

//#include "sdsl/bit_vector.hpp"
//#include "../include/sdsl/bp_tree.hpp"

using namespace sdsl;
using namespace std;

void getXRange(bp_tree<>& bp, size_t node) {
    size_t left = bp.preorder_rank(node)-2;
    size_t bp_close =  bp.close(node);
    size_t right = bp_close -  bp.preorder_rank(bp.close(node));
    cout << "left = " << left << endl;
    cout << "right = " << right << endl;
    cout << "-----" << endl;
}
void getX(bp_tree<>& bp, size_t node) {
    size_t bp_close =  bp.close(node);
    cout <<  bp_close - bp.preorder_rank(bp_close) << endl;

}

void traverse(treap_grid<>& tg, size_t node) {
    if (node != std::numeric_limits<size_t>::max()) {
        bp_tree<> bp = tg.getTree();
        getX(bp,node);
        traverse(tg, tg.getLeft(node));
        traverse(tg, tg.getRight(node));
    }
}

int main(int argc, char** argv) {
//// tree extracted from Navarro and Sadakane's paper
//    string bvs = "((()()())(()()))";
    bit_vector bp(16,0);
//    for (size_t i = 0 ; i < bvs.size(); i++) {
//        if (bvs[i] == '(')
//            bp[i] = 1;
//        else
//            bp[i] = 0;
//    }
//
//    bp_tree<> bpTree(&bp);
//    cout << "nodes = \t" << bpTree.nodes() << endl;
//    cout << "root = \t" << bpTree.root() << endl;
//    cout << "first child = \t" << bpTree.first_child(0) << endl;
//    cout << "subtree_size = \t" << bpTree.subtree_size(0) << endl;
//    cout << "subtree_size = \t" << bpTree.subtree_size(1) << endl;
//    cout << "subtree_size = \t" << bpTree.subtree_size(2) << endl;
//    cout << "next_sibling of first child = \t" << bpTree.next_sibling(1) << endl;
//    cout << "subtree_size = \t" << bpTree.subtree_size (9) << endl;

    int_vector<> x = {0,1,2,3,4,5,6,7,8,9,10};
    int_vector<> y = {0,0,0,1,2,3,4,0,0,0,1};
    int_vector<> w = {1,1,1,1,1,1,1,10,1,1,9};

    int_vector_buffer<> buf_w("test.w", std::ios::out);
    int_vector_buffer<> buf_y("test.y", std::ios::out);
    int_vector_buffer<> buf_x("test.x", std::ios::out);

    for (size_t i = 0 ; i < x.size(); i++) {
        buf_w.push_back(w[i]);
        buf_y.push_back(y[i]);
        buf_x.push_back(x[i]);
    }
    buf_w.close(false);
    buf_y.close(false);
    buf_x.close(false);
    treap_grid<> tg(buf_x,buf_y,buf_w);


    bp_tree<> bpTree = tg.getTree();
    tg.print_preorder();

    cout << "nodes = \t" << bpTree.nodes() << endl;
    cout << "root = \t" << bpTree.root() << endl;
    cout << "first child = \t" << bpTree.first_child(0) << endl;
    cout << "subtree_size = \t" << bpTree.subtree_size(0) << endl;
    cout << "subtree_size = \t" << bpTree.subtree_size(2) << endl;
    cout << "next_sibling of first child = \t" << bpTree.next_sibling(1) << endl;
    cout << "next_sible of next_sibling of first child = \t" << bpTree.next_sibling(bpTree.next_sibling(1)) << endl;
    cout << "first child = \t" << bpTree.first_child(1) << endl;
    cout << "first child = \t" << bpTree.first_child(2) << endl;


    getXRange(bpTree,1);
    getX(bpTree,1);
    getXRange(bpTree,17);
    getX(bpTree,17);
    getXRange(bpTree,21);
    getX(bpTree,21);
    cout << "======" << endl;
    getX(bpTree,2);
    cout << "traversal..." << endl;
    traverse(tg,1);




//    tg.print_inorder();
//    tg.getX(1);

}