#include <cstdlib>
#include <cstring>
#include <iostream>
#include <sstream>
#include <fstream>
#include <random>

#include "WGraph.h"

using namespace std;

int main(int argc, char* argv[]) {
    string input_graph = argv[1];
    string input_format = argv[2];
    string graph_type = argv[3];
    string type = argv[4];

    bool isBin;
    if (input_format == "bin") {
        isBin = true;
    }
    else if (input_format == "txt") {
        isBin = false;
    }
    else {
        cout << "Incorrect input format" << endl;
        exit(0);
    }

    WGraph g(input_graph, isBin);

    if (graph_type == "road") {
        g.tree_decomp();
    }
    else if (graph_type == "social") {
        g.tree_decomp_partial(0.7);
    }
    else {
        cout << "Don't perform tree_decomp" << endl;
    }

    g.get_graph_statistics();

    g.build_index(type, graph_type);
    g.get_index_size(type);

    cout << "Complete!" << endl;
    return 0;
}