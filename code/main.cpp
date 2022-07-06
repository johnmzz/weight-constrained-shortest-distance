#include <cstdlib>
#include <cstring>
#include <iostream>
#include <sstream>
#include <fstream>
#include <random>

#include "WGraph.h"

using namespace std;

int main(int argc, char* argv[]) {
    string data_path = "/import/vldb01/2/scratch/mazhuo/data/";
    string input_graph = argv[1];
    string input_format = argv[2];
    string graph_type = argv[3];
    string type = argv[4];
    int threshold = stoi(argv[5]);

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

    WGraph g(data_path, input_graph, isBin);

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

    g.build_index(type, threshold);
    g.get_index_size(type);

    cout << "Complete!" << endl;
    return 0;
}