Weight Constrained Shortest Path
=============================================

Implementation of WCSD..

Please cite the paper as follows:

Compilation and Usage
---------------------

We tested our program on Debian 4.19 amd64, and it requires the following
packages: `cmake`,`boost`:

Compilation:
```
git clone <this repository>
cd <repository_name>
make
```

Usage:
```
$ wcsd [path to input graph] [input file type] [graph type] [version]
```

**Example.** Process road network "NY" graph in "bin" format with "wc-index-plus":
```
$ ./wcsd NY bin road wc-index-plus
```

Functions
---------------------
Index construction:
* wc-index: vertex_prioritized_indexing();
* wc-index-plus: vertex_prioritized_indexing_plus();
* naive-index: build_naive_index();

Online solutions:
* constrained_shortest_distance_naive();
* constrained_shortest_distance_dijkstra();
* constrained_shortest_distance_plus();
