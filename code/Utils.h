#pragma once
#ifndef _Utils_H
#define _Utils_H

#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include <iterator>
#include <algorithm>
#include <set>
#include <cmath>
#include <string>
#include <tuple>
#include <boost/dynamic_bitset.hpp>
#include <queue>
#include <limits.h>
#include <chrono>

using namespace std;

typedef vector<vector<uint32_t> > AdjList;
typedef vector<vector<uint8_t>> WeightList;
// typedef vector<uint32_t> WeightList;
#define INVALID_VALUE 0xffffffff
#define INVALID_VALUE_8bit 0xff
#define INVALID_VALUE_16bit 0xffff
#define INF 0xffffffff

struct descWeightSort
{
    bool operator()(pair<uint32_t, uint32_t> const& a, pair<uint32_t, uint32_t> const& b) {
        if (a.second != b.second) {
            return (a.second < b.second);
        }
        else {
            return (a.first > b.first);
        }
    }
};

struct IndexStruct
{
    vector<vector<uint32_t> > labeling_v;
    vector<vector<uint32_t> > offset;
    vector<vector<uint32_t> > labeling_d;
    vector<vector<uint32_t> > labeling_w;
};

inline uint32_t query_with_index_dist(IndexStruct& myidx, uint32_t s, uint32_t t, uint32_t r, double& qtime) {
    auto t1 = std::chrono::high_resolution_clock::now();
    uint32_t min_dist = INVALID_VALUE;
    uint32_t i = 0, j = 0;

    while (i < myidx.labeling_v[s].size() && j < myidx.labeling_v[t].size()) {
        uint32_t curr_u = myidx.labeling_v[s][i];
        uint32_t curr_w = myidx.labeling_v[t][j];
        if (curr_u < curr_w) {
            i++;
        }
        else if (curr_u > curr_w) {
            j++;
        }
        else {
            uint32_t d1 = INVALID_VALUE, d2 = INVALID_VALUE;
            for (uint32_t p = myidx.offset[s][i]; p < myidx.offset[s][i + 1]; p++) {
                if (myidx.labeling_w[s][p] >= r) {
                    d1 = myidx.labeling_d[s][p];
                    break;
                }
            }
            for (uint32_t q = myidx.offset[t][j]; q < myidx.offset[t][j + 1]; q++) {
                if (myidx.labeling_w[t][q] >= r) {
                    d2 = myidx.labeling_d[t][q];
                    break;
                }
            }

            if (d1 != INVALID_VALUE && d2 != INVALID_VALUE) {
                if (d1 + d2 < min_dist) { min_dist = d1 + d2; }
            }
            i++;
            j++;
        }
    }

    while (i < myidx.labeling_v[s].size()) {
        uint32_t curr_u = myidx.labeling_v[s][i];
        if (curr_u == t) {
            uint32_t d1 = INVALID_VALUE;
            for (uint32_t p = myidx.offset[s][i]; p < myidx.offset[s][i + 1]; p++) {
                if (myidx.labeling_w[s][p] >= r) {
                    d1 = myidx.labeling_d[s][p];
                    break;
                }
            }
            if (d1 != INVALID_VALUE && d1 < min_dist) {
                min_dist = d1;
            }
        }
        i++;
    }

    while (j < myidx.labeling_v[t].size()) {
        uint32_t curr_w = myidx.labeling_v[t][j];
        if (curr_w == s) {
            uint32_t d2 = INVALID_VALUE;
            for (uint32_t q = myidx.offset[t][j]; q < myidx.offset[t][j + 1]; q++) {
                if (myidx.labeling_w[t][q] >= r) {
                    d2 = myidx.labeling_d[t][q];
                    break;
                }
            }
            if (d2 != INVALID_VALUE && d2 < min_dist) {
                min_dist = d2;
            }
        }
        j++;
    }
    auto t2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> fp_ms = t2 - t1;
    qtime += fp_ms.count();
    return min_dist;
}
#endif
