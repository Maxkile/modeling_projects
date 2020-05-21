#ifndef PARALLEL_HPP
#define PARALLEL_HPP

#pragma once

#include "VariableSizeMeshContainer.hpp"
#include "stdafx.hpp"
#include "toposBuild.hpp"
#include <algorithm>
#include <map>
#include <set>

using namespace std;

namespace parallel {
void build_list_of_neighbors(map<int, int> &list_of_neighbors, const vector<int> &part, int self_id);

VariableSizeMeshContainer<int> build_topoNN_2(VariableSizeMeshContainer<int> &topoNN, const vector<int> &L2G,
                                              unsigned n_own);

void build_list_send_recv(VariableSizeMeshContainer<int> &topoNN_2, map<int, int> &G2L, vector<int> &L2G,
                          vector<int> &part, map<int, int> &list_of_neighbors, vector<set<int>> &send,
                          vector<set<int>> &recv, unsigned n_own, unsigned self_id);
} // namespace parallel

#endif // PARALLEL_HPP
