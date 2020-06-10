#ifndef SOLVER_HPP
#define SOLVER_HPP

#pragma once
#include "stdafx.hpp"

#include <cstdlib>
#include <cstring>

#include "IO.hpp"
#include "Sparse.hpp"
#include "toposBuild.hpp"
#include "vmo.hpp"

namespace solver {
int solveFromTopoNN(const VariableSizeMeshContainer<int> &topoNN, const char *type_matr, NodesInfo *nodesinfo,
                    map<int, int> &list_of_neighbors, vector<set<int>> &send, vector<set<int>> &recv, vector<int> &L2G,
                    int n_own, size_t totalSize, size_t threadsNumber = 4);
}

#endif
