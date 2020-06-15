#ifndef SOLVER_HPP
#define SOLVER_HPP

#pragma once
#include "stdafx.hpp"

#include <cstdlib>
#include <cstring>

#include "IO.hpp"
#include "Sparse.hpp"
#include "solverInfo.hpp"
#include "toposBuild.hpp"
#include "vmo.hpp"

namespace solver {
template <typename T> void defaultInitB(vector<double> &b, const vector<T> &source);

int solveFromTopoNN(const VariableSizeMeshContainer<int> &topoNN, NodesInfo *nodesinfo,
                    map<int, int> &list_of_neighbors, vector<set<int>> &send, vector<set<int>> &recv, vector<int> &L2G,
                    int n_own, size_t totalSize, double &timeSpent, SolverInfo *solverInfo);
} // namespace solver

#endif
