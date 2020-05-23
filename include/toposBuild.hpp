#ifndef TOPOS_HPP
#define TOPOS_HPP

#pragma once
#include "stdafx.hpp"

#include <iostream>
#include <map>
#include <set>

#include "FixedSizeMeshContainer.hpp"
#include "VariableSizeMeshContainer.hpp"
#include "decomposition.hpp"
#include "vmo.hpp"

namespace topos {

using mapping = std::map<int, int>;

constexpr int TRIANGLE_NODES = 3;
constexpr int SQUARE_NODES = 4;

// coords
void build_coord(FixedSizeMeshContainer<double> &C, int Lx, int Ly, int Nx, int Ny);

// defines how many triangles and squares left in mesh after skipping elements
int computeMeshFiguresNumberLeft(int figCount1, int figCount2, int skippedElemsCount, int curMeshFigureStructure);

VariableSizeMeshContainer<int> toLocalIndexes(const VariableSizeMeshContainer<int> &originEN, const mapping &G2L);

VariableSizeMeshContainer<int> toGlobalIndexes(const VariableSizeMeshContainer<int> &topoNN, const vector<int> &L2G,
                                               size_t n_own);

vector<int> toLocalIndexes(const vector<int> &origin, const mapping &G2L);

vector<int> toGlobalIndexes(const vector<int> &origin, const vector<int> &L2G);

// topoEN
// Also generates G2L,L2G,halo's,interfaces vectors required for topologies
// generating
VariableSizeMeshContainer<int> build_topoEN(int Nx, int Ny, int k3, int k4, size_t submesh_id,
                                            const vector<pair<size_t, vector<int>>> &submeshes, mapping &G2L,
                                            vector<int> &L2G, vector<int> &nodes, vector<int> &part, size_t &n_own);

// topoSN
VariableSizeMeshContainer<int> build_topoSN(int Nx, int Ny, int k3, int k4);

// reversed topologyconst
VariableSizeMeshContainer<int> build_reverse_topo(const VariableSizeMeshContainer<int> &topo);

// topoBNS
VariableSizeMeshContainer<int> build_topoBNS(const VariableSizeMeshContainer<int> &topo);

// topoBSN
VariableSizeMeshContainer<int> build_topoBSN(int Nx, int Ny);

// topoNN_1
VariableSizeMeshContainer<int> build_topoNN_from_topoSN(const VariableSizeMeshContainer<int> &topoSN);

// topoNN_2
VariableSizeMeshContainer<int> build_topoNN_from_topoEN(const VariableSizeMeshContainer<int> &topoEN);

} // namespace topos

#endif
