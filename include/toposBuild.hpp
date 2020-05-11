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
constexpr int TRIANGLE_NODES = 3;
constexpr int SQUARE_NODES = 4;

// coords
void build_coord(FixedSizeMeshContainer<double> &C, int Lx, int Ly, int Nx,
                 int Ny);

// defines how many triangles and squares left in mesh after skipping elements
int computeMeshFiguresNumberLeft(int figCount1, int figCount2,
                                 int skippedElemsCount,
                                 int curMeshFigureStructure);

VariableSizeMeshContainer<int>
toLocalIndexes(const VariableSizeMeshContainer<int> &originEN,
               map<int, int> &G2L);

vector<int> toLocalIndexes(const vector<int> &origin, map<int, int> &G2L);

// topoEN
// Also generates G2L,L2G,halo's,interfaces vectors required for topologies
// generating
VariableSizeMeshContainer<int>
build_topoEN(int Nx, int Ny, int k3, int k4, size_t submesh_id,
             const vector<pair<size_t, vector<int>>> &submeshes,
             map<int, int> &G2L, vector<int> &L2G, vector<int> &inner,
             vector<int> &interface, vector<int> &haloes);

// topoSN
VariableSizeMeshContainer<int> build_topoSN(int Nx, int Ny, int k3, int k4);

// reversed topology
VariableSizeMeshContainer<int>
build_reverse_topo(const VariableSizeMeshContainer<int> &topo);

// topoBNS
VariableSizeMeshContainer<int>
build_topoBNS(const VariableSizeMeshContainer<int> &topo);

// topoBSN
VariableSizeMeshContainer<int> build_topoBSN(int Nx, int Ny);

// topoNN_1
VariableSizeMeshContainer<int>
build_topoNN_1(const VariableSizeMeshContainer<int> &topoSN);

// topoNN_2
VariableSizeMeshContainer<int>
build_topoNN_2(const VariableSizeMeshContainer<int> &topoEN);
  
} // namespace topos
