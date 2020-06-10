#ifndef DECOMP_HPP
#define DECOMP_HPP

#pragma once

#include <map>
#include <set>

#include "FixedSizeMeshContainer.hpp"
#include "stdafx.hpp"
#include "vmo.hpp"

namespace decomp {

/*
 * Mesh decomposing(only decart meshes are supported yet),returns pair - submesh
 * id -> submesh coords
 */
pair<size_t, vector<int>> decomposeMesh(int Nx, int Ny, int Px, int Py, int px, int py);

/*
 *  Mesh decomposing(only decart meshes are supported yet). Get all "ibeg, iend;
 * jbeg, jend" for all indexes [0,Px - 1],[0, Py - 1]
 */
vector<pair<size_t, vector<int>>> decomposeMesh(int Nx, int Ny, int Px, int Py);

/*
 * Get submmesh global id among submeshes
 */
size_t getSubmeshIdByCoords(int x, int y, const vector<pair<size_t, vector<int>>> &submeshes, int Nx, int Ny);

/*
 * Is node of 'id' submesh is interface node of 'current_id' submesh
 */
bool isInterface(int x, int y, const vector<pair<size_t, vector<int>>> &submeshes, const size_t target_submesh_id,
                 int Nx, int Ny);

/*
 * Is node of 'id' submesh is halo node of 'current_id' submesh
 */
bool isHalo(int x, int y, const vector<pair<size_t, vector<int>>> &submeshes, const size_t target_submesh_id, int Nx,
            int Ny);

/*
 * Finds for vector with 'node_id' as a key and inserts node. If no such key was found, creates it and expands
 * vector.
 */
void insertHalo(vector<pair<size_t, vector<int>>> &haloes, int node_pos, size_t node_id);

/*
 * Global mesh nodes indexing
 */
void getGlobalIndexes(map<int, int> &G2L, vector<int> &global);

/*
 * Local mesh nodes indexing
 */
void getLocalIndexes(map<int, int> &G2L, vector<int> &local);
} // namespace decomp

#endif
