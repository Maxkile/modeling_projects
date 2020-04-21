#pragma once

#include "FixedSizeMeshContainer.hpp"
#include "stdafx.hpp"
#include <map>

namespace decomp
{
    /*
     * Mesh decomposing(only decart meshes are supported yet),returns pair - submesh id -> submesh coords
    */
    pair<size_t,vector<int>> decomposeMesh(int Nx, int Ny, int Px, int Py, int px, int py);

    /*
     *  Mesh decomposing(only decart meshes are supported yet). Get all "ibeg, iend; jbeg, jend" for all indexes [0,Px - 1],[0, Py - 1]
    */
    vector<pair<size_t,vector<int>>> decomposeMesh(int Nx, int Ny, int Px, int Py);
    /*
     * Get submmesh global id among submeshes
    */
    size_t getSubmeshIdByCoords(int x,int y, const vector<pair<size_t,vector<int>>>& submeshes, int Nx, int Ny);

    /*
     * Is node of 'id' submesh is interface node of 'current_id' submesh
    */
    bool isInterface(int x, int y, const vector<pair<size_t,vector<int>>>& submeshes, int Nx, int Ny, size_t current_id);

    /*
     * Is node of 'id' submesh is halo node of 'current_id' submesh
    */
    bool isHalo(int x, int y, const vector<pair<size_t,vector<int>>>& submeshes, int Nx, int Ny, size_t current_id);

    /*
     * Global mesh nodes indexing
    */
    void getGlobalIndexes(map<int,int>& G2L, vector<int>& global);

    /*
     * Local mesh nodes indexing
    */
    void getLocalIndexes(map<int,int>& G2L, vector<int>& local);
}
