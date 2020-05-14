#pragma once
#include "stdafx.hpp"

#include <cstdlib>
#include <cstring>

#include "IO.hpp"
#include "Sparse.hpp"
#include "toposBuild.hpp"
#include "vmo.hpp"

int solveFromTopoNN(const VariableSizeMeshContainer<int> &topoNN, const char *type_matr, NodesInfo *nodesinfo,
                    size_t threadsNumber = 4);
