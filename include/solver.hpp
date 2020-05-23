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
                    size_t threadsNumber = 4);
}

#endif
