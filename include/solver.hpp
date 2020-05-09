#pragma once
#include "stdafx.hpp"

#include <cstdlib>
#include <cstring>

#include "IO.hpp"
#include "Sparse.hpp"
#include "toposBuild.hpp"
#include "vmo.hpp"

int solveFromTopoNN(VariableSizeMeshContainer<int> &topoNN,
                    const char *type_matr, int threads = 4);
