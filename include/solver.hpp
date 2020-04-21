#pragma once
#include "stdafx.hpp"

#include <cstring>
#include <cstdlib>

#include "vmo.hpp"
#include "toposBuild.hpp"
#include "Sparse.hpp"
#include "IO.hpp"


int solveFromTopoNN(VariableSizeMeshContainer<int>& topoNN, const char* type_matr, int threads = 4);
