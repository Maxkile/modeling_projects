#pragma once
#include "stdafx.hpp"

#include <cstring>
#include <cstdlib>
#include <omp.h>

#include "vmo.hpp"
#include "toposBuild.hpp"
#include "Sparse.hpp"
#include "IO.hpp"


int solveFromTopoNN(VariableSizeMeshContainer<int>& topoNN, const char* type_matr, int threads = 4) {
    
    std::vector<double> x;
	std::ofstream fout;
    size_t n;

    Sparse<double>* matrix = nullptr;

    if (threads < 1)
    {
        std::cout << "Threads number can't be lower than 1!" << std::endl;
        return 1;
    }
    else if (!strcmp("-csr", type_matr))
    {
        matrix = new SparseCSR<double>(topoNN);
        std::cout << std::endl << "Type matrix: CSR\n" << endl;
    }
    else if (!strcmp("-ellp",type_matr))
    {
        matrix = new SparseELL<double>(topoNN);
        std::cout << std::endl << "Type matrix: ELLPACK\n" << endl;
    }
    else if (!strcmp("-coo",type_matr))
    {
        matrix = new SparseCOO<double>(topoNN);
        std::cout << std::endl << "Type matrix: COORD\n" << endl;
    }
    else
    {
        std::cout << "Wrong sparse matrix format chosen!" << std::endl;
        return 1;
    }

    n = matrix->getDenseRows();

    vector<double> b(n);

    for (size_t i = 1; i < n + 1; i++)
        b[i - 1] = i;

    x = vmo::conGradSolver(*matrix,b,threads);

    fout.open("decision.txt");
    std::cout << "Writing decision to file..." << std::endl;

    for (size_t i = 0; i < x.size(); i++)
    {
        fout << x[i] << endl;
    }
    std::cout << "Writing to file completed" << std::endl;

    fout.close();

    delete matrix;
	
	return 0;
}
