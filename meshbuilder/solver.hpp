#pragma once
#include "stdafx.hpp"

#include <cstring>
#include <cstdlib>
#include <omp.h>

#include "vmo.hpp"
#include "toposBuild.hpp"
#include "Sparse.hpp"
#include "IO.hpp"

int solver(VariableSizeMeshContainer<int>& topoNN, const char* type_matr, int threads = 4) {
    double start, end;
    std::vector<double> x;
	std::ofstream fout;
    size_t n;

    if (threads < 1)
    {
        std::cout << "Threads number can't be lower than 1!" << std::endl;
        return 1;
    }

    if (!strcmp("-csr", type_matr))
    {
        //csr
        SparseCSR<double> matrix(topoNN);
        n = matrix.getDenseRows();
        vector<double> b(n);

		std::cout << std::endl << "Type matrix: CSR\n" << endl;

        for (size_t i = 1; i < n + 1; i++)
            b[i - 1] = i;

        start = omp_get_wtime();
        x = vmo::conGradSolver(matrix,b,threads);
        end = omp_get_wtime();

        std::cout << "\nTime solution:\n\t" << end - start << " sec" << std::endl;
		
    }
    else if (!strcmp("-ellp",type_matr))
    {
        //ellpack
		SparseELL<double> matrix(topoNN);

        n = matrix.getDenseRows();
        std::vector<double> b(n);

		std::cout << std::endl << "Type matrix: ELLPACK\n" << std::endl;

        for (size_t i = 1; i < n + 1; i++)
            b[i - 1] = i;

        start = omp_get_wtime();
        x = vmo::conGradSolver(matrix,b,threads);
        end = omp_get_wtime();

        std::cout << "\nTime solution:\n\t" << end - start << " sec" << std::endl;

    }
    else if (!strcmp("-coo",type_matr))
    {
        //coord
        SparseCOO<double> matrix(topoNN);

        n = matrix.getDenseRows();
        std::vector<double> b(n);

		std::cout << std::endl << "Type matrix: COO\n" << std::endl;
		        
        for (size_t i = 1; i < n + 1; i++)
            b[i - 1] = i;

        start = omp_get_wtime();
        x = vmo::conGradSolver(matrix,b,threads);
        end = omp_get_wtime();

        std::cout << "\nTime solution:\n\t" << end - start << " sec" << std::endl;
    }
    else
    {
        std::cout << "Wrong sparse matrix format chosen!" << std::endl;
        return 1;
    }

    fout.open("decision.txt");
    std::cout << "Writing decision to file..." << std::endl;
    
    for (size_t i = 0; i < x.size(); i++)
    {
        fout << x[i] << endl;
    }
    std::cout << "Completed" << std::endl;

    fout.close();
	
	return 0;
}
