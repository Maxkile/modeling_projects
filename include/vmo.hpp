#ifndef VMO_HPP
#define VMO_HPP

#pragma once

#include "Sparse.hpp"
#include "VariableSizeMeshContainer.hpp"
#include "nodeType.hpp"
#include "parallel.hpp"
#include "platformDependencies.hpp"
#include "stdafx.hpp"

#include <cmath>
#include <iomanip>
#include <omp.h>

#include <mpi.h>

extern int proc_id;     // process ID
extern int proc_number; // number of processes
extern MPI_Comm MCW;    // default communicator
extern int MASTER_ID;   // master process ID

// vector matrix operations: consider we are having same template arguments
namespace vmo {

#define LINE_SEPARATOR "-------------------------------------------------------------------------------\n"

static double timeDot = 0;
static double timeMultiply = 0;
static double timeAxpby = 0;
static size_t i = 0;

template <typename VT, typename ST>
void axpby(std::vector<VT> &x, const std::vector<VT> &y, const ST &a, const ST &b, NodesInfo *nodesinfo,
           size_t threadsNumber = 4) {
    size_t size = x.size();
    if (size != y.size()) {
        std::cerr << "Incompatible sizes in axpby!" << std::endl;
    } else {
        double start = omp_get_wtime();
#pragma omp parallel for num_threads(threadsNumber)
        for (OPENMP_INDEX_TYPE i = 0; i < size; ++i) {
            x[i] = a * x[i] + b * y[i];
        }
        timeAxpby += omp_get_wtime() - start;
    }
}

template <typename VT>
VT dot(const std::vector<VT> &x, const std::vector<VT> &y, int n_own, NodesInfo *nodesinfo, size_t threadsNumber = 4);

template <typename VT>
void multiply(const std::vector<VT> &x, std::vector<VT> &result, VT scalar, size_t threadsNumber = 4) {
    double start = omp_get_wtime();
#pragma omp parallel for num_threads(threadsNumber)
    for (OPENMP_INDEX_TYPE i = 0; i < x.size(); ++i) {
        result[i] = scalar * x[i];
    }
    timeMultiply += omp_get_wtime() - start;
}

// Compute eucledian norm
double norm(const std::vector<double> &vec, int n_own, NodesInfo *nodesinfo);

template <typename VT> void join(vector<VT> &target, const vector<VT> &arg) {
    size_t size = arg.size();
    target.reserve(size);
    for (size_t i = 0; i < size; ++i) {
        target.push_back(arg[i]);
    }
}

template <typename VT, typename... Args>
void join(vector<VT> &target, const vector<VT> &arg, const vector<Args> &... args) {
    join(target, arg);
    join(target, args...);
}

// Think that A = A^(T) > 0
template <typename M, typename V>
std::vector<double> conGradSolver(Sparse<M> &A, const std::vector<V> &b, map<int, int> &list_of_neighbors,
                                  vector<set<int>> &send, vector<set<int>> &recv, int n_own, NodesInfo *nodesinfo,
                                  size_t threadsNum = 4, size_t n_max = 100, double eps = 0.0000001);
} // namespace vmo

#endif
