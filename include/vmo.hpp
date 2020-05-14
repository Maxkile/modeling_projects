#pragma once
#include "stdafx.hpp"

#include <cmath>
#include <iomanip>
#include <omp.h>

#include "Sparse.hpp"
#include "nodeType.hpp"
#include "platformDependencies.hpp"

// vector matrix operations: consider we are having same template arguments
namespace vmo {

static double timeDot = 0;
static double timeMultiply = 0;
static double timeAxpby = 0;

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
VT dot(const std::vector<VT> &x, const std::vector<VT> &y, NodesInfo *nodesinfo, size_t threadsNumber = 4) {
    VT result = 0;
    size_t size = x.size();

    if (size != y.size()) {
        std::cerr << "Incompatible sizes in dot!" << std::endl;
        return result; // incorrent
    } else {
        double start = omp_get_wtime();
#pragma omp parallel for num_threads(threadsNumber) reduction(+ : result)
        for (OPENMP_INDEX_TYPE i = 0; i < size; ++i) {
            result += x[i] * y[i];
        }
        timeDot += omp_get_wtime() - start;
        return result;
    }
}

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
double norm(const std::vector<double> &vec, NodesInfo *nodesinfo);

static size_t i = 0;

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
std::vector<double> conGradSolver(Sparse<M> &A, const std::vector<V> &b, NodesInfo *nodesinfo, size_t threadsNum = 4,
                                  size_t n_max = 100, double eps = 0.0000001) {
    std::vector<double> x_prev(A.getDenseColumns());
    double timeSPmv = 0;

    if (A.getDenseColumns() != b.size()) {
        std::cerr << "Incompatible sizes in Solver!" << std::endl;
    } else {
        size_t size = A.getDenseColumns();

        std::vector<M> diagonal = A.getDiagonal();
        std::vector<M> rev_M(diagonal.size());

        // forming preconditioning matrix - M
        for (size_t i = 0; i < size; ++i) {
            rev_M[i] = 1 / (diagonal[i] + eps);
        }

        SparseELL<double> reverse_M(rev_M); // preconditioner from diagonal vector

        // initializang parameters for algorithm

        std::vector<double> x_cur;

        std::vector<double> r_prev = b; // we think that x0 = (0,0,0,0,0,0 ... 0)^(T) -> b - Ax = b

        std::vector<double> p_cur; // p and p+1 are conjugated relative to A
        std::vector<double> p_prev;

        std::vector<double> q(size);
        std::vector<double> z(size);
        std::vector<double> multiplyResult(size);

        double delta_cur, delta_prev;
        double alpha, beta;
        double b_norm = norm(b, nodesinfo);

        size_t k = 1;

        std::cout << "| Iteration | Norm value |" << endl;

        double conGradSolverTime = omp_get_wtime();

        do {
            reverse_M.spmv(r_prev, z, timeSPmv, nodesinfo, threadsNum);
            delta_cur = dot(r_prev, z, nodesinfo, threadsNum);

            if (k == 1) {
                p_cur = z;
            } else {
                beta = delta_cur / delta_prev;
                multiply(p_prev, multiplyResult, beta, threadsNum);
                axpby(multiplyResult, z, 1, 1, nodesinfo, threadsNum);
                p_cur = multiplyResult;
            }

            A.spmv(p_cur, q, timeSPmv, nodesinfo, threadsNum);
            alpha = delta_cur / dot(p_cur, q, nodesinfo);

            multiply(p_cur, multiplyResult, alpha, threadsNum);
            axpby(x_prev, multiplyResult, 1, 1, nodesinfo, threadsNum);

            multiply(q, multiplyResult, alpha, threadsNum);
            axpby(r_prev, multiplyResult, 1, -1, nodesinfo, threadsNum);

            std::cout << "|" << setw(11) << k << "|" << setw(12) << fixed << setprecision(7)
                      << norm(r_prev, nodesinfo) / b_norm << "|" << endl;

            p_prev = p_cur;
            delta_prev = delta_cur;

            k += 1;

        } while ((norm(r_prev, nodesinfo) / b_norm >= eps) && (k <= n_max));

        conGradSolverTime = omp_get_wtime() - conGradSolverTime;

        std::cout << std::endl;
        std::cout << "Dot functions total time: " << timeDot << std::endl;
        std::cout << "Axpby functions total time: " << timeAxpby << std::endl;
        std::cout << "Multiply functions total time: " << timeMultiply << std::endl;
        std::cout << "Spmv functions total time: " << timeSPmv << std::endl;
        std::cout << std::endl;
        std::cout << "Conjugate gradient method solution total time: " << conGradSolverTime << std::endl;
        std::cout << std::endl;
    }

    return x_prev;
}
} // namespace vmo
