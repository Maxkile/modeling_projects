#include "vmo.hpp"

template <>
double vmo::dot<double>(const std::vector<double> &x, const std::vector<double> &y, int n_own, NodesInfo *nodesinfo,
                        size_t threadsNumber) {
    double result = 0;
    size_t size = n_own;
    double *buf_send = new double[1]; // for All reduce
    double *buf_recv = new double[1]; // for All reduce

    double start = omp_get_wtime();
#pragma omp parallel for num_threads(threadsNumber) reduction(+ : result)
    for (OPENMP_INDEX_TYPE i = 0; i < size; ++i) {
        result += x[i] * y[i];
    }
    timeDot += omp_get_wtime() - start;

    buf_send[0] = result;

    int mpi_res = MPI_Allreduce(buf_send, buf_recv, 1, MPI_DOUBLE, MPI_SUM, MCW);
    if (mpi_res != MPI_SUCCESS) {
        parallel::printf_master(proc_id, "All reduce (dot) error...\n");
        MPI_Finalize();
        exit(1);
    }

    result = buf_recv[0];

    delete[] buf_recv;
    delete[] buf_send;

    return result;
}

// Compute eucledian norm
double vmo::norm(const std::vector<double> &vec, int n_own, NodesInfo *nodesinfo) {
    return sqrt(dot(vec, vec, n_own, nodesinfo));
}

// Think that A = A^(T) > 0
template <>
std::vector<double> vmo::conGradSolver<double, double>(Sparse<double> &A, const std::vector<double> &b,
                                                       map<int, int> &list_of_neighbors, vector<set<int>> &send,
                                                       vector<set<int>> &recv, int n_own, NodesInfo *nodesinfo,
                                                       size_t threadsNum, size_t n_max, double eps) {
    unsigned count = A.getDenseColumns();
    std::vector<double> x_prev(count);
    double timeSPmv = 0;
    double r_prev_norm;
    int mpi_res;

    if (count != b.size()) {
        parallel::printf_master(proc_id, "Incompatible sizes in Solver!\n");
        MPI_Finalize();
        exit(1);
    } else {
        size_t size = count;

        std::vector<double> diagonal = A.getDiagonal();
        std::vector<double> rev_M(diagonal.size());

        // forming preconditioning matrix - M
        for (size_t i = 0; i < size; ++i) {
            rev_M[i] = 1 / (diagonal[i] + eps);
        }

        SparseELL<double> reverse_M(rev_M); // preconditioner from diagonal vector

        // initializing parameters for algorithm
        std::vector<double> r_prev = b; // we think that x0 = (0,0,0,0,0,0 ... 0)^(T) -> b - Ax = b
        std::vector<double> p_cur;      // p and p+1 are conjugated relative to A
        std::vector<double> p_prev;
        std::vector<double> q(size);
        std::vector<double> z(size);
        std::vector<double> multiplyResult(size);

        double delta_cur, delta_prev;
        double alpha, beta;
        double b_norm = sqrt(dot(b, b, n_own, nodesinfo, threadsNum));

        size_t k = 1;

        parallel::printf_master(proc_id, "| Iteration | Norm value |\n");

        do {
            parallel::update_halo(r_prev, list_of_neighbors, send, recv, proc_id, MCW);
            reverse_M.spmv(r_prev, z, timeSPmv, nodesinfo, threadsNum);
            delta_cur = dot(r_prev, z, n_own, nodesinfo, threadsNum);

            if (k == 1) {
                p_cur = z;
            } else {
                beta = delta_cur / delta_prev;
                multiply(p_prev, multiplyResult, beta, threadsNum);
                axpby(multiplyResult, z, 1, 1, nodesinfo, threadsNum);
                p_cur = multiplyResult;
            }

            parallel::update_halo(p_cur, list_of_neighbors, send, recv, proc_id, MCW);
            A.spmv(p_cur, q, timeSPmv, nodesinfo, threadsNum);
            alpha = delta_cur / dot(p_cur, q, n_own, nodesinfo);

            multiply(p_cur, multiplyResult, alpha, threadsNum);
            axpby(x_prev, multiplyResult, 1, 1, nodesinfo, threadsNum);

            multiply(q, multiplyResult, alpha, threadsNum);
            axpby(r_prev, multiplyResult, 1, -1, nodesinfo, threadsNum);

            r_prev_norm = sqrt(dot(r_prev, r_prev, n_own, nodesinfo, threadsNum));

            if (proc_id == MASTER_ID) {
                cout << "|" << setw(11) << k << "|" << setw(12) << fixed << setprecision(7) << r_prev_norm / b_norm
                     << "|" << endl;
            }

            p_prev = p_cur;
            delta_prev = delta_cur;

            k += 1;

        } while (((r_prev_norm / b_norm) >= eps) && (k <= n_max));

        parallel::barrier();
        cout.flush();

        for (int p = 0; p < proc_number; p++) {
            parallel::barrier();

            if (p == proc_id) {
                cout << endl;
                printf(LINE_SEPARATOR);
                cout << "\t\t Process: " << proc_id << endl;
                printf(LINE_SEPARATOR);
                cout << endl;
                cout << "Dot functions total time: " << timeDot << endl;
                cout << "Axpby functions total time: " << timeAxpby << endl;
                cout << "Multiply functions total time: " << timeMultiply << endl;
                cout << "Spmv functions total time: " << timeSPmv << endl;
                cout << endl;
            }
        }
    }

    return x_prev;
}
