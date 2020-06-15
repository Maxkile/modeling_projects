#include "solver.hpp"

extern int proc_id;     // process ID
extern int proc_number; // number of processes
extern MPI_Comm MCW;    // default communicator
extern int MASTER_ID;   // master process ID

template <> void solver::defaultInitB<int>(vector<double> &b, const vector<int> &source) {
    b.clear();
    size_t n = source.size();
    b.resize(n);
    for (size_t i = 0; i < n; i++) {
        b[i] = 1 + source[i];
    }
}

int solver::solveFromTopoNN(const VariableSizeMeshContainer<int> &topoNN, NodesInfo *nodesinfo,
                            map<int, int> &list_of_neighbors, vector<set<int>> &send, vector<set<int>> &recv,
                            vector<int> &L2G, int n_own, size_t totalSize, double &timeSpent, SolverInfo *solverInfo) {
    std::vector<double> x;

    double startTime;

    Sparse<double> *matrix = nullptr;

    parallel::printf_master(proc_id, "Threads number: ", solverInfo->threadsNumber);

    switch (solverInfo->sparseMatrixType) {
    case SparseType::CSR:
        matrix = new SparseCSR<double>(topoNN, L2G);
        parallel::printf_master(proc_id, "Type matrix: CSR\n\n");
        break;
    case SparseType::ELLPACK:
        matrix = new SparseELL<double>(topoNN, L2G);
        parallel::printf_master(proc_id, "Type matrix: ELLPACK\n\n");
        break;
    case SparseType::COO:
        matrix = new SparseCOO<double>(topoNN, L2G);
        parallel::printf_master(proc_id, "Type matrix: COO\n\n");
        break;
    }

    vector<double> b;
    defaultInitB(b, L2G);

    startTime = omp_get_wtime();
    x = vmo::conGradSolver(*matrix, b, list_of_neighbors, send, recv, n_own, nodesinfo, solverInfo->threadsNumber);
    timeSpent = omp_get_wtime() - startTime;

    cout.flush();
    parallel::barrier();

    if (solverInfo->outputStrategy == OutputStrategy::GATHER) {
        IO::mpi_gather_write(solverInfo->outputType, x, n_own, L2G, totalSize);
    } else {
        IO::mpi_self_write(solverInfo->outputType, x, n_own, L2G);
    }

    parallel::barrier();
    delete matrix;
    return 0;
}
