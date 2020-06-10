#include "solver.hpp"

extern int proc_id;     // process ID
extern int proc_number; // number of processes
extern MPI_Comm MCW;    // default communicator
extern int MASTER_ID;   // master process ID

int solver::solveFromTopoNN(const VariableSizeMeshContainer<int> &topoNN, const char *type_matr, NodesInfo *nodesinfo,
                            map<int, int> &list_of_neighbors, vector<set<int>> &send, vector<set<int>> &recv,
                            vector<int> &L2G, int n_own, size_t totalSize, size_t threadsNumber) {
    std::vector<double> x;
    size_t n;

    Sparse<double> *matrix = nullptr;

    if (threadsNumber < 1) {
        std::cout << "Threads number can't be lower than 1!" << std::endl;
        return 1;
    } else if (!strcmp("-csr", type_matr)) {
        matrix = new SparseCSR<double>(topoNN, L2G);
        if (proc_id == MASTER_ID) {
            std::cout << std::endl << "Type matrix: CSR\n" << endl;
        }
    } else if (!strcmp("-ellp", type_matr)) {
        matrix = new SparseELL<double>(topoNN, L2G);
        if (proc_id == MASTER_ID) {
            std::cout << std::endl << "Type matrix: ELLPACK\n" << endl;
        }
    } else if (!strcmp("-coo", type_matr)) {
        matrix = new SparseCOO<double>(topoNN, L2G);
        if (proc_id == MASTER_ID) {
            std::cout << std::endl << "Type matrix: COORD\n" << endl;
        }
    } else {
        if (proc_id == MASTER_ID) {
            std::cout << "Wrong sparse matrix format chosen!" << std::endl;
            return 1;
        }
    }

    n = matrix->getDenseRows();

    vector<double> b(n);

    for (size_t i = 0; i < n; i++) {
        b[i] = 1 + L2G[i];
    }

    x = vmo::conGradSolver(*matrix, b, list_of_neighbors, send, recv, n_own, nodesinfo, threadsNumber);

    cout.flush();
    parallel::barrier();
    IO::MPI_self_write("solution.txt", x, n_own, L2G);
    parallel::barrier();
    IO::MPI_gather_write("solution_gather.txt", x, n_own, L2G, totalSize);
    parallel::barrier();

    delete matrix;

    return 0;
}
