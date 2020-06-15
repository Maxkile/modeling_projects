#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iterator>
#include <map>
#include <omp.h>

#include "FixedSizeMeshContainer.hpp"
#include "IO.hpp"
#include "VariableSizeMeshContainer.hpp"
#include "decomposition.hpp"
#include "parallel.hpp"
#include "solver.hpp"
#include "toposBuild.hpp"
#include "vtkGenerator.hpp"

using namespace std;

constexpr size_t PROCESSOR_NAME_SIZE = 128;
constexpr size_t DEFAULT_THREADS_NUM = 4;

int mpi_initialized = 0;       // flag that MPI is initialized
int proc_id = 0;               // process ID
int proc_number = 1;           // number of processes
MPI_Comm MCW = MPI_COMM_WORLD; // default communicator
int MASTER_ID = 0;             // master process ID

int main(int argc, char **argv) {

    int Nx, Ny, k3, k4;
    int Lx, Ly;
    int Px, Py;

    VariableSizeMeshContainer<int> topoEN;
    VariableSizeMeshContainer<int> topoNE;
    VariableSizeMeshContainer<int> topoNN;
    vector<pair<size_t, vector<int>>> submeshes;

    vector<int> part;
    vector<int> L2G;
    vector<set<int>> send;
    vector<set<int>> recv;
    map<int, int> neighbors;
    size_t n_own = 0; // number of own nodes in 'nodes' vector(offset to haloes)
    vector<int> nodes;

    int type = 0;

    NodesInfo *nodesInfo = nullptr; // used only in terminal input. Contains info about nodes to operate with
    SolverInfo* solverInfo = nullptr;

    double start, end, start_temp, end_temp;

    //       Èíèöèàëèçàöèÿ ïàðàëëåëüíîé ÷àñòè
    //----------------------------------------------//
    int mpi_res;
    mpi_res = MPI_Init(&argc, &argv);
    if (mpi_res != MPI_SUCCESS) {
        cerr << "MPI_Init failed!" << endl;
        exit(1);
    }

    mpi_res = MPI_Comm_rank(MCW, &proc_id);
    if (mpi_res != MPI_SUCCESS) {
        cerr << "Communicator rank failed!"
             << "Code: " << mpi_res << endl;
        exit(1);
    }

    mpi_res = MPI_Comm_size(MCW, &proc_number);
    if (mpi_res != MPI_SUCCESS) {
        cerr << "Communicator size failed!"
             << "Code: " << mpi_res << endl;
        exit(1);
    }
    char proc_name[PROCESSOR_NAME_SIZE];
    int proc_name_size = PROCESSOR_NAME_SIZE;
    mpi_res = MPI_Get_processor_name(proc_name, &proc_name_size);
    if (mpi_res != MPI_SUCCESS) {
        cerr << "Get processor name failed!"
             << "Code: " << mpi_res << endl;
        exit(1);
    }
    mpi_initialized = 1;
    //----------------------------------------------//

    if ((argc == 1) ||
        !((argc == 11) || (argc == 12) || (argc == 13) || (argc == 14) || (argc == 15) || (strcmp(argv[1], "-g")))) {
        if (proc_id == MASTER_ID) {

            cout << "Concurrent parallel solver(the conjugate gradient method) v1.2" << endl;
            cout << "Usage: -g <Lx> <Ly> <Nx> <Ny> <k3> <k4> <Px> <Py> <sparse> (-n <threads>) (<method>) (<output>)"
                 << endl;
            cout << "Commands:" << endl;
            cout << "\t-g                Generate mesh key" << endl;
            cout << "\t-t               Number of threads key" << endl;
            cout << "Options:" << endl;
            cout << "\t<Lx>                      Mesh X length" << endl;
            cout << "\t<Ly>                      Mesh Y length" << endl;
            cout << "\t<Nx>                      Number of X nodes" << endl;
            cout << "\t<Ny>                      Number of Y nodes" << endl;
            cout << "\t<k3>                      Number of squares in sequence" << endl;
            cout << "\t<k4>                      Number of triangles in sequence" << endl;
            cout << "\t<Px>                      Decomposition parts along OX axis" << endl;
            cout << "\t<Py>                      Decomposition parts along OY axis" << endl;
            cout << "\t<sparse>                  Sparse matrix format. Choose between "
                    "\"-csr\",\"-ellp\" and \"-coo\"."
                 << endl;
            cout << "\t<threads>                 Threads to parallel vector matrix operation. By default the value is 4"
                 << endl;
            cout << "\t<method>                  Which way to collect output decision vector. By default processes are "
                    "synchronized and output is blocking."
                    "\n\t\t\t\t  Type \"-gather\"  to gather output vectors to main process and sort before output"
                 << endl;
            cout << "\t<output>                  Where to write.By default writing to stdout. Type \"-file\" to write "
                    "to \"solution.txt\""
                 << endl;
        }
        MPI_Finalize();
        return 1;
    }

    try {
        // creating mesh
        Lx = atoi(argv[2]);
        Ly = atoi(argv[3]);
        Nx = atoi(argv[4]);
        Ny = atoi(argv[5]);

        k3 = atoi(argv[6]);
        k4 = atoi(argv[7]);

        Px = atoi(argv[8]);
        Py = atoi(argv[9]);

        if ((Lx <= 0) || (Ly <= 0) || (k3 < 0) || (k4 < 0)) {
            parallel::printf_master(proc_id, "\nWrong sizes: <Lx>,<Ly>,<Nx>,<Ny> can't be less or equal to "
                                             "zero! <k3>,<k4> must be greater than one! Decomposition parts "
                                             "must be greater than\n");
            MPI_Finalize();
            return 1;
        } else if ((Px <= 0) || (Py <= 0) || (Px > Nx - 1) || (Py > Ny - 1)) {
            parallel::printf_master(proc_id, "\nWrong sizes: <Px>,<Py> can't be less or equal to zero! "
                                             "<Px>,<Py> can't be greater or equal than <Nx> and <Ny>!\n");
            MPI_Finalize();
            return 1;
        }
        if (k3 * k3 + k4 * k4 == 0) {
            parallel::printf_master(proc_id, "\nWrong sizes: <k3> <k4>\n");
            MPI_Finalize();
            return 1;
        }

        if (proc_number != Px * Py) {
            parallel::printf_master(
                proc_id, "\nNumber or processes \"-n\" must fit the expression: Px * Py! But %d != %d * %d\n",
                proc_number, Px, Py);
            MPI_Finalize();
            return 1;
        }
        map<int, int> G2L;

        submeshes = decomp::decomposeMesh(Nx, Ny, Px, Py);

        start_temp = omp_get_wtime();
        topoEN = topos::build_topoEN(Nx, Ny, k3, k4, proc_id, submeshes, G2L, L2G, nodes, part, n_own);
        end_temp = omp_get_wtime();

        // Local mapping
        nodes = topos::toLocalIndexes(nodes, G2L);
        topoEN = topos::toLocalIndexes(topoEN, G2L);

        /////////////////////////////////////////////////Logging
        {
            fflush(stdout);
            parallel::barrier();
            for (int p = 0; p < proc_number; p++) {
                parallel::barrier();

                if (proc_id != p)
                    continue;

                cout << endl;
                printf(LINE_SEPARATOR);
                cout << "\t\t Process: " << proc_id << endl;
                printf(LINE_SEPARATOR);

                cout << "Time:" << endl;
                cout << "\ttopoEN: " << end_temp - start_temp << " sec" << endl;

                start = omp_get_wtime();
                topoNN = topos::build_topoNN_from_topoEN(topoEN);
                end = omp_get_wtime();

                cout << "\ttopoNN: " << end - start << " sec" << endl;
                cout << "\nMemory:" << endl;
                cout << "\ttopoEN: " << (topoEN).getTotalSize() * sizeof(int) << " bytes" << endl;
                cout << "\ttopoNN: " << (topoNN).getTotalSize() * sizeof(int) << " bytes" << endl;

                parallel::build_list_of_neighbors(neighbors, part, proc_id);
                parallel::build_list_send_recv(topoNN, part, neighbors, send, recv, proc_id);

                fflush(stdout);
                parallel::barrier();
            }
        }
        /////////////////////////////////////////////////////

        // Deallocating
        G2L.clear();

        // Allocating nodesInfo struct
        nodesInfo = new NodesInfo(nodes.size(), NodesType::ALL); // default

        solverInfo = new SolverInfo();

        // Default params
        solverInfo->threadsNumber = DEFAULT_THREADS_NUM;
        solverInfo->outputType = OutputType::STDOUT;
        solverInfo->outputStrategy = OutputStrategy::SEPARATE;

        // Solver itself

        if (!strcmp("-csr", argv[10])) {
            solverInfo->sparseMatrixType = SparseType::CSR;
        } else if (!strcmp("-ellp", argv[10])) {
            solverInfo->sparseMatrixType = SparseType::ELLPACK;
        } else if (!strcmp("-coo", argv[10])) {
            solverInfo->sparseMatrixType = SparseType::COO;
        } else {
            parallel::printf_master(proc_id, "Wrong sparse matrix format chosen!");
            delete nodesInfo;
            return 1;
        }

        if (argc == 12) { // no threads number parameter(it needs -t key)
            if (!strcmp(argv[11], "-gather")) {
                solverInfo->outputStrategy = OutputStrategy::GATHER;
            } else if (!strcmp(argv[11], "-file")) {
                solverInfo->outputType = OutputType::FILE;
            } else {
                parallel::printf_master(proc_id, "Wrong threads number or thread key usage!");
                delete nodesInfo;
                return 1;
            }
        } else if (argc == 13) {
            if (!strcmp(argv[11], "-gather")) {
                solverInfo->outputStrategy = OutputStrategy::GATHER;
                if (!strcmp(argv[12], "-file")) {
                    solverInfo->outputType = OutputType::FILE;
                } else {
                    parallel::printf_master(proc_id, "Wrong output type parameter!");
                    delete nodesInfo;
                    return 1;
                }
            } else if (!strcmp(argv[11], "-file")) {
                solverInfo->outputType = OutputType::FILE;
                if (!strcmp(argv[12], "-gather")) {
                    solverInfo->outputStrategy = OutputStrategy::GATHER;
                } else {
                    parallel::printf_master(proc_id, "Wrong output strategy type parameter!");
                    delete nodesInfo;
                    return 1;
                }
            } else if (!strcmp(argv[11], "-t")) {
                solverInfo->threadsNumber = atoi(argv[12]);
            } else {
                parallel::printf_master(proc_id, "Wrong threads number or thread key usage!");
                delete nodesInfo;
                return 1;
            }
        } else if (argc == 14) {
            if (!strcmp(argv[11], "-t")) {
                solverInfo->threadsNumber = atoi(argv[12]);
                if (!strcmp(argv[13], "-gather")) {
                    solverInfo->outputStrategy = OutputStrategy::GATHER;
                } else if (!strcmp(argv[13], "-file")) {
                    solverInfo->outputType = OutputType::FILE;
                } else {
                    parallel::printf_master(proc_id, "Wrong parameter: \"%s\"!", argv[13]);
                    delete nodesInfo;
                    return 1;
                }
            } else if (!strcmp(argv[11], "-gather")) {
                solverInfo->outputStrategy = OutputStrategy::GATHER;
                if (!strcmp(argv[12], "-t")) {
                    solverInfo->threadsNumber = atoi(argv[13]);
                } else {
                    parallel::printf_master(proc_id, "Wrong threads number or thread key usage!");
                    delete nodesInfo;
                    return 1;
                }
            } else if (!strcmp(argv[11], "-file")) {
                solverInfo->outputType = OutputType::FILE;
                if (!strcmp(argv[12], "-t")) {
                    solverInfo->threadsNumber = atoi(argv[13]);
                } else {
                    parallel::printf_master(proc_id, "Wrong threads number or thread key usage!");
                    delete nodesInfo;
                    return 1;
                }
            } else {
                parallel::printf_master(proc_id, "Wrong parameter: \"%s\"!", argv[11]);
                delete nodesInfo;
                return 1;
            }
        } else if (argc == 15) {
            if (!strcmp(argv[11], "-t")) {
                solverInfo->threadsNumber = atoi(argv[12]);
                if (!strcmp(argv[13], "-gather")) {
                    solverInfo->outputStrategy = OutputStrategy::GATHER;
                    if (!strcmp(argv[14], "-file")) {
                        solverInfo->outputType = OutputType::FILE;
                    } else {
                        parallel::printf_master(proc_id, "Wrong output type parameter!");
                        delete nodesInfo;
                        return 1;
                    }
                } else if (!strcmp(argv[13], "-file")) {
                    solverInfo->outputType = OutputType::FILE;
                    if (!strcmp(argv[14], "-gather")) {
                        solverInfo->outputStrategy = OutputStrategy::GATHER;
                    } else {
                        parallel::printf_master(proc_id, "Wrong output strategy type parameter!");
                        delete nodesInfo;
                        return 1;
                    }
                }
            } else if (!strcmp(argv[11], "-gather")) {
                solverInfo->outputStrategy = OutputStrategy::GATHER;
                if (!strcmp(argv[12], "-file")) {
                    solverInfo->outputType = OutputType::FILE;
                    if (!strcmp(argv[13], "-t")) {
                        solverInfo->threadsNumber = atoi(argv[14]);
                    } else {
                        parallel::printf_master(proc_id, "Wrong threads number or thread key usage!");
                        delete nodesInfo;
                        return 1;
                    }
                } else if (!strcmp(argv[12], "-t")) {
                    solverInfo->threadsNumber = atoi(argv[13]);
                    if (!strcmp(argv[14], "-file")) {
                        solverInfo->outputType = OutputType::FILE;
                    } else {
                        parallel::printf_master(proc_id, "Wrong output type parameter!");
                        delete nodesInfo;
                        return 1;
                    }
                }
            } else if (!strcmp(argv[11], "-file")) {
                solverInfo->outputType = OutputType::FILE;
                if (!strcmp(argv[12], "-gather")) {
                    solverInfo->outputStrategy = OutputStrategy::GATHER;
                    if (!strcmp(argv[13], "-t")) {
                        solverInfo->threadsNumber = atoi(argv[14]);
                    } else {
                        parallel::printf_master(proc_id, "Wrong threads number or thread key usage!");
                        delete nodesInfo;
                        return 1;
                    }
                } else if (!strcmp(argv[12], "-t")) {
                    solverInfo->threadsNumber = atoi(argv[13]);
                    if (!strcmp(argv[14], "-gather")) {
                        solverInfo->outputStrategy = OutputStrategy::GATHER;
                    } else {
                        parallel::printf_master(proc_id, "Wrong output type parameter!");
                        delete nodesInfo;
                        return 1;
                    }
                }
            } else {
                parallel::printf_master(proc_id, "Wrong parameter: \"%s\"!", argv[11]);
                delete nodesInfo;
                return 1;
            }
        }

        parallel::barrier();
        parallel::printf_master(proc_id, LINE_SEPARATOR, "\t\t Gather writing", LINE_SEPARATOR);
        parallel::printf_master(proc_id, "Number of processes: %d\n", proc_number);
        parallel::printf_master(proc_id, "Master process %d started on %s\n", MASTER_ID, proc_name);
        parallel::barrier();

        topoEN.clear();
        topoNE.clear();

        if (solverInfo->threadsNumber < 1) {
            cout.flush();
            parallel::printf_master(proc_id, "\nThreads number can't be lower than 1!\n");
            MPI_Finalize();
            return 1;
        }

        parallel::barrier();

        size_t totalSize = static_cast<size_t>(Nx) * Ny;
        solver::solveFromTopoNN(topoNN, nodesInfo, neighbors, send, recv, L2G, n_own, totalSize, end, solverInfo);
        parallel::barrier();
        cout.flush();
        parallel::printf_master(proc_id, "\nSolution time: %.6g\n", end);

        delete solverInfo;
        delete nodesInfo;
    } catch (exception &exc) {
        cout << exc.what() << endl;
        delete solverInfo;
        delete nodesInfo;
        return 1;
    }

#ifdef _WIN32
    if (proc_id == MASTER_ID) {
        cout << endl;
        system("pause");
    }
#elif _WIN64
    if (proc_id == MASTER_ID) {
        cout << endl;
        system("pause");
    }
#endif

    MPI_Finalize();

    return 0;
}
