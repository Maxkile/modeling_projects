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

int mpi_initialized = 0;       // flag that MPI is initialized
int proc_id = 0;               // process ID
int proc_number = 1;           // number of processes
MPI_Comm MCW = MPI_COMM_WORLD; // default communicator
int MASTER_ID = 0;             // master process ID

int main(int argc, char **argv) {

    int Nx, Ny, k3, k4;
    int Lx, Ly;
    int Px, Py;

    FixedSizeMeshContainer<double> C; // coordinates
    VariableSizeMeshContainer<int> topoEN;

    VariableSizeMeshContainer<int> topoNE;
    VariableSizeMeshContainer<int> topoSN;
    VariableSizeMeshContainer<int> topoNS;
    VariableSizeMeshContainer<int> topoNN;
    vector<pair<size_t, vector<int>>> submeshes;
    vector<int> part;
    vector<int> L2G;
    vector<set<int>> send;
    vector<set<int>> recv;
    map<int, int> neighbors;
    size_t n_own; // number of own nodes in 'nodes' vector(offset to haloes)
    vector<int> nodes;

    int type = 0;

    NodesInfo *nodesInfo = nullptr; // used only in terminal input. Contains info about nodes to operate with

    double start, end, start_temp, end_temp;

    //       Инициализация параллельной части
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

    if ((argc == 1) || !((argc == 5) || (argc == 11) || (argc == 12) || (argc == 7) || (argc == 13)) ||
        (strcmp(argv[argc - 1], "--print") && strcmp(argv[argc - 1], "--out") && strcmp(argv[argc - 2], "--out") &&
         strcmp(argv[argc - 2], "--vtk") && strcmp(argv[argc - 3], "--solver"))) {

        if (proc_id == MASTER_ID) {
            cout << "Mesh builder and solver v0.2 beta" << endl;
            cout << "Usage:  --gen <Lx> <Ly> <Nx> <Ny> <k3> <k4> <Px> <Py>| --file "
                    "\n\t--print | --out (<path>) | --vtk <filename> | --solver "
                    "<format> <threads>"
                 << endl;
            cout << "Commands:" << endl;
            cout << "\t--gen                     Generate mesh" << endl;
            cout << "\t--file                    Read mesh from file(in current "
                    "directory)"
                 << endl;
            cout << "\t--print                   Print mesh to stdout" << endl;
            cout << "\t--out                     Print mesh to file" << endl;
            cout << "\t--vtk                     Print mesh in \" .vtk \" format to file" << endl;
            cout << "\t--solver                  An example of solving the  Ax = b "
                    "using the conjugate gradient method "
                 << endl;
            cout << "Options:" << endl;
            cout << "\t<Lx>                      Mesh X length" << endl;
            cout << "\t<Ly>                      Mesh Y length" << endl;
            cout << "\t<Nx>                      Number of X nodes" << endl;
            cout << "\t<Ny>                      Number of Y nodes" << endl;
            cout << "\t<k3>                      Number of squares in sequence" << endl;
            cout << "\t<k4>                      Number of triangles in sequence" << endl;
            cout << "\t<Px>                      Decomposition parts along OX axis" << endl;
            cout << "\t<Py>                      Decomposition parts along OY axis" << endl;
            cout << "\t<path>                    Full output path, not necessary. By "
                    "default, write to current directory."
                 << endl;
            cout << "\t<filename>                Vtk output file name" << endl;
            cout << "\t<format>                  Sparse matrix format. Choose between "
                    "\"-csr\",\"-ellp\" and \"-coo\"."
                 << endl;
            cout << "\t<threads>                 Threads to parallel SPMV operation." << endl;
        }
        MPI_Finalize();
        return 1;
    }

    try {
        // creating mesh
        if (!strcmp(argv[1], "--gen")) {
            Lx = atoi(argv[2]);
            Ly = atoi(argv[3]);
            Nx = atoi(argv[4]);
            Ny = atoi(argv[5]);

            k3 = atoi(argv[6]);
            k4 = atoi(argv[7]);

            Px = atoi(argv[8]);
            Py = atoi(argv[9]);

            if ((Lx <= 0) || (Ly <= 0) || (k3 < 0) || (k4 < 0)) {
                if (proc_id == MASTER_ID) {
                    cout << "Wrong sizes: <Lx>,<Ly>,<Nx>,<Ny> can't be less or equal to "
                            "zero! <k3>,<k4> must be greater than one! Decomposition parts "
                            "must be greater than"
                         << endl;
                }
                MPI_Finalize();
                return 1;
            } else if ((Px <= 0) || (Py <= 0) || (Px > Nx - 1) || (Py > Ny - 1)) {
                if (proc_id == MASTER_ID) {
                    cout << "Wrong sizes: <Px>,<Py> can't be less or equal to zero! "
                            "<Px>,<Py> can't be greater or equal than <Nx> and <Ny>!"
                         << endl;
                }
                MPI_Finalize();
                return 1;
            }
            if (k3 * k3 + k4 * k4 == 0) {
                if (proc_id == MASTER_ID) {
                    cout << "Wrong sizes: <k3> <k4>" << endl;
                }
                MPI_Finalize();
                return 1;
            }

            C.setBlockSize(2);

            map<int, int> G2L;

            submeshes = decomp::decomposeMesh(Nx, Ny, Px, Py);

            if (proc_id == MASTER_ID) {
                // start = omp_get_wtime();
                // topos::build_coord(C, Lx, Ly, Nx, Ny);
                // end = omp_get_wtime();
                // cout << "\nTime:" << endl;
                // cout << "\tC:      " << end - start << " sec" << endl;
                // cout << "\nMemory:" << endl;
                // cout << "\tC:      " << (C).getTotalSize() * sizeof(double) << " bytes" << endl;
            }

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

            // Deallocating
            G2L.clear();

            // Allocating nodesInfo struct
            nodesInfo = new NodesInfo(nodes.size(), NodesType::ALL); // default

        } else if (!strcmp(argv[1], "--file")) {
            if (IO::read_file(C, topoEN, topoSN)) {
                cout << "File can't be read!" << endl;
                return 1;
            }
            type = 1;
            topoNE = topos::build_reverse_topo(topoEN);
            topoNS = topos::build_reverse_topo(topoSN);
            topoNN = topos::build_topoNN_from_topoSN(topoSN);
        } else {
            cout << "Expected \"--gen\" or \"--file\"" << endl;
            return 1;
        }

        // output mesh
        if (!strcmp(argv[argc - 2], "--vtk")) {
            if (argv[argc - 1] == nullptr) {
                cout << "Expected <filename>" << endl;
                return 1;
            } else {
                string path = argv[argc - 1];
                vtkGenerator<double, int> vtk(path);
                vtk.printInVTK(C, topoEN);
            }
        } else if (!strcmp(argv[argc - 2], "--out")) {
            IO::write_file(C, topoEN, topoSN, argv[argc - 1]);
        } else if (!strcmp(argv[argc - 1], "--out")) {
            IO::write_file(C, topoEN, topoSN);
        } else if (!strcmp(argv[argc - 1], "--print")) {
            fflush(stdout);
            parallel::barrier();
            for (int p = 0; p < proc_number; p++) {
                parallel::barrier();

                if (p != proc_id)
                    continue;

                cout << endl;
                printf(LINE_SEPARATOR);
                cout << "\t\t Process: " << proc_id << endl;
                printf(LINE_SEPARATOR);

                if (proc_id == MASTER_ID) {
                    cout << endl << "Decomposition into: " << endl;
                    for (auto it = submeshes.begin(); it != submeshes.end(); ++it) {
                        cout << it->first << ": " << it->second[0] << " " << it->second[1] << " " << it->second[2]
                             << " " << it->second[3] << " " << endl;
                    }
                    cout << endl;
                }

                cout << "Global nodes vector:" << endl;
                cout << "Own nodes: " << n_own << ". Haloes nodes: " << nodes.size() - n_own << endl;
                for (size_t i = 0; i < nodes.size(); ++i) {
                    cout << nodes[i] << " ";
                }
                cout << endl << endl;

                cout << "Part: " << endl;
                for (auto i = part.begin(); i != part.end(); ++i) {
                    cout << *i << " ";
                }
                cout << endl << endl;

                cout << "Local nodes vector:" << endl;
                cout << "Own nodes: " << n_own << ". Haloes nodes: " << nodes.size() - n_own << endl;
                for (size_t i = 0; i < nodes.size(); ++i) {
                    cout << nodes[i] << " ";
                }
                cout << endl << endl;

                cout << "L2G: " << endl;
                for (size_t i = 0; i < L2G.size(); ++i) {
                    cout << L2G[i] << " ";
                }
                cout << endl << endl;

                cout << endl << "List of neighbors: " << endl;
                for (map<int, int>::iterator it = neighbors.begin(); it != neighbors.end(); ++it) {
                    cout << it->first << ": " << it->second << endl;
                }

                cout << endl << "List of send: " << endl;
                for (unsigned i = 0; i < send.capacity(); ++i) {
                    cout << i << ": ";
                    for (set<int>::iterator it = send[i].begin(); it != send[i].end(); ++it) {
                        cout << *it << " ";
                    }
                    cout << endl;
                }

                cout << endl << "List of recv: " << endl;
                for (unsigned i = 0; i < recv.capacity(); ++i) {
                    cout << i << ": ";
                    for (set<int>::iterator it = recv[i].begin(); it != recv[i].end(); ++it) {
                        cout << *it << " ";
                    }
                    cout << endl;
                }

                //            cout << "\nCoordinates:\n" << endl;
                //            C.printContainer();
                cout << "\nTopoEN(local):\n" << endl;
                topoEN.printContainer();
                //            cout << "\nTopoNE:\n" << endl;
                //            topoNE.printContainer();
                //            cout << "\nTopoSN:\n" << endl;
                //            topoSN.printContainer();
                //            cout << "\nTopoNS:\n" << endl;
                //            topoNS.printContainer();
                cout << "\ntopoNN:\n" << endl;
                topoNN.printContainer();
                //            cout << "\nTopoNN_2:\n" << endl;
                //            topoNN_2.printContainer();

                fflush(stdout);
                parallel::barrier();
            }

            if (proc_id == MASTER_ID) {
                if (!type)
                    IO::draw_mesh(Nx - 1, Ny - 1, k3, k4);
            }
        } else if (!strcmp(argv[argc - 3], "--solver")) {
            parallel::barrier();
            if (proc_id == MASTER_ID) {
                cout << endl;
                printf(LINE_SEPARATOR);
                cout << "\t\t SOLUTION " << endl;
                printf(LINE_SEPARATOR);
            }
            parallel::barrier();

            parallel::printf_master(proc_id, "Number of processes: %d\n", proc_number);
            parallel::printf_master(proc_id, "Master process 0 started on %s\n", proc_name);

            topoEN.clear();
            topoNE.clear();
            topoSN.clear();
            topoNS.clear();

            if (nodesInfo == nullptr) {
                parallel::printf_master(proc_id, "Please, use terminal input of topologies for solver!");
                MPI_Finalize();
                return 1;
            }

            parallel::barrier();
            start = omp_get_wtime();

            size_t totalSize = Nx * Ny;
            solver::solveFromTopoNN(topoNN, argv[argc - 2], nodesInfo, neighbors, send, recv, L2G, n_own, totalSize,
                                    atoi(argv[argc - 1]));
            parallel::barrier();
            end = omp_get_wtime();
            parallel::printf_master(proc_id, "\nSolution time: %.6g\n", end - start);

            delete nodesInfo;

        } else {
            cout << "Expected \"--vtk\" <filename> ,\"--file\" or \"--print\" or "
                    "\"--solver <format> <threads>\""
                 << endl;
            return 1;
        }

    } catch (exception &exc) {
        cout << exc.what() << endl;
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
