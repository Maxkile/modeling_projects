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
#include "solver.hpp"
#include "toposBuild.hpp"
#include "vtkGenerator.hpp"

using namespace std;

int main(int argc, char **argv) {

    int Nx, Ny, k3, k4;
    int Lx, Ly;
    int Px, Py;

    FixedSizeMeshContainer<double> C; // coordinates
    VariableSizeMeshContainer<int> topoEN;

    VariableSizeMeshContainer<int> topoNE;
    VariableSizeMeshContainer<int> topoSN;
    VariableSizeMeshContainer<int> topoNS;
    VariableSizeMeshContainer<int> topoNN_1;
    VariableSizeMeshContainer<int> topoNN_2;
    int type = 0;

    NodesInfo *nodesInfo = nullptr; // used only in terminal input. Contains info about nodes to operate with

    double start, end;

    if ((argc == 1) || !((argc == 5) || (argc == 11) || (argc == 12) || (argc == 7) || (argc == 13)) ||
        (strcmp(argv[argc - 1], "--print") && strcmp(argv[argc - 1], "--out") && strcmp(argv[argc - 2], "--out") &&
         strcmp(argv[argc - 2], "--vtk") && strcmp(argv[argc - 3], "--solver"))) {

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
                cout << "Wrong sizes: <Lx>,<Ly>,<Nx>,<Ny> can't be less or equal to "
                        "zero! <k3>,<k4> must be greater than one! Decomposition parts "
                        "must be greater than"
                     << endl;
                return 1;
            } else if ((Px <= 0) || (Py <= 0) || (Px > Nx - 1) || (Py > Ny - 1)) {
                cout << "Wrong sizes: <Px>,<Py> can't be less or equal to zero! "
                        "<Px>,<Py> can't be greater or equal than <Nx> and <Ny>!"
                     << endl;
                return 1;
            }
            if (k3 * k3 + k4 * k4 == 0) {
                cout << "Wrong sizes: <k3> <k4>" << endl;
                return 1;
            }
            C.setBlockSize(2);

            cout << "Time:" << endl;

            start = omp_get_wtime();
            topos::build_coord(C, Lx, Ly, Nx, Ny);
            end = omp_get_wtime();
            cout << "\tC:      " << end - start << " sec" << endl;

            map<int, int> G2L;
            vector<int> L2G;

            vector<int> part;

            vector<int> nodes;

            vector<pair<size_t, vector<int>>> submeshes = decomp::decomposeMesh(Nx, Ny, Px, Py);
            size_t n_own; // number of own nodes in 'nodes' vector(offset to haloes)

            start = omp_get_wtime();
            topoEN = topos::build_topoEN(Nx, Ny, k3, k4, 1, submeshes, G2L, L2G, nodes, n_own);
            end = omp_get_wtime();

            decomp::formPart(part, nodes, submeshes, Nx, Ny);

            // Local mapping
            nodes = topos::toLocalIndexes(nodes, G2L);
            topoEN = topos::toLocalIndexes(topoEN, G2L);
            for (auto i = G2L.cbegin(); i != G2L.cend(); ++i) {
                cout << i->first << " ";
            }
            cout << endl;

            /////////////////////////////////////////////////Logging
            {
                cout << "Part: " << endl;
                for (auto i = part.begin(); i != part.end(); ++i) {
                    std::cout << *i << " ";
                }
                std::cout << std::endl;

                cout << "Nodes vector in local numeration:" << endl;
                cout << "Own nodes: " << n_own << ". Haloes nodes: " << nodes.size() - n_own << endl;
                for (size_t i = 0; i < nodes.size(); ++i) {
                    std::cout << nodes[i] << " ";
                }
                std::cout << std::endl;

                cout << "G2L: " << endl;
                for (auto i = G2L.cbegin(); i != G2L.cend(); ++i) {
                    std::cout << i->first << "->" << i->second << " ";
                }
                std::cout << std::endl;

                cout << "L2G: " << endl;
                for (size_t i = 0; i < L2G.size(); ++i) {
                    std::cout << L2G[i] << " ";
                }
                std::cout << std::endl;

                cout << "\ttopoEN: " << end - start << " sec" << endl;
            }
            ////////////////////////////////////////////////////////////

            // Deallocating
            G2L.clear();
            topoEN.printContainer();
            // Allocating nodesInfo struct
            nodesInfo = new NodesInfo(nodes.size(), NodesType::ALL); // default

            //            start = omp_get_wtime();
            //            topoNE = topos::build_reverse_topo(topoEN);
            //            end = omp_get_wtime();
            cout << "\ttopoNE: " << end - start << " sec" << endl;

            start = omp_get_wtime();
            topoSN = topos::build_topoSN(Nx, Ny, k3, k4);
            end = omp_get_wtime();
            cout << "\ttopoSN: " << end - start << " sec" << endl;

            start = omp_get_wtime();
            topoNS = topos::build_reverse_topo(topoSN);
            end = omp_get_wtime();
            cout << "\ttopoNS: " << end - start << " sec" << endl;

            start = omp_get_wtime();
            topoNN_1 = topos::build_topoNN_from_topoSN(topoSN);
            end = omp_get_wtime();
            cout << "\ttopoNN_1: " << end - start << " sec" << endl;

            //            start = omp_get_wtime();
            //            topoNN_2 = topos::build_topoNN_from_topoEN(topoEN);
            //            end = omp_get_wtime();
            cout << "\ttopoNN_2: " << end - start << " sec" << endl;

        } else if (!strcmp(argv[1], "--file")) {
            if (read_file(C, topoEN, topoSN)) {
                cout << "File can't be read!" << endl;
                return 1;
            }
            type = 1;
            topoNE = topos::build_reverse_topo(topoEN);
            topoNS = topos::build_reverse_topo(topoSN);
            topoNN_1 = topos::build_topoNN_from_topoSN(topoSN);
            topoNN_2 = topos::build_topoNN_from_topoEN(topoEN);
        } else {
            cout << "Expected \"--gen\" or \"--file\"" << endl;
            return 1;
        }

        cout << "\nMemory:" << endl;
        cout << "\tC:      " << (C).getTotalSize() * sizeof(double) << " bytes" << endl;
        cout << "\ttopoEN: " << (topoEN).getTotalSize() * sizeof(int) << " bytes" << endl;
        cout << "\ttopoNE: " << (topoNE).getTotalSize() * sizeof(int) << " bytes" << endl;
        cout << "\ttopoSN: " << (topoSN).getTotalSize() * sizeof(int) << " bytes" << endl;
        cout << "\ttopoNS: " << (topoNS).getTotalSize() * sizeof(int) << " bytes" << endl;
        cout << "\ttopoNN_1: " << (topoNN_1).getTotalSize() * sizeof(int) << " bytes" << endl;
        cout << "\ttopoNN_2: " << (topoNN_2).getTotalSize() * sizeof(int) << " bytes" << endl;

        // output mesh
        if (!strcmp(argv[argc - 2], "--vtk")) {
            if (argv[argc - 1] == nullptr) {
                cout << "Expected <filename>" << endl;
                return 1;
            } else {
                std::string path = argv[argc - 1];
                vtkGenerator<double, int> vtk(path);
                vtk.printInVTK(C, topoEN);
            }
        } else if (!strcmp(argv[argc - 2], "--out")) {
            write_file(C, topoEN, topoSN, argv[argc - 1]);
        } else if (!strcmp(argv[argc - 1], "--out")) {
            write_file(C, topoEN, topoSN);
        } else if (!strcmp(argv[argc - 1], "--print")) {
            cout << "\nCoordinates:\n" << endl;
            C.printContainer();
            cout << "\nTopoEN(local):\n" << endl;
            topoEN.printContainer();
            cout << "\nTopoNE:\n" << endl;
            topoNE.printContainer();
            cout << "\nTopoSN:\n" << endl;
            topoSN.printContainer();
            cout << "\nTopoNS:\n" << endl;
            topoNS.printContainer();
            cout << "\nTopoNN_1:\n" << endl;
            topoNN_1.printContainer();
            cout << "\nTopoNN_2:\n" << endl;
            topoNN_2.printContainer();

            if (!type)
                draw_mesh(Nx - 1, Ny - 1, k3, k4);
        } else if (!strcmp(argv[argc - 3], "--solver")) {

            topoEN.clear();
            topoNE.clear();
            topoSN.clear();
            topoNS.clear();

            if (nodesInfo == nullptr) {
                cout << "Please, use terminal input of topologies for solver!" << endl;
                return 1;
            }

            solveFromTopoNN(topoNN_2, argv[argc - 2], nodesInfo, atoi(argv[argc - 1]));
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
    system("pause");
#endif

    return 0;
}
