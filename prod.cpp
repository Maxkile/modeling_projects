#include <fstream>
#include <iterator>
#include <cstdlib>
#include <cstring>
#include <omp.h>
#include <map>

#include "toposBuild.hpp"
#include "decomposition.hpp"
#include "VariableSizeMeshContainer.hpp"
#include "FixedSizeMeshContainer.hpp"
#include "vtkGenerator.hpp"
#include "IO.hpp"
#include "solver.hpp"

using namespace std;

int main(int argc, char **argv) {
    
    int Nx, Ny, k3, k4;
    int Lx, Ly;
    int Px, Py;

    FixedSizeMeshContainer<double> C;//coordinates
    VariableSizeMeshContainer<int> topoEN;

    VariableSizeMeshContainer<int> localTopoEN;

    VariableSizeMeshContainer<int> topoNE;
    VariableSizeMeshContainer<int> topoSN;
    VariableSizeMeshContainer<int> topoNS;
    VariableSizeMeshContainer<int> topoNN;

    int type = 0;

    double start, end;

    if ((argc == 1) || !((argc == 5) || (argc == 11) || (argc == 12) || (argc == 7) || (argc == 13)) || (strcmp(argv[argc-1], "--print") && strcmp(argv[argc-1], "--out") && strcmp(argv[argc-2], "--out") && strcmp(argv[argc-2], "--vtk") && strcmp(argv[argc-3], "--solver"))){

        cout << "Mesh builder and solver v0.2 beta" << endl;
        cout << "Usage:  --gen <Lx> <Ly> <Nx> <Ny> <k3> <k4> | --file \n\t--print | --out (<path>) | --vtk <filename> | --solver <format> <threads>" <<endl;
        cout << "Commands:" << endl;
        cout << "\t--gen                     Generate mesh" << endl;
        cout << "\t--file                    Read mesh from file(in current directory)" << endl;
        cout << "\t--print                   Print mesh to stdout" << endl;
        cout << "\t--out                     Print mesh to file" << endl;
        cout << "\t--vtk                     Print mesh in \" .vtk \" format to file" << endl;
        cout << "\t--solver                  An example of solving the  Ax = b using the conjugate gradient method " << endl;
        cout << "Options:" << endl;
        cout << "\t<Lx>                      Mesh X length" << endl;
        cout << "\t<Ly>                      Mesh Y length" << endl;
        cout << "\t<Nx>                      Number of X nodes" << endl;
        cout << "\t<Ny>                      Number of Y nodes" << endl;
        cout << "\t<k3>                      Number of squares in sequence" << endl;
        cout << "\t<k4>                      Number of triangles in sequence" << endl;
        cout << "\t<Px>                      Decomposition parts along OX axis" << endl;
        cout << "\t<Py>                      Decomposition parts along OY axis" << endl;
       	cout << "\t<path>                    Full output path, not necessary. By default, write to current directory." << endl;
        cout << "\t<filename>                Vtk output file name" << endl;
        cout << "\t<format>                  Sparse matrix format. Choose between \"-csr\",\"-ellp\" and \"-coo\"." << endl;
        cout << "\t<threads>                 Threads to parallel SPMV operation." << endl;
        return 1;
    }
    try {
        //creating mesh
        if (!strcmp(argv[1], "--gen")){
            Lx = atoi(argv[2]);
            Ly = atoi(argv[3]);
            Nx = atoi(argv[4]);
            Ny = atoi(argv[5]);

            k3 = atoi(argv[6]);
            k4 = atoi(argv[7]);

            Px = atoi(argv[8]);
            Py = atoi(argv[9]);

            if ((Lx <= 0) || (Ly <= 0) || (k3 < 0) || (k4 < 0))
            {
                cout << "Wrong sizes: <Lx>,<Ly>,<Nx>,<Ny> can't be less or equal to zero! <k3>,<k4> must be greater than one! Decomposition parts must be greater than" << endl;
                return 1;
            }
            else if ((Px <= 0) || (Py <= 0) || (Px > Nx - 1) || (Py > Ny - 1))
            {
                cout << "Wrong sizes: <Px>,<Py> can't be less or equal to zero! <Px>,<Py> can't be greater than <Nx> and <Ny>!" << endl;
                return 1;
            }
            if (k3*k3 + k4*k4 == 0) {
                cout << "Wrong sizes: <k3> <k4>" << endl;
                return 1;
            }
            C.setBlockSize(2);

            cout << "Time:" << endl;

            start = omp_get_wtime();
            topos::build_coord(C, Lx, Ly, Nx, Ny);
            end = omp_get_wtime();
            cout << "\tC:      " << end - start << " sec" << endl;

            map<int,int> G2L;
            vector<int> part;

            FixedSizeMeshContainer<int> submeshes = decomp::decomposeMesh(Nx,Ny,Px,Py);

            // vector<int> submeshIndexes_0_0 = decomp::decomposeMesh(Nx, Ny, Px, Py, 0, 0);
            // std::cout << submeshIndexes_0_0[0] << submeshIndexes_0_0[1] << submeshIndexes_0_0[2] << submeshIndexes_0_0[3] << endl;
            // vector<int> submeshIndexes_0_1 = decomp::decomposeMesh(Nx, Ny, Px, Py, 0, 1);
            // std::cout << submeshIndexes_0_1[0] << submeshIndexes_0_1[1] << submeshIndexes_0_1[2] << submeshIndexes_0_1[3] << endl;
            // vector<int> submeshIndexes_1_0 = decomp::decomposeMesh(Nx, Ny, Px, Py, 1, 0);
            // std::cout << submeshIndexes_1_0[0] << submeshIndexes_1_0[1] << submeshIndexes_1_0[2] << submeshIndexes_1_0[3] << endl;
            // vector<int> submeshIndexes_1_1 = decomp::decomposeMesh(Nx, Ny, Px, Py, 1, 1);
            // std::cout << submeshIndexes_1_1[0] << submeshIndexes_1_1[1] << submeshIndexes_1_1[2] << submeshIndexes_1_1[3] << endl;

            // FixedSizeMeshContainer<int> submeshIndexes = decomp::decomposeMesh(Nx,Ny,2,2);


            int px,py;

            px = 2;
            py = 1;

            start = omp_get_wtime();
            topoEN = topos::build_topoEN(Nx, Ny, k3, k4, submeshes[7][0], submeshes[7][1], submeshes[7][2], submeshes[7][3], px, py, G2L, part);
            end = omp_get_wtime();

            cout << "\ttopoEN: " << end - start << " sec" << endl;

            localTopoEN = topos::toLocalIndexesTopoEN(topoEN,G2L);

            start = omp_get_wtime();
            topoNE = topos::build_reverse_topo(localTopoEN);
            end = omp_get_wtime();
            cout << "\ttopoNE: " << end - start << " sec" << endl;

            start = omp_get_wtime();
            topoSN = topos::build_topoSN(Nx,Ny,k3,k4);
            end = omp_get_wtime();
            cout << "\ttopoSN: " << end - start << " sec" << endl;

            start = omp_get_wtime();
            topoNS = topos::build_reverse_topo(topoSN);
            end = omp_get_wtime();
            cout << "\ttopoNS: " << end - start << " sec" << endl;

            start = omp_get_wtime();
            topoNN = topos::build_topoNN(topoSN);
            end = omp_get_wtime();
            cout << "\ttopoNN: " << end - start << " sec" << endl;
        }
        else if (!strcmp(argv[1], "--file")){
            if (read_file(C, topoEN,topoSN))
            {
                cout << "File can't be read!" << endl;
                return 1;
            }
            type = 1;
            topoNE = topos::build_reverse_topo(topoEN);
            topoNS = topos::build_reverse_topo(topoSN);
            topoNN = topos::build_topoNN(topoSN);
        }
        else
        {
            cout << "Expected \"--gen\" or \"--file\"" << endl;
            return 1;
        }

        cout << "\nMemory:" << endl;
        cout << "\tC:      " <<  (C).getTotalSize() * sizeof(double) << " bytes" << endl;
        cout << "\ttopoEN: " <<  (topoEN).getTotalSize() * sizeof(int) << " bytes" << endl;
        cout << "\ttopoNE: " <<  (topoNE).getTotalSize() * sizeof(int) << " bytes" << endl;
        cout << "\ttopoSN: " <<  (topoSN).getTotalSize() * sizeof(int) << " bytes" << endl;
        cout << "\ttopoNS: " <<  (topoNS).getTotalSize() * sizeof(int) << " bytes" << endl;
        cout << "\ttopoNN: " <<  (topoNN).getTotalSize() * sizeof(int) << " bytes" << endl;


        //output mesh
        if (!strcmp(argv[argc-2], "--vtk"))
        {
            if (argv[argc - 1] == nullptr)
            {
                cout << "Expected <filename>" << endl;
                return 1;
            }
            else
            {
                std::string path = argv[argc-1];
                vtkGenerator<double,int> vtk(path);
                vtk.printInVTK(C,topoEN);
            }
        }
        else if (!strcmp(argv[argc-2], "--out")){
            write_file(C,topoEN,topoSN,argv[argc-1]);
        }
        else if (!strcmp(argv[argc-1], "--out")){
            write_file(C,topoEN,topoSN);
        }
        else if (!strcmp(argv[argc-1], "--print")){
            cout << "\nCoordinates:\n" << endl;
            C.printContainer();
            cout << "\nTopoEN:\n" << endl;
            topoEN.printContainer();
            cout << "\nTopoEN(local):\n" << endl;
            localTopoEN.printContainer();
            cout << "\nTopoNE:\n" << endl;
            topoNE.printContainer();
            cout << "\nTopoSN:\n" << endl;
            topoSN.printContainer();
            cout << "\nTopoNS:\n" << endl;
            topoNS.printContainer();
            cout << "\nTopoNN:\n" << endl;
            topoNN.printContainer();

            if (!type) draw_grid(Nx - 1, Ny - 1, k3, k4);
        }
        else if (!strcmp(argv[argc-3], "--solver")){

            topoEN.clear();
            topoNE.clear();
            topoSN.clear();
            topoNS.clear();

            solveFromTopoNN(topoNN, argv[argc-2], atoi(argv[argc-1]));
        }
        else
        {
            cout << "Expected \"--vtk\" <filename> ,\"--file\" or \"--print\" or \"--solver <format> <threads>\"" << endl;
            return 1;
        }
    }
    catch(exception& exc){
        cout << exc.what() << endl;
        return 1;
    }

	#ifdef _WIN32
		system("pause");
	#endif

    return 0;
}
