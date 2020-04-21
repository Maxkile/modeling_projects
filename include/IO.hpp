#pragma once
#include "stdafx.hpp"

#include "platformDependencies.hpp"
#include "FixedSizeMeshContainer.hpp"
#include "VariableSizeMeshContainer.hpp"


void draw_mesh(int Nx, int Ny, int k3, int k4);

int num_elem(int Nx, int Ny, int k3, int k4);

template <typename T1,typename T2>
void write_file(FixedSizeMeshContainer<T1>& C, VariableSizeMeshContainer<T2>& topoEN, VariableSizeMeshContainer<T2>& topoSN,const std::string& path = ""){
    ofstream fout;

	if (path != "")
	{
        if (mkdir(path.c_str(),0777) == EEXIST)
        {
            std::cerr << "Directory exists!" << std::endl;
        }
	}

    fout.open((path + "/" "mesh.txt").c_str());
    fout << "Nn " << C.getBlockNumber() << '\n';
    fout << "Nt " << topoEN.getBlockNumber() << '\n';
    fout << "NFaceBC " << topoSN.getBlockNumber() << '\n';
    fout << "NumCoords " << C.getBlockSize() << '\n';
    fout.close();

    fout.open((path + "/" + "coordinate.msh").c_str());
    for (size_t i = 0; i < C.getBlockNumber(); i++)
        fout << C[i][0] << " " << C[i][1] << '\n';
    fout.close();

    fout.open((path + "/" + "topo.msh").c_str());
    for (size_t i = 0; i < topoEN.getBlockNumber(); i++) {
        fout << topoEN.getBlockSize(i);
        for (size_t j = 0; j < topoEN.getBlockSize(i); j++)
            fout << " " << topoEN[i][j];
        fout << '\n';
    }
    fout.close();

    fout.open((path + "/" + "bctopo.msh").c_str());
    for (size_t i = 0; i < topoSN.getBlockNumber(); i++){
        fout << topoSN.getBlockSize(i) - 1;
        for (size_t j = 1; j < topoSN.getBlockSize(i); j++)
            fout << " " << topoSN[i][j];
        fout << " " << topoSN[i][0];
        fout << '\n';
    }
    fout.close();

}

template <typename T1,typename T2>
int read_file(FixedSizeMeshContainer<T1> &C, VariableSizeMeshContainer<T2> &topoEN, VariableSizeMeshContainer<T2> &topoSN){
    ifstream fin;
    int nN, nE, NFaceBC, NumCoords, count_node;
    char str;
    vector<int> BlockSize;


    fin.open("mesh.txt");
    if (fin.is_open()){
        for (int i = 0; i < 2; i++) fin >> str;
        fin >> nN;
        for (int i = 0; i < 2; i++) fin >> str;
        fin >> nE;
        for (int i = 0; i < 7; i++) fin >> str;
        fin >> NFaceBC;
        for (int i = 0; i < 9; i++) fin >> str;
        fin >> NumCoords;
    }
    else {
        cout << '\n' << "File \"mesh.txt\" not found" << '\n';
        return 1;
    }
    fin.close();


    fin.open("coordinate.msh");
    C.setBlockSize(NumCoords);
    if (fin.is_open()){
        vector<double> temp;
        double coord;

        for (int i = 0; i < nN; i++){
            for (int j = 0; j < NumCoords; j++){
                fin >> coord;
                temp.push_back(coord);
            }
            C.add(temp);
            temp.clear();
        }
    }
    else {
        cout << '\n' << "File \"coordinate.msh\" not found" << '\n';
        return 1;
    }
    fin.close();

    fin.open("topo.msh");
    if (fin.is_open()){
        vector<int> temp;
        BlockSize.push_back(0);
        topoEN.add(temp, BlockSize);
        BlockSize.clear();
        int coord;

        for (int i = 0; i < nE; i++){
            fin >> count_node;
            BlockSize.push_back(count_node);

            for (int j = 0; j < count_node; j++){
                fin >> coord;
                temp.push_back(coord);
            };

            topoEN.add(temp, BlockSize);

            temp.clear();
            BlockSize.clear();
        }
    }
    else {
        cout << '\n' << "File \"topo.msh\" not found" << '\n';
        return 1;
    }
    fin.close();

    fin.open("bctopo.msh");
    if (fin.is_open()){
        vector<int> temp;
        BlockSize.push_back(0);
        topoSN.add(temp, BlockSize);
        BlockSize.clear();
        int number_edge;
        int coord;

        for (int i = 0; i < NFaceBC; i++){
            fin >> count_node;
            BlockSize.push_back(count_node+1);

            temp.push_back(0);
            for (int j = 0; j < count_node; j++){
                fin >> coord;
                temp.push_back(coord);
            };
            fin >> number_edge;
            temp[0] = number_edge;

            topoSN.add(temp, BlockSize);

            temp.clear();
            BlockSize.clear();
        }

    }
    else {
        cout << '\n' << "File \"bctopo.msh\" not found" << '\n';
        return 1;
    }
    fin.close();


    return 0;
}
