#include "FixedSizeMeshContainer.hpp"
#include "VariableSizeMeshContainer.hpp"
#include <omp.h>

#pragma once

namespace topos
{

    //coords
    template<typename T>
    void build_coord(FixedSizeMeshContainer<T>& C, T Lx, T Ly, int Nx, int Ny){
        vector<T> temp;

        for (int i = 0; i < Ny; i++) {
            for (int j = 0; j < Nx; j++) {
                temp.push_back((Lx / (Nx - 1)) * j);
                temp.push_back((Ly / (Ny - 1)) * i);
                C.add(temp);
                temp.clear();
            }
        }
    }
    
    //topoEN
    VariableSizeMeshContainer<int> build_topoEN(int Nx, int Ny, int k3, int k4, int nE){
        vector<int> BlockSize;
        vector<int> temp;
        VariableSizeMeshContainer<int> topoEN(temp, BlockSize);

        bool figure = k3 > 0 ? false : true; // 0 - triangle, 1 - square
        int count_figure = figure ? k4 : k3;
        int EN_i = 0;
        int temp_nE = nE;

        while (temp_nE>0) {
            if (!figure) {
                temp.push_back(EN_i);
                temp.push_back(EN_i+1);
                temp.push_back(EN_i+Nx);

                BlockSize.push_back(3);


                temp.push_back(EN_i+1);
                temp.push_back(EN_i+1+Nx);
                temp.push_back(EN_i+Nx);

                BlockSize.push_back(3);


                temp_nE-=2;
            }
            else 
            {
                temp.push_back(EN_i);
                temp.push_back(EN_i+1);
                temp.push_back(EN_i+1+Nx);
                temp.push_back(EN_i+Nx);
                
                BlockSize.push_back(4);

                temp_nE--;
            }
            count_figure--;
            if (count_figure == 0) {
                figure = !figure;
                count_figure = figure ? k4 : k3;
                if (count_figure == 0) {
                    figure = !figure;
                    count_figure = figure ? k4 : k3;
                }
            }
            EN_i++;
            if (EN_i % Nx == Nx - 1)
                EN_i++;
        }

        topoEN.add(temp, BlockSize);
        BlockSize.clear();
        temp.clear();

        temp.clear();
        BlockSize.clear();

        return topoEN;
    }

    //topoNE
    template<typename T>
    VariableSizeMeshContainer<T> build_topoNE(const VariableSizeMeshContainer<T>& topoEN){
        vector<int> BlockSize;
        vector<T> temp;
        vector<int> count_i_mas;
        vector<int> count_mass;

        VariableSizeMeshContainer<T> topoNE(temp, BlockSize);
        int nE = topoEN.getBlockNumber();
        
        int nN = 0;
        for (int i = 0; i < nE; i++){
            for (int j = 0; j < topoEN.getBlockSize(i); j++)
                if (topoEN[i][j]>nN)
                    nN = topoEN[i][j];
        }
        nN++;

        count_i_mas.reserve(nN);
        count_mass.reserve(nN);
        
        for (int i = 0; i < nN; i++){
            count_mass[i] = 0;
            count_i_mas[i] = 0;
        }

        for (int i = 0; i < nE; i++){
            for (int j = 0; j < topoEN.getBlockSize(i); j++) {
                count_mass[topoEN[i][j]]++;
                count_i_mas[topoEN[i][j]]++;
            }
        }
        for (int i = 0; i < nN; i++){
            BlockSize.push_back(count_mass[i]);
            for (int j = 0; j < count_mass[i]; j++)
                temp.push_back(0);
        }
         
        topoNE.add(temp,BlockSize);

        for (int i = 0; i < nE; i++){
            for (int j = 0; j < topoEN.getBlockSize(i); j++){
                topoNE[topoEN[i][j]][count_mass[topoEN[i][j]] - count_i_mas[topoEN[i][j]]] = i;
                count_i_mas[topoEN[i][j]]--;
            }
        }

        return topoNE;
    }

    //topoBSN
    VariableSizeMeshContainer<int> build_topoBSN(int Nx, int Ny){
        vector<int> BlockSize;
        vector<int> temp;
        VariableSizeMeshContainer<int> topoBSN(temp, BlockSize);

        for(int i = 0; i < Nx - 1; ++i) {
            temp.push_back(0);
            temp.push_back(i);
            temp.push_back(i + 1);
            BlockSize.push_back(3);

            topoBSN.add(temp, BlockSize);

            temp.clear();
            BlockSize.clear();
        }

        for(int i = 0; i < Ny - 1; ++i) {
            temp.push_back(1);
            temp.push_back(Nx - 1 + i);
            temp.push_back(Nx - 1 + i + 1);
            BlockSize.push_back(3);

            topoBSN.add(temp, BlockSize);

            temp.clear();
            BlockSize.clear();
        }

        for(int i = 0; i < Nx - 1; ++i) {
            temp.push_back(2);
            temp.push_back(Nx + Ny - 2 + i);
            temp.push_back(Nx + Ny - 2 + i + 1);
            BlockSize.push_back(3);

            topoBSN.add(temp, BlockSize);

            temp.clear();
            BlockSize.clear();
        }

        for(int i = 0; i < Ny - 1; ++i) {
            temp.push_back(3);
            temp.push_back(Nx + Nx + Ny - 3 + i);

            if (i + 1 == Ny - 1)
                temp.push_back(0);
            else
                temp.push_back(Nx + Nx + Ny - 3 + i + 1);
            BlockSize.push_back(3);

            topoBSN.add(temp, BlockSize);

            temp.clear();
            BlockSize.clear();
        }

        return topoBSN;
    }

    //topoBNS
    template<typename T>
    VariableSizeMeshContainer<T> build_topoBNS(const VariableSizeMeshContainer<T>& topoBSN){
        vector<int> BlockSize;
        vector<T> temp;
        int Nx, Ny;
        VariableSizeMeshContainer<T> topoBNS(temp, BlockSize);

        std::cout << "TOPO_BNS" << std::endl;
        std::cout << topoBSN.getTotalSize() << std::endl;     

        for(Nx = 1; topoBSN[Nx - 1][0] == 0; ++Nx) {}

        for(Ny = 1; topoBSN[Nx + Ny - 2][0] == 1; ++Ny) {}

        for(int i = 0; i < (Nx + Ny - 2) * 2; ++i) {
            if (i == 0) 
                temp.push_back(2*(Nx + Ny - 2) - 1);
            else 
                temp.push_back(i - 1);  

            temp.push_back(i);
            BlockSize.push_back(2);   

            topoBNS.add(temp, BlockSize);

            temp.clear();        
            BlockSize.clear();
        }

        return topoBNS;
    }

    //topoSN
    VariableSizeMeshContainer<int> build_topoSN(int Nx, int Ny, int k3, int k4){
        vector<int> BlockSize;
        vector<int> temp;
        int k = 0;
        VariableSizeMeshContainer<int> topoSN(temp, BlockSize);    

        for(int i = 0; i < (Nx * Ny); ++i) {
            if (i % Nx != Nx - 1) {
                temp.push_back(i);
                temp.push_back(i + 1);
                BlockSize.push_back(2);

                topoSN.add(temp, BlockSize);

                temp.clear();
                BlockSize.clear();
            }

            if (i >= Nx) {
                temp.push_back(i);
                temp.push_back(i - Nx);
                BlockSize.push_back(2);

                topoSN.add(temp, BlockSize);

                temp.clear();
                BlockSize.clear();  
            }

            if ((i % Nx != Nx - 1) && (i >= Nx)){
                if (k < k3) {
                    temp.push_back(i);
                    temp.push_back(i - Nx + 1);
                    BlockSize.push_back(2);

                    topoSN.add(temp, BlockSize);

                    temp.clear();
                    BlockSize.clear();    
                }

                ++k;

                if (k == k3 + k4)
                    k = 0;  
            }
        }

        return topoSN;
    } 

    //topoNN
    template<typename T>
    VariableSizeMeshContainer<T> build_topoNN(const VariableSizeMeshContainer<T>& topoSN){
        vector<int> BlockSize;
        vector<T> temp;
        vector<int> count_i_mas;
        vector<int> count_mass;
        int k;

        VariableSizeMeshContainer<T> topoNN(temp, BlockSize);
        int nS = topoSN.getBlockNumber();
        
        int nN = 0;
        for (int i = 0; i < nS; i++){
            for (int j = 0; j < topoSN.getBlockSize(i); j++)
                if (topoSN[i][j]>nN)
                    nN = topoSN[i][j];
        }
        nN++;

        count_i_mas.reserve(nN);
        count_mass.reserve(nN);
        
        for (int i = 0; i < nN; i++){
            count_mass[i] = 0;
            count_i_mas[i] = 0;
        }

        for (int i = 0; i < nS; i++){
            for (int j = 0; j < topoSN.getBlockSize(i); j++) {
                count_mass[topoSN[i][j]]++;
                count_i_mas[topoSN[i][j]]++;
            }
        }
        for (int i = 0; i < nN; i++){
            BlockSize.push_back(count_mass[i]);
            for (int j = 0; j < count_mass[i]; j++)
                temp.push_back(0);
        }
         
        topoNN.add(temp,BlockSize);

        for (int i = 0; i < nS; i++){
            for (int j = 0; j < topoSN.getBlockSize(i); j++){
                k = count_mass[topoSN[i][j]] - count_i_mas[topoSN[i][j]];
                topoNN[topoSN[i][j]][k] = j ? topoSN[i][0] : topoSN[i][1];
                count_i_mas[topoSN[i][j]]--;
            }
        }

        return topoNN;
    }
}
