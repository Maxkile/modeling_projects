 #pragma once
#include "stdafx.hpp"

#include <map>
#include <omp.h>

#include "FixedSizeMeshContainer.hpp"
#include "VariableSizeMeshContainer.hpp"
#include "decomposition.hpp"

namespace topos
{
	constexpr int TRIANGLE_NODES = 3;
    constexpr int SQUARE_NODES = 4;

    //coords
    template<typename T>
    void build_coord(FixedSizeMeshContainer<T>& C, int Lx, int Ly, int Nx, int Ny){
        vector<T> temp;

        for (int i = 0; i < Ny; i++)
        {
            for (int j = 0; j < Nx; j++)
            {
                temp.push_back(static_cast<T>((Lx / (Nx - 1)) * j));
                temp.push_back(static_cast<T>((Ly / (Ny - 1)) * j));
            }
        }
		C.add(temp);
    }

    int computeMeshFiguresNumberLeft(int figCount1, int figCount2, int skippedElemsCount, int curMeshFigureStructure)//defines how many triangles and squares left in mesh after skipping elements
    {
        int elemsLeft = skippedElemsCount > curMeshFigureStructure ? (figCount1 + figCount2) - ((skippedElemsCount - curMeshFigureStructure) % (figCount1 + figCount2)) : (figCount1 + figCount2) - curMeshFigureStructure;

        //curMeshFigureStructure - left figures to put on mesh till current time 
        //consider curMeshFigureStructure is ALWAYS <= figCount1 + figCount2

        if (skippedElemsCount > curMeshFigureStructure)
        {
            elemsLeft = (figCount1 + figCount2) - ((skippedElemsCount - curMeshFigureStructure) % (figCount1 + figCount2));
        }
        else if (skippedElemsCount == curMeshFigureStructure)
        {
            elemsLeft = (figCount1 + figCount2) - (skippedElemsCount - curMeshFigureStructure);
        }
        else
        {
            elemsLeft = curMeshFigureStructure - skippedElemsCount;
        }

        return elemsLeft;
    }


    VariableSizeMeshContainer<int> toLocalIndexesTopoEN(VariableSizeMeshContainer<int>& originEN, map<int,int>& G2L)
    {
        vector<int> BlockSize;
        vector<int> temp;
        VariableSizeMeshContainer<int> local(temp,BlockSize);

        for(size_t i = 0; i < originEN.getBlockNumber(); ++i)
        {
            size_t blockSize = originEN.getBlockSize(i);

            for(size_t j = 0; j < blockSize; ++j)
            {
                temp.push_back(G2L[originEN[i][j]]);
            }
            BlockSize.push_back(blockSize);
        }

        local.add(temp,BlockSize);

        return local;
    }

    //topoEN
    //also generates G2L,L2G,halo's,interfaces
    VariableSizeMeshContainer<int> build_topoEN(int Nx, int Ny, int k3, int k4, size_t submesh_id, const vector<pair<size_t,vector<int>>>& submeshes, map<int,int>& G2L, vector<int>& L2G, vector<int>& part,vector<pair<size_t,int>>& haloes,vector<pair<size_t,int>>& interfaces){//Px, px - current submesh params
        
        G2L.clear();
        L2G.clear();
        part.clear();
        haloes.clear();
        interfaces.clear();

        vector<int> BlockSize;
        vector<int> temp;
        VariableSizeMeshContainer<int> topoEN(temp, BlockSize);


        int beg_i = submeshes[submesh_id].second[0];
        int end_i = submeshes[submesh_id].second[1];
        int beg_j = submeshes[submesh_id].second[2];
        int end_j = submeshes[submesh_id].second[3];

        size_t mesh_id;

        double part_time = 0;
        double start;

        if ((beg_i > Nx) || (end_i > Nx) || (beg_j > Ny) || (end_j > Ny) || (end_i <= 0) || (end_i <= 0) || (end_i <= 0) || (end_i <= 0))
        {
            cerr << "Wrong submesh parameters!" << endl;
            return topoEN;
        }
        else if (submesh_id >= submeshes.size())
        {
            cerr << "Wrong submesh id!" << endl;
            return topoEN;
        }

        //including interface elements
        // beg_i = (beg_i > 0) ? beg_i - 1 : beg_i;
        // end_i = (end_i < Nx - 1) ? end_i + 1 : end_i;
        // beg_j = (beg_j > 0) ? beg_j - 1 : beg_j;
        // end_j = (end_j < Ny - 1) ? end_j + 1 : end_j;
        
        int fullElementsSkipped = (Ny - 1) * beg_i + beg_j;//total number of elements is k3 + k4. So two triangles is one element itself
        int meshFigureStructureCur = computeMeshFiguresNumberLeft(k3,k4,fullElementsSkipped,0);
       
        int local_j = 0;
        int cur_i = beg_i;
        int cur_j = beg_j;

        while(cur_i < end_i)
        {
            while(cur_j < end_j)
            {
                G2L.insert(pair<int,int>(Ny * cur_i + cur_j, local_j));
                L2G.push_back(Ny * cur_i + cur_j);

                //forming part(inner nodes)

                start = omp_get_wtime();
                part.push_back(submeshes[submesh_id].first);
                part_time += omp_get_wtime() - start;

                if (meshFigureStructureCur > k4)//triangle
                {
                    temp.push_back(Ny * cur_i + cur_j);
                    temp.push_back(Ny * cur_i + cur_j + 1);
                    temp.push_back(Ny * (cur_i + 1) + cur_j);
                    BlockSize.push_back(TRIANGLE_NODES);

                    temp.push_back(Ny * cur_i + cur_j + 1);
                    temp.push_back(Ny * (cur_i + 1) + cur_j + 1);
                    temp.push_back(Ny * (cur_i + 1) + cur_j);
                    BlockSize.push_back(TRIANGLE_NODES);                    
                }
                else if (meshFigureStructureCur <= k4)//square
                {
                    temp.push_back(Ny * cur_i + cur_j);
                    temp.push_back(Ny * cur_i + cur_j + 1);
                    temp.push_back(Ny * (cur_i + 1) + cur_j + 1);
                    temp.push_back(Ny * (cur_i + 1) + cur_j);

                    BlockSize.push_back(SQUARE_NODES);
                }

                meshFigureStructureCur--; 
                if (meshFigureStructureCur == 0)
                {
                    meshFigureStructureCur = k3 + k4;
                } 

                cur_j++;
                local_j++;
            }

            fullElementsSkipped = (Ny - (cur_j + 1)) + beg_j;
            meshFigureStructureCur = computeMeshFiguresNumberLeft(k3, k4, fullElementsSkipped,meshFigureStructureCur);
                
            //changing y coord
            G2L.insert(pair<int,int>(Ny * cur_i + cur_j, local_j));
            L2G.push_back(Ny * cur_i + cur_j);


            //forming halo and part
            mesh_id = decomp::getSubmeshIdByCoords(cur_i,cur_j,submeshes,Nx,Ny);
            start = omp_get_wtime();
            part.push_back(mesh_id);
            part_time += omp_get_wtime() - start;
            haloes.push_back(std::make_pair(mesh_id,Ny * cur_i + cur_j));

            local_j++;
            cur_j = beg_j;
            cur_i++;
        }

        for(cur_j = beg_j; cur_j <= end_j; ++cur_j,++local_j)//last y, we haven't visited it yet, but have to store them too
        {
            G2L.insert(pair<int,int>(Ny * end_i + cur_j, local_j));
            L2G.push_back(Ny * end_i + cur_j);

            //forming halo and part
            mesh_id = decomp::getSubmeshIdByCoords(cur_i,cur_j,submeshes,Nx,Ny);
            start = omp_get_wtime();
            part.push_back(mesh_id);
            part_time += omp_get_wtime() - start;
            haloes.push_back(std::make_pair(mesh_id,Ny * end_i + cur_j));

        }
        std::cout << "Total part forming time: " << part_time << " sec" <<  std::endl;

        topoEN.add(temp, BlockSize);
        return topoEN;
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
            }

            if (i >= Nx) {
                temp.push_back(i);
                temp.push_back(i - Nx);
                BlockSize.push_back(2);
            }

            if ((i % Nx != Nx - 1) && (i >= Nx)){
                if (k < k3) {
                    temp.push_back(i);
                    temp.push_back(i - Nx + 1);
                    BlockSize.push_back(2);
                }

                ++k;

                if (k == k3 + k4)
                    k = 0;  
            }
        }
		topoSN.add(temp, BlockSize);

        return topoSN;
    } 

    //reversed topology
    template<typename T>
    VariableSizeMeshContainer<T> build_reverse_topo(const VariableSizeMeshContainer<T>& topo)
    {
        vector<int> BlockSize;
        vector<T> temp;
        vector<int> nt_i_mas;
        vector<int> count_mass;
        vector<int> count_i_mas;

        VariableSizeMeshContainer<T> reverse_topo(temp, BlockSize);
        size_t nE = topo.getBlockNumber();
        
        int nN = 0;
        for (size_t i = 0; i < nE; i++){
            for (size_t j = 0; j < topo.getBlockSize(i); j++)
                if (topo[i][j]>nN)
                    nN = topo[i][j];
        }
        nN++;
		
        count_i_mas.reserve(nN);
        count_mass.reserve(nN);
    
		
                
        for (int i = 0; i < nN; i++){
            count_mass.push_back(0);
            count_i_mas.push_back(0);
		}


        for (size_t i = 0; i < nE; i++){
            for (size_t j = 0; j < topo.getBlockSize(i); j++) {
                count_mass[topo[i][j]]++;
                count_i_mas[topo[i][j]]++;
            }
        }
        for (int i = 0; i < nN; i++){
            BlockSize.push_back(count_mass[i]);
            for (int j = 0; j < count_mass[i]; j++)
                temp.push_back(0);
        }
         

        reverse_topo.add(temp,BlockSize);

        for (size_t i = 0; i < nE; i++){
            for (size_t j = 0; j < topo.getBlockSize(i); j++){
                reverse_topo[topo[i][j]][count_mass[topo[i][j]] - count_i_mas[topo[i][j]]] = i;
                count_i_mas[topo[i][j]]--;
            }
        }

        return reverse_topo;
    }

    // topoBSN
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
        }

        for(int i = 0; i < Nx - 1; ++i) {
            temp.push_back(2);
            temp.push_back(Nx + Ny - 2 + i);
            temp.push_back(Nx + Ny - 2 + i + 1);
            BlockSize.push_back(3);
        }

        for(int i = 0; i < Ny - 1; ++i) {
            temp.push_back(3);
            temp.push_back(Nx + Nx + Ny - 3 + i);

            if (i + 1 == Ny - 1)
                temp.push_back(0);
            else
                temp.push_back(Nx + Nx + Ny - 3 + i + 1);
            BlockSize.push_back(3);
        }
        topoBSN.add(temp, BlockSize);

        return topoBSN;
    }

    // topoBNS
    template<typename T>
    VariableSizeMeshContainer<T> build_topoBNS(const VariableSizeMeshContainer<T>& topo){
        vector<int> BlockSize;
        vector<T> temp;
        int Nx, Ny;
        VariableSizeMeshContainer<T> topoBNS(temp, BlockSize);       
        for(Nx = 1; topo[Nx - 1][0] == 0; ++Nx) {}

        for(Ny = 1; topo[Nx + Ny - 2][0] == 1; ++Ny) {}


        for(int i = 0; i < (Nx + Ny - 2) * 2; ++i) {
            if (i == 0) 
                temp.push_back(2*(Nx + Ny - 2) - 1);
            else 
                temp.push_back(i - 1);  

            temp.push_back(i);
            BlockSize.push_back(2);           
        }

        topoBNS.add(temp, BlockSize);

        return topoBNS;
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
            for (size_t j = 0; j < topoSN.getBlockSize(i); j++)
                if (topoSN[i][j]>nN)
                    nN = topoSN[i][j];
        }
        nN++;
       
        count_i_mas.reserve(nN);
        count_mass.reserve(nN);

        for (int i = 0; i < nN; i++){
            count_mass.push_back(0);
            count_i_mas.push_back(0);
        }

        for (int i = 0; i < nS; i++){
            for (size_t j = 0; j < topoSN.getBlockSize(i); j++) {
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
            for (size_t j = 0; j < topoSN.getBlockSize(i); j++){
                k = count_mass[topoSN[i][j]] - count_i_mas[topoSN[i][j]];
                topoNN[topoSN[i][j]][k] = j ? topoSN[i][0] : topoSN[i][1];
                count_i_mas[topoSN[i][j]]--;
            }
        }
        

        return topoNN;
    }

}
