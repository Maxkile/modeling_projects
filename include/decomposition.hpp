#pragma once

#include "stdafx.hpp"
#include <map>

namespace decomp
{
	 // Mesh decomposing(only decart meshes are supported yet)
    vector<int> decomposeMesh(int Nx, int Ny, int Px, int Py, int px, int py){
        vector<int> coord(4);// ibeg, iend; jbeg, jend indexes

        if ((px >= Px) || (py >= Py) || (px < 0) || (py < 0))
        {
            cerr << "Decomposition error: Wrong \'px\' or \'py\' values!" << endl;
            return coord;
        }
        else if ((px >= Px) || (py >= Py))
        {
            cerr << "Decomposition error: Part index is more than part number!" << endl;
            return coord;
        }
        else if ((Px > Nx - 1) || (Py > Ny - 1))
        {
            cerr << "Decomposition error: Number of parts is more than elements themselves!" << endl;
            return coord;
        }
        else
        {
            int xPartSize = (Nx - 1) / Px;//Nx nodes -> (Nx - 1) elements on OX
            int yPartSize = (Ny - 1) / Py;//Ny nodes -> (Ny - 1) elements on OY

            int xLeft = (Nx - 1) % Px;
            int yLeft = (Ny - 1) % Py;

            if ((px < xLeft) && (py < yLeft))//add left element to coord[x] and coord[y]
            {
                coord[0] = xPartSize * px + px;
                coord[1] = xPartSize * (px + 1) + px + 1;
                coord[2] = yPartSize * py + py;
                coord[3] = yPartSize * (py + 1) + py + 1;
            }
            else if (py < yLeft)//add left element to coord[y]
            {
                coord[0] = xPartSize * px + xLeft;
                coord[1] = xPartSize * (px + 1) + xLeft;
                coord[2] = yPartSize * py + py;
                coord[3] = yPartSize * (py + 1) + py + 1;
            }
            else if (px < xLeft)//add left element to coord[x]
            {
                coord[0] = xPartSize * px + px;
                coord[1] = xPartSize * (px + 1) + px + 1;
                coord[2] = yPartSize * py + yLeft;
                coord[3] = yPartSize * (py + 1) + yLeft;
            }
            else//no extra element added
            {
                coord[0] = xPartSize * px + xLeft;
                coord[1] = xPartSize * (px + 1) + xLeft;
                coord[2] = yPartSize * py + yLeft;
                coord[3] = yPartSize * (py + 1) + yLeft;
            }

            return coord;

        }
    }

    // Mesh decomposing(only decart meshes are supported yet). Get all "ibeg, iend; jbeg, jend" for all indexes [0,Px - 1],[0, Py - 1]
    FixedSizeMeshContainer<int> decomposeMesh(int Nx, int Ny, int Px, int Py){

        FixedSizeMeshContainer<int> coords(4);
        coords.reserve(4 * Px * Py);

		if ((Px > Nx - 1) || (Py > Ny - 1))
        {
            cerr << "Decomposition error: Number of parts is more than elements themselves!" << endl;
            return coords;
        }
        else
        {
        	for(int px = 0; px < Px; ++px)
	        {
	            for(int py = 0; py < Py; ++py)
	            {
	                coords.add(decomposeMesh(Nx, Ny, Px, Py, px, py));
	            }
	        }
        	coords.printContainer();

        	return coords;
        }
    }


    void getGlobalIndexes(map<int,int>& G2L, vector<int>& global)
    {
        global.clear();

        for(map<int,int>::iterator map_iter = G2L.begin(); map_iter != G2L.end(); ++map_iter)//forming global
        {
            global.push_back(map_iter->first);
        }
    }


    void getLocalIndexes(map<int,int>& G2L, vector<int>& local)
    {
        local.clear();

        for(map<int,int>::iterator map_iter = G2L.begin(); map_iter != G2L.end(); ++map_iter)//forming local
        {
            local.push_back(map_iter->second);
        }
    }

}
