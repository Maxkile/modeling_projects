#include "decomposition.hpp"

/*
 * Mesh decomposing(only decart meshes are supported yet),returns pair - submesh id -> submesh coords
*/
pair<size_t,vector<int>> decomp::decomposeMesh(int Nx, int Ny, int Px, int Py, int px, int py){

    pair<size_t,vector<int>> submesh_params;
    size_t id;
    vector<int> coord(4);// ibeg, iend; jbeg, jend indexes

    if ((px >= Px) || (py >= Py) || (px < 0) || (py < 0))
    {
        cerr << "Decomposition error: Wrong \'px\' or \'py\' values!" << endl;
        return submesh_params;
    }
    else if ((px >= Px) || (py >= Py))
    {
        cerr << "Decomposition error: Part index is more than part number!" << endl;
        return submesh_params;
    }
    else if ((Px > Nx - 1) || (Py > Ny - 1))
    {
        cerr << "Decomposition error: Number of parts is more than elements themselves!" << endl;
        return submesh_params;
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

        id = px * Py + py;
        return std::make_pair(id,coord);

    }
}

/*
 *  Mesh decomposing(only decart meshes are supported yet). Get all "ibeg, iend; jbeg, jend" for all indexes [0,Px - 1],[0, Py - 1]
*/
vector<pair<size_t,vector<int>>> decomp::decomposeMesh(int Nx, int Ny, int Px, int Py){

    vector<pair<size_t,vector<int>>> submeshes;

    if ((Px > Nx - 1) || (Py > Ny - 1))
    {
        cerr << "Decomposition error: Number of parts is more than elements themselves!" << endl;
        return submeshes;
    }
    else
    {
        for(int px = 0; px < Px; ++px)
        {
            for(int py = 0; py < Py; ++py)
            {
                submeshes.emplace_back(decomposeMesh(Nx, Ny, Px, Py, px, py));
            }
        }
        return submeshes;
    }
}

/*
 * Get submmesh global id among submeshes
*/
size_t decomp::getSubmeshIdByCoords(int x,int y, const vector<pair<size_t,vector<int>>>& submeshes, int Nx, int Ny)
{
    size_t submesh_id;
    for(auto iter = submeshes.begin(); iter != submeshes.end(); ++iter)
    {
        int beg_x = iter->second[0];
        int end_x = iter->second[1];
        int beg_y = iter->second[2];
        int end_y = iter->second[3];

        if ((end_x == Nx - 1) && (end_y == Ny - 1))//border
        {
            if ((x >= beg_x) && (x <= end_x) && (y >= beg_y) && (y <= end_y))
            {
                submesh_id = iter->first;
                break;
            }
        }
        else if (end_x == Nx - 1)//border
        {
            if ((x >= beg_x) && (x <= end_x) && (y >= beg_y) && (y < end_y))
            {
                submesh_id = iter->first;
                break;
            }
        }
        else if (end_y == Ny - 1)//border
        {
            if ((x >= beg_x) && (x < end_x) && (y >= beg_y) && (y <= end_y))
            {
                submesh_id = iter->first;
                break;
            }
        }
        else
        {
            if ((x >= beg_x) && (x < end_x) && (y >= beg_y) && (y < end_y))
            {
                submesh_id = iter->first;
                break;
            }
        }
    }

    return submesh_id;
}

/*
 * Is node of with (x,y) coords is interface node of 'current_id' submesh
*/
bool decomp::isInterface(int x, int y, const vector<pair<size_t,vector<int>>>& submeshes, int Nx, int Ny, size_t current_id)
{
    return !(decomp::getSubmeshIdByCoords(x,y,submeshes,Nx,Ny) == current_id);
}


/*
 * Is node of 'id' submesh is halo node of 'current_id' submesh
*/
bool decomp::isHalo(int x, int y, const vector<pair<size_t,vector<int>>>& submeshes, int Nx, int Ny, size_t current_id)
{
    if (decomp::getSubmeshIdByCoords(x,y,submeshes,Nx,Ny) == current_id)//not interface
    {
        if(!(decomp::getSubmeshIdByCoords(x + 1,y,submeshes,Nx,Ny) == current_id) || !(decomp::getSubmeshIdByCoords(x,y + 1,submeshes,Nx,Ny) == current_id)
                ||(!(decomp::getSubmeshIdByCoords(x + 1,y + 1,submeshes,Nx,Ny) == current_id)) || (!(decomp::getSubmeshIdByCoords(x - 1,y,submeshes,Nx,Ny) == current_id))
                || (!(decomp::getSubmeshIdByCoords(x,y - 1,submeshes,Nx,Ny) == current_id)) || (!(decomp::getSubmeshIdByCoords(x - 1,y - 1,submeshes,Nx,Ny) == current_id))
                || (!(decomp::getSubmeshIdByCoords(x - 1,y + 1,submeshes,Nx,Ny) == current_id)) || (!(decomp::getSubmeshIdByCoords(x + 1,y - 1,submeshes,Nx,Ny) == current_id)))
        {
            return true;
        }
    }
    return false;
}

/*
 * Global mesh nodes indexing
*/
void decomp::getGlobalIndexes(map<int,int>& G2L, vector<int>& global)
{
    global.clear();

    for(map<int,int>::iterator map_iter = G2L.begin(); map_iter != G2L.end(); ++map_iter)//forming global
    {
        global.push_back(map_iter->first);
    }
}

/*
 * Local mesh nodes indexing
*/
void decomp::getLocalIndexes(map<int,int>& G2L, vector<int>& local)
{
    local.clear();

    for(map<int,int>::iterator map_iter = G2L.begin(); map_iter != G2L.end(); ++map_iter)//forming local
    {
        local.push_back(map_iter->second);
    }
}
