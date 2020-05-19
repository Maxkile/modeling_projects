#include "decomposition.hpp"

/*
 * Mesh decomposing(only decart meshes are supported yet),returns pair - submesh
 * id -> submesh coords
 */
pair<size_t, vector<int>> decomp::decomposeMesh(int Nx, int Ny, int Px, int Py, int px, int py) {

    pair<size_t, vector<int>> submesh_params;
    size_t id;
    vector<int> coord(4); // ibeg, iend; jbeg, jend indexes

    if ((px >= Px) || (py >= Py) || (px < 0) || (py < 0)) {
        cerr << "Decomposition error: Wrong \'px\' or \'py\' values!" << endl;
        return submesh_params;
    } else if ((px >= Px) || (py >= Py)) {
        cerr << "Decomposition error: Part index is more than part number!" << endl;
        return submesh_params;
    } else if ((Px > Nx - 1) || (Py > Ny - 1)) {
        cerr << "Decomposition error: Number of parts is more than elements "
                "themselves!"
             << endl;
        return submesh_params;
    } else {
        int xPartSize = (Nx - 1) / Px; // Nx nodes -> (Nx - 1) elements on OX
        int yPartSize = (Ny - 1) / Py; // Ny nodes -> (Ny - 1) elements on OY

        int xLeft = (Nx - 1) % Px;
        int yLeft = (Ny - 1) % Py;

        if ((px < xLeft) && (py < yLeft)) // add left element to coord[x] and coord[y]
        {
            coord[0] = xPartSize * px + px;
            coord[1] = xPartSize * (px + 1) + px + 1;
            coord[2] = yPartSize * py + py;
            coord[3] = yPartSize * (py + 1) + py + 1;
        } else if (py < yLeft) // add left element to coord[y]
        {
            coord[0] = xPartSize * px + xLeft;
            coord[1] = xPartSize * (px + 1) + xLeft;
            coord[2] = yPartSize * py + py;
            coord[3] = yPartSize * (py + 1) + py + 1;
        } else if (px < xLeft) // add left element to coord[x]
        {
            coord[0] = xPartSize * px + px;
            coord[1] = xPartSize * (px + 1) + px + 1;
            coord[2] = yPartSize * py + yLeft;
            coord[3] = yPartSize * (py + 1) + yLeft;
        } else // no extra element added
        {
            coord[0] = xPartSize * px + xLeft;
            coord[1] = xPartSize * (px + 1) + xLeft;
            coord[2] = yPartSize * py + yLeft;
            coord[3] = yPartSize * (py + 1) + yLeft;
        }

        id = px * Py + py;
        return std::make_pair(id, coord);
    }
}

/*
 *  Mesh decomposing(only decart meshes are supported yet). Get all "ibeg, iend;
 * jbeg, jend" for all indexes [0,Px - 1],[0, Py - 1]
 */
vector<pair<size_t, vector<int>>> decomp::decomposeMesh(int Nx, int Ny, int Px, int Py) {

    vector<pair<size_t, vector<int>>> submeshes;

    if ((Px > Nx - 1) || (Py > Ny - 1)) {
        cerr << "Decomposition error: Number of parts is more than elements "
                "themselves!"
             << endl;
        return submeshes;
    } else {
        for (int px = 0; px < Px; ++px) {
            for (int py = 0; py < Py; ++py) {
                submeshes.emplace_back(decomposeMesh(Nx, Ny, Px, Py, px, py));
            }
        }
        return submeshes;
    }
}

/*
 * Get submmesh global id among submeshes
 */
size_t decomp::getSubmeshIdByCoords(int x, int y, const vector<pair<size_t, vector<int>>> &submeshes, int Nx, int Ny) {
    size_t submesh_id;
    for (auto iter = submeshes.begin(); iter != submeshes.end(); ++iter) {
        int beg_x = iter->second[0];
        int end_x = iter->second[1];
        int beg_y = iter->second[2];
        int end_y = iter->second[3];

        if ((end_x == Nx - 1) && (end_y == Ny - 1)) // border
        {
            if ((x >= beg_x) && (x <= end_x) && (y >= beg_y) && (y <= end_y)) {
                submesh_id = iter->first;
                break;
            }
        } else if (end_x == Nx - 1) // border
        {
            if ((x >= beg_x) && (x <= end_x) && (y >= beg_y) && (y < end_y)) {
                submesh_id = iter->first;
                break;
            }
        } else if (end_y == Ny - 1) // border
        {
            if ((x >= beg_x) && (x < end_x) && (y >= beg_y) && (y <= end_y)) {
                submesh_id = iter->first;
                break;
            }
        } else {
            if ((x >= beg_x) && (x < end_x) && (y >= beg_y) && (y < end_y)) {
                submesh_id = iter->first;
                break;
            }
        }
    }

    return submesh_id;
}

/*
 * Is node of 'id' submesh is just interface neighbour node of 'current_id'
 * submesh
 *
 */

bool decomp::isInterfaceNeighbour(int x, int y, const vector<pair<size_t, vector<int>>> &submeshes, int Nx, int Ny,
                                  size_t current_id) {
    if (x >= 0 && y >= 0 && x < Nx && y < Ny &&
        !(decomp::getSubmeshIdByCoords(x, y, submeshes, Nx, Ny) == current_id)) {
        return true;
    } else {
        return false;
    }
}

/*
 * Is node of 'id' submesh is halo node of 'current_id' submesh
 */
bool decomp::isInterface(int x, int y, const vector<pair<size_t, vector<int>>> &submeshes, int Nx, int Ny, int k3,
                         int k4, size_t current_id) {
    if (!isInterfaceNeighbour(x, y, submeshes, Nx, Ny, current_id) ||
        !(x >= 0 && y >= 0 && x < Nx && y < Ny)) // not neigthbour
    {
        if (isInterfaceNeighbour(x + 1, y, submeshes, Nx, Ny, current_id) ||
            isInterfaceNeighbour(x, y + 1, submeshes, Nx, Ny, current_id) ||
            isInterfaceNeighbour(x - 1, y, submeshes, Nx, Ny, current_id) ||
            isInterfaceNeighbour(x, y - 1, submeshes, Nx, Ny, current_id)) {
            return true;
        } else {
            // Separately checking nodes (x + 1,y - 1) and (x - 1,y + 1)
            int figuresLeft;
            // checking nodes (x,y - 1) and (x - 1,y) as if they are triangle edges
            if (isInterfaceNeighbour(x - 1, y + 1, submeshes, Nx, Ny, current_id)) {
                figuresLeft = ((x - 1) * Ny + y) % (k3 + k4);
            } else if (isInterfaceNeighbour(x + 1, y - 1, submeshes, Nx, Ny, current_id)) {
                figuresLeft = (x * Ny + (y - 1)) % (k3 + k4);
            }
            return (figuresLeft >= 0 && figuresLeft < k4);
        }
    }
    return false;
}

/*
 * Considering that (cur_i,cur_j) node is an interface node check and add all
 * neighbour nodes to halo set. Also adds them to part set. Mesh is
 * cartesian-like
 *
 * - - - - -
 * - # # # -
 * - # * # -
 * - # # # -
 * - - - - -
 *
 * For '*' node '#' will be checked if they are halo or not.
 */
void decomp::addHaloNodes(int x, int y, const vector<pair<size_t, vector<int>>> &submeshes, int Nx, int Ny, int k3,
                          int k4, size_t current_id, set<int> &haloes) {
    for (int i = -1; i <= 1; ++i) {
        for (int j = -1; j <= 1; ++j) {
            int newX = x + i;
            int newY = y + j;
            // Skipping these nodes...
            if ((newX == x + 1 && newY == y + 1) || (newX == x - 1 && newY == y - 1)) {
                continue;
            } else if (newX >= 0 && newY >= 0 && newX < Nx && newY < Ny &&
                       isInterfaceNeighbour(newX, newY, submeshes, Nx, Ny, current_id)) {
                int figuresLeft;

                // Separately checking nodes (x + 1,y - 1) and (x - 1,y + 1)
                // checking nodes (x,y - 1) and (x - 1,y) as if they are triangle edges
                if (newX == x - 1 && newY == y + 1) {
                    figuresLeft = ((x - 1) * Ny + y) % (k3 + k4);
                } else if (newX == x + 1 && newY == y - 1) {
                    figuresLeft = (x * Ny + (y - 1)) % (k3 + k4);
                } else {
                    haloes.insert(newX * Ny + newY);
                    continue;
                }
                if (figuresLeft >= 0 && figuresLeft < k4) {
                    haloes.insert(newX * Ny + newY);
                }
            }
        }
    }
}

/*
 * Forming part vector from vector of nodes
 */
void decomp::formPart(vector<int> &part, const vector<int> &nodes, const vector<pair<size_t, vector<int>>> &submeshes,
                      int Nx, int Ny) {
    part.clear();
    size_t size = nodes.size();

    part.resize(size);

    for (size_t i = 0; i < size; ++i) {
        int x = nodes[i] / Ny;
        int y = nodes[i] % Ny;
        part[i] = decomp::getSubmeshIdByCoords(x, y, submeshes, Nx, Ny);
    }
}

/*
 * Global mesh nodes indexing
 */
void decomp::getGlobalIndexes(map<int, int> &G2L, vector<int> &global) {
    global.clear();

    for (map<int, int>::iterator map_iter = G2L.begin(); map_iter != G2L.end(); ++map_iter) // forming global
    {
        global.push_back(map_iter->first);
    }
}

/*
 * Local mesh nodes indexing
 */
void decomp::getLocalIndexes(map<int, int> &G2L, vector<int> &local) {
    local.clear();

    for (map<int, int>::iterator map_iter = G2L.begin(); map_iter != G2L.end(); ++map_iter) // forming local
    {
        local.push_back(map_iter->second);
    }
}
