#include "toposBuild.hpp"

#include <omp.h>

// coords
void topos::build_coord(FixedSizeMeshContainer<double> &C, int Lx, int Ly, int Nx, int Ny) {
    vector<double> temp;
    temp.reserve(2 * static_cast<size_t>(Nx) * Ny);

    for (int i = 0; i < Ny; i++) {
        for (int j = 0; j < Nx; j++) {
            temp.push_back((Lx / static_cast<double>(Nx) - 1) * i);
            temp.push_back((Ly / static_cast<double>(Ny) - 1) * j);
        }
    }
    C.add(temp);
}

/*
 * Defines how many triangles and squares left in mesh after skipping elements
 */
int topos::computeMeshFiguresNumberLeft(int figCount1, int figCount2, int skippedElemsCount,
                                        int curMeshFigureStructure) {
    int elemsLeft =
        skippedElemsCount > curMeshFigureStructure
            ? (figCount1 + figCount2) - ((skippedElemsCount - curMeshFigureStructure) % (figCount1 + figCount2))
            : (figCount1 + figCount2) - curMeshFigureStructure;

    // curMeshFigureStructure - left figures to put on mesh till current time
    // consider curMeshFigureStructure is ALWAYS <= figCount1 + figCount2

    if (skippedElemsCount > curMeshFigureStructure) {
        elemsLeft = (figCount1 + figCount2) - ((skippedElemsCount - curMeshFigureStructure) % (figCount1 + figCount2));
    } else if (skippedElemsCount == curMeshFigureStructure) {
        elemsLeft = (figCount1 + figCount2) - (skippedElemsCount - curMeshFigureStructure);
    } else {
        elemsLeft = curMeshFigureStructure - skippedElemsCount;
    }

    return elemsLeft;
}

VariableSizeMeshContainer<int> topos::toLocalIndexes(const VariableSizeMeshContainer<int> &originEN,
                                                     const mapping &G2L) {
    vector<int> BlockSize;
    vector<int> temp;
    VariableSizeMeshContainer<int> local(temp, BlockSize);

    for (size_t i = 0; i < originEN.getBlockNumber(); ++i) {
        size_t blockSize = originEN.getBlockSize(i);

        for (size_t j = 0; j < blockSize; ++j) {
            temp.push_back(G2L.find(originEN[i][j])->second);
        }
        BlockSize.push_back(blockSize);
    }

    local.add(temp, BlockSize);

    return local;
}

VariableSizeMeshContainer<int> topos::toGlobalIndexes(const VariableSizeMeshContainer<int> &topoNN,
                                                      const vector<int> &L2G, size_t n_own) {
    vector<int> BlockSize;
    vector<int> temp;
    VariableSizeMeshContainer<int> topoNN_2(temp, BlockSize);

    for (size_t i = 0; i < n_own; ++i) {
        size_t blockSize = topoNN.getBlockSize(i);

        for (size_t j = 0; j < blockSize; ++j) {
            temp.push_back(L2G[topoNN[i][j]]);
        }

        BlockSize.push_back(blockSize);
    }

    topoNN_2.add(temp, BlockSize);

    return topoNN_2;
}

vector<int> topos::toLocalIndexes(const vector<int> &origin, const mapping &G2L) {
    vector<int> local(origin.size());
    for (size_t i = 0; i < origin.size(); ++i) {
        local[i] = G2L.find(origin[i])->second;
    }
    return local;
}

vector<int> topos::toGlobalIndexes(const vector<int> &origin, const vector<int> &L2G, size_t n_own) {
    size_t size = origin.size();
    vector<int> local(size);
    for (size_t i = 0; i < size; ++i) {
        local[i] = origin[L2G[i]];
    }
    return local;
}

/*
 * Also generates G2L,L2G,halo's,interfaces vectors required for topologies
 * generating.
 * @nodes Global numeration for inner, halo and interface submesh nodes
 */
VariableSizeMeshContainer<int> topos::build_topoEN(int Nx, int Ny, int k3, int k4, size_t submesh_id,
                                                   const vector<pair<size_t, vector<int>>> &submeshes, mapping &G2L,
                                                   vector<int> &L2G, vector<int> &nodes, vector<int> &part,
                                                   size_t &n_own) {
    G2L.clear();
    L2G.clear();
    nodes.clear();
    part.clear();

    vector<int> inner;
    vector<int> interface;

    vector<pair<size_t, vector<int>>> haloes;

    // for sorting halo by owners
    vector<int> haloes_left;
    vector<int> haloes_up;
    vector<int> haloes_right;
    vector<int> haloes_down;

    vector<int> BlockSize;
    vector<int> temp;
    VariableSizeMeshContainer<int> topoEN(temp, BlockSize);

    size_t beg_i = submeshes[submesh_id].second[0];
    size_t end_i = submeshes[submesh_id].second[1];
    size_t beg_j = submeshes[submesh_id].second[2];
    size_t end_j = submeshes[submesh_id].second[3];

    if ((beg_i > Nx) || (end_i > Nx) || (beg_j > Ny) || (end_j > Ny) || (beg_i < 0) || (beg_j < 0) || (end_i <= 0) ||
        (end_j <= 0)) {
        cerr << "Wrong submesh parameters!" << endl;
        return topoEN;
    }

    int fullElementsSkipped = (Ny - 1) * beg_i + beg_j; // total number of elements is k3 + k4. So two
                                                        // triangles are one element themselves
    int meshFigureStructureCur = computeMeshFiguresNumberLeft(k3, k4, fullElementsSkipped, 0);

    size_t cur_i = beg_i, cur_j = beg_j, node_id = 0, haloes_size = 0;
    while (cur_i <= end_i) {
        while (cur_j <= end_j) {
            // Forming vectors
            if (decomp::isHalo(cur_i, cur_j, submeshes, submesh_id, Nx, Ny)) {
                decomp::insertHalo(haloes, Ny * cur_i + cur_j, node_id);
                haloes_size++;
            } else if (decomp::isInterface(cur_i, cur_j, submeshes, submesh_id, Nx, Ny)) {
                interface.push_back(Ny * cur_i + cur_j);
            } else {
                inner.push_back(Ny * cur_i + cur_j); // inner
            }

            if (cur_j == end_j || cur_i == end_i) { // on the border, skip building topology element(on border)
                cur_j++;
                continue;
            } else {

                // Building topology
                if (meshFigureStructureCur > k4) // triangle
                {
                    temp.push_back(Ny * cur_i + cur_j);
                    temp.push_back(Ny * cur_i + cur_j + 1);
                    temp.push_back(Ny * (cur_i + 1) + cur_j);
                    BlockSize.push_back(TRIANGLE_NODES);

                    temp.push_back(Ny * cur_i + cur_j + 1);
                    temp.push_back(Ny * (cur_i + 1) + cur_j + 1);
                    temp.push_back(Ny * (cur_i + 1) + cur_j);
                    BlockSize.push_back(TRIANGLE_NODES);
                } else if (meshFigureStructureCur <= k4) // square
                {
                    temp.push_back(Ny * cur_i + cur_j);
                    temp.push_back(Ny * cur_i + cur_j + 1);
                    temp.push_back(Ny * (cur_i + 1) + cur_j + 1);
                    temp.push_back(Ny * (cur_i + 1) + cur_j);

                    BlockSize.push_back(SQUARE_NODES);
                }

                meshFigureStructureCur--;
                if (meshFigureStructureCur == 0) {
                    meshFigureStructureCur = k3 + k4;
                }
            }
            cur_j++;
        }

        fullElementsSkipped = (Ny - cur_j) + beg_j;
        meshFigureStructureCur = computeMeshFiguresNumberLeft(k3, k4, fullElementsSkipped, meshFigureStructureCur);

        cur_j = beg_j;
        cur_i++;
    }

    // Forming own nodes number(n_own)
    n_own = inner.size() + interface.size();
    size_t n_all = n_own + haloes_size;

    // Forming nodes vector
    nodes.reserve(n_all);
    vmo::join(nodes, inner, interface);
    for (auto iter = haloes.cbegin(); iter != haloes.cend(); ++iter) {
        vmo::join(nodes, iter->second);
    }

    // Forming G2L,L2G and part
    L2G.reserve(n_all);
    part.reserve(n_all);
    for (size_t i = 0; i < nodes.size(); ++i) {
        int item = nodes[i];
        G2L.insert(pair<int, int>(item, i));
        L2G.push_back(item);
        part.push_back(decomp::getSubmeshIdByCoords(item / Ny, item % Ny, submeshes, Nx, Ny));
    }

    topoEN.add(temp, BlockSize);
    return topoEN;
}

// topoSN
VariableSizeMeshContainer<int> topos::build_topoSN(int Nx, int Ny, int k3, int k4) {
    vector<int> BlockSize;
    vector<int> temp;
    int k = 0;
    VariableSizeMeshContainer<int> topoSN(temp, BlockSize);

    for (int i = 0; i < (Nx * Ny); ++i) {
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

        if ((i % Nx != Nx - 1) && (i >= Nx)) {
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

// reversed topology
VariableSizeMeshContainer<int> topos::build_reverse_topo(const VariableSizeMeshContainer<int> &topo) {
    vector<int> BlockSize;
    vector<int> temp;
    vector<int> nt_i_mas;
    vector<int> count_mass;
    vector<int> count_i_mas;

    VariableSizeMeshContainer<int> reverse_topo(temp, BlockSize);
    size_t nE = topo.getBlockNumber();

    int nN = 0;
    for (size_t i = 0; i < nE; i++) {
        for (size_t j = 0; j < topo.getBlockSize(i); j++)
            if (topo[i][j] > nN)
                nN = topo[i][j];
    }
    nN++;

    count_i_mas.reserve(nN);
    count_mass.reserve(nN);

    for (int i = 0; i < nN; i++) {
        count_mass.push_back(0);
        count_i_mas.push_back(0);
    }

    for (size_t i = 0; i < nE; i++) {
        for (size_t j = 0; j < topo.getBlockSize(i); j++) {
            count_mass[topo[i][j]]++;
            count_i_mas[topo[i][j]]++;
        }
    }
    for (int i = 0; i < nN; i++) {
        BlockSize.push_back(count_mass[i]);
        for (int j = 0; j < count_mass[i]; j++)
            temp.push_back(0);
    }

    reverse_topo.add(temp, BlockSize);

    for (size_t i = 0; i < nE; i++) {
        for (size_t j = 0; j < topo.getBlockSize(i); j++) {
            reverse_topo[topo[i][j]][count_mass[topo[i][j]] - count_i_mas[topo[i][j]]] = i;
            count_i_mas[topo[i][j]]--;
        }
    }

    return reverse_topo;
}

// topoBSN
VariableSizeMeshContainer<int> topos::build_topoBSN(int Nx, int Ny) {
    vector<int> BlockSize;
    vector<int> temp;
    VariableSizeMeshContainer<int> topoBSN(temp, BlockSize);

    for (int i = 0; i < Nx - 1; ++i) {
        temp.push_back(0);
        temp.push_back(i);
        temp.push_back(i + 1);
        BlockSize.push_back(3);

        topoBSN.add(temp, BlockSize);

        temp.clear();
        BlockSize.clear();
    }

    for (int i = 0; i < Ny - 1; ++i) {
        temp.push_back(1);
        temp.push_back(Nx - 1 + i);
        temp.push_back(Nx - 1 + i + 1);
        BlockSize.push_back(3);
    }

    for (int i = 0; i < Nx - 1; ++i) {
        temp.push_back(2);
        temp.push_back(Nx + Ny - 2 + i);
        temp.push_back(Nx + Ny - 2 + i + 1);
        BlockSize.push_back(3);
    }

    for (int i = 0; i < Ny - 1; ++i) {
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
VariableSizeMeshContainer<int> topos::build_topoBNS(const VariableSizeMeshContainer<int> &topo) {
    vector<int> BlockSize;
    vector<int> temp;
    int Nx, Ny;
    VariableSizeMeshContainer<int> topoBNS(temp, BlockSize);
    for (Nx = 1; topo[static_cast<size_t>(Nx) - 1][0] == 0; ++Nx) {
    }

    for (Ny = 1; topo[static_cast<size_t>(Nx) + Ny - 2][0] == 1; ++Ny) {
    }

    for (int i = 0; i < (Nx + Ny - 2) * 2; ++i) {
        if (i == 0)
            temp.push_back(2 * (Nx + Ny - 2) - 1);
        else
            temp.push_back(i - 1);

        temp.push_back(i);
        BlockSize.push_back(2);
    }

    topoBNS.add(temp, BlockSize);

    return topoBNS;
}

// topoNN_1
VariableSizeMeshContainer<int> topos::build_topoNN_from_topoSN(const VariableSizeMeshContainer<int> &topoSN) {
    vector<int> BlockSize;
    vector<int> temp;
    vector<int> count_i_mas;
    vector<int> count_mass;
    int k;

    VariableSizeMeshContainer<int> topoNN(temp, BlockSize);
    size_t nS = topoSN.getBlockNumber();

    int nN = 0;
    for (size_t i = 0; i < nS; i++) {
        for (size_t j = 0; j < topoSN.getBlockSize(i); j++)
            if (topoSN[i][j] > nN)
                nN = topoSN[i][j];
    }
    nN++;

    count_i_mas.reserve(nN);
    count_mass.reserve(nN);

    for (int i = 0; i < nN; i++) {
        count_mass.push_back(0);
        count_i_mas.push_back(0);
    }

    for (int i = 0; i < nS; i++) {
        for (size_t j = 0; j < topoSN.getBlockSize(i); j++) {
            count_mass[topoSN[i][j]]++;
            count_i_mas[topoSN[i][j]]++;
        }
    }

    for (int i = 0; i < nN; i++) {
        BlockSize.push_back(count_mass[i]);
        for (int j = 0; j < count_mass[i]; j++)
            temp.push_back(0);
    }

    topoNN.add(temp, BlockSize);

    for (int i = 0; i < nS; i++) {
        for (size_t j = 0; j < topoSN.getBlockSize(i); j++) {
            k = count_mass[topoSN[i][j]] - count_i_mas[topoSN[i][j]];
            topoNN[topoSN[i][j]][k] = j ? topoSN[i][0] : topoSN[i][1];
            count_i_mas[topoSN[i][j]]--;
        }
    }

    return topoNN;
}

// topoNN_2
VariableSizeMeshContainer<int> topos::build_topoNN_from_topoEN(const VariableSizeMeshContainer<int> &topoEN) {
    int k, max_node = 0, t;

    k = topoEN.getBlockNumber();
    // max_node = topoEN[k - 1][topoEN.getBlockSize(k - 1) - 2] + 1;

    for (int i = 0; i < k; i++) {
        t = topoEN.getBlockSize(i);

        for (int j = 0; j < t; j++) {
            if (max_node < topoEN[i][j])
                max_node = topoEN[i][j];
        }
    }

    max_node++;

    vector<map<int, int>> arr(max_node);

    for (int i = 0; i < k; i++) {
        t = topoEN.getBlockSize(i);

        for (int j = 0; j < t; j++) {
            arr[topoEN[i][j]][topoEN[i][(j + 1) % t]] = 0;
            arr[topoEN[i][j]][topoEN[i][(j - 1 + t) % t]] = 0;
        }
    }

    vector<int> BlockSize(max_node, 0);
    vector<int> temp;

    for (int i = 0; i < max_node; i++) {
        for (map<int, int>::iterator it = arr[i].begin(); it != arr[i].end(); it++) {
            ++BlockSize[i];
            temp.push_back(it->first);
        }
    }

    VariableSizeMeshContainer<int> topoNN(temp, BlockSize);

    return topoNN;
}
