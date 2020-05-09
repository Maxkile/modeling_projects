#include "toposBuild.hpp"

#include <omp.h>

// coords
void topos::build_coord(FixedSizeMeshContainer<double> &C, int Lx, int Ly,
                        int Nx, int Ny) {
  vector<double> temp;
  temp.reserve(2 * Nx * Ny);

  for (int i = 0; i < Ny; i++) {
    for (int j = 0; j < Nx; j++) {
      temp.push_back((Lx / static_cast<double>(Nx - 1)) * i);
      temp.push_back((Ly / static_cast<double>(Ny - 1)) * j);
    }
  }
  C.add(temp);
}

// defines how many triangles and squares left in mesh after skipping elements
int topos::computeMeshFiguresNumberLeft(int figCount1, int figCount2,
                                        int skippedElemsCount,
                                        int curMeshFigureStructure) {
  int elemsLeft = skippedElemsCount > curMeshFigureStructure
                      ? (figCount1 + figCount2) -
                            ((skippedElemsCount - curMeshFigureStructure) %
                             (figCount1 + figCount2))
                      : (figCount1 + figCount2) - curMeshFigureStructure;

  // curMeshFigureStructure - left figures to put on mesh till current time
  // consider curMeshFigureStructure is ALWAYS <= figCount1 + figCount2

  if (skippedElemsCount > curMeshFigureStructure) {
    elemsLeft = (figCount1 + figCount2) -
                ((skippedElemsCount - curMeshFigureStructure) %
                 (figCount1 + figCount2));
  } else if (skippedElemsCount == curMeshFigureStructure) {
    elemsLeft =
        (figCount1 + figCount2) - (skippedElemsCount - curMeshFigureStructure);
  } else {
    elemsLeft = curMeshFigureStructure - skippedElemsCount;
  }

  return elemsLeft;
}

VariableSizeMeshContainer<int>
topos::toLocalIndexes(const VariableSizeMeshContainer<int> &originEN,
                      map<int, int> &G2L) {
  vector<int> BlockSize;
  vector<int> temp;
  VariableSizeMeshContainer<int> local(temp, BlockSize);

  for (size_t i = 0; i < originEN.getBlockNumber(); ++i) {
    size_t blockSize = originEN.getBlockSize(i);

    for (size_t j = 0; j < blockSize; ++j) {
      temp.push_back(G2L[originEN[i][j]]);
    }
    BlockSize.push_back(blockSize);
  }

  local.add(temp, BlockSize);

  return local;
}

vector<int> topos::toLocalIndexes(const vector<int> &origin,
                                  map<int, int> &G2L) {
  vector<int> local(origin.size());
  for (size_t i = 0; i < origin.size(); ++i) {
    local[i] = G2L[origin[i]];
  }
  return local;
}

/*
 * Also generates G2L,L2G,halo's,interfaces vectors required for topologies
 * generating.
 * @nodes Global numeration for inner, halo and interface submesh nodes
 */
VariableSizeMeshContainer<int>
topos::build_topoEN(int Nx, int Ny, int k3, int k4, size_t submesh_id,
                    const vector<pair<size_t, vector<int>>> &submeshes,
                    map<int, int> &G2L, vector<int> &L2G, vector<int> &inner,
                    vector<int> &interface, vector<int> &haloes) {
  G2L.clear();
  L2G.clear();

  inner.clear();
  interface.clear();
  haloes.clear();

  vector<int> BlockSize;
  vector<int> temp;
  VariableSizeMeshContainer<int> topoEN(temp, BlockSize);

  int beg_i = submeshes[submesh_id].second[0];
  int end_i = submeshes[submesh_id].second[1];
  int beg_j = submeshes[submesh_id].second[2];
  int end_j = submeshes[submesh_id].second[3];

  set<int> haloes_set;

  if ((beg_i > Nx) || (end_i > Nx) || (beg_j > Ny) || (end_j > Ny) ||
      (end_i <= 0) || (end_i <= 0) || (end_i <= 0) || (end_i <= 0)) {
    cerr << "Wrong submesh parameters!" << endl;
    return topoEN;
  } else if (submesh_id >= submeshes.size()) {
    cerr << "Wrong submesh id!" << endl;
    return topoEN;
  }

  int fullElementsSkipped =
      (Ny - 1) * beg_i + beg_j; // total number of elements is k3 + k4. So two
                                // triangles is one element itself
  int meshFigureStructureCur =
      computeMeshFiguresNumberLeft(k3, k4, fullElementsSkipped, 0);

  int local_j = 0; // local numeration
  int cur_i = beg_i;
  int cur_j = beg_j;

  while (cur_i < end_i) {
    while (cur_j < end_j) {
      G2L.insert(pair<int, int>(Ny * cur_i + cur_j, local_j));
      L2G.push_back(Ny * cur_i + cur_j);

      if (decomp::isHalo(cur_i, cur_j, submeshes, Nx, Ny, submesh_id)) {
        haloes_set.insert(Ny * cur_i + cur_j);
      } else if (decomp::isInterface(cur_i, cur_j, submeshes, Nx, Ny,
                                     submesh_id)) {
        interface.push_back(Ny * cur_i + cur_j);
        decomp::addHaloNodes(cur_i, cur_j, submeshes, Nx, Ny, submesh_id,
                             haloes_set);
      } else {
        inner.push_back(Ny * cur_i + cur_j);
      }

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

      cur_j++;
      local_j++;
    }

    fullElementsSkipped = (Ny - (cur_j + 1)) + beg_j;
    meshFigureStructureCur = computeMeshFiguresNumberLeft(
        k3, k4, fullElementsSkipped, meshFigureStructureCur);

    // changing y coord
    G2L.insert(pair<int, int>(Ny * cur_i + cur_j, local_j));
    L2G.push_back(Ny * cur_i + cur_j);

    if (decomp::isHalo(cur_i, cur_j, submeshes, Nx, Ny, submesh_id)) {
      haloes_set.insert(Ny * cur_i + cur_j);
    } else if (decomp::isInterface(cur_i, cur_j, submeshes, Nx, Ny,
                                   submesh_id)) {
      interface.push_back(Ny * cur_i + cur_j);
      decomp::addHaloNodes(cur_i, cur_j, submeshes, Nx, Ny, submesh_id,
                           haloes_set);
    } else {
      inner.push_back(Ny * cur_i + cur_j);
    }

    local_j++;
    cur_j = beg_j;
    cur_i++;
  }

  for (cur_j = beg_j; cur_j <= end_j;
       ++cur_j, ++local_j) // last y, we haven't visited it yet, but have to
                           // store them too
  {
    G2L.insert(pair<int, int>(Ny * end_i + cur_j, local_j));
    L2G.push_back(Ny * end_i + cur_j);

    if (decomp::isHalo(cur_i, cur_j, submeshes, Nx, Ny, submesh_id)) {
      haloes_set.insert(Ny * cur_i + cur_j);
    } else if (decomp::isInterface(cur_i, cur_j, submeshes, Nx, Ny,
                                   submesh_id)) // if interface -> check and add
                                                // all border nodes as haloes
    {
      interface.push_back(Ny * cur_i + cur_j);
      decomp::addHaloNodes(cur_i, cur_j, submeshes, Nx, Ny, submesh_id,
                           haloes_set);
    } else {
      inner.push_back(Ny * cur_i + cur_j);
    }
  }

  haloes.reserve(haloes_set.size());
  for (auto it = haloes_set.cbegin(); it != haloes_set.cend(); ++it) {
    haloes.push_back(*it);
  }

  topoEN.add(temp, BlockSize);
  return topoEN;
}

// topoSN
VariableSizeMeshContainer<int> topos::build_topoSN(int Nx, int Ny, int k3,
                                                   int k4) {
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
VariableSizeMeshContainer<int>
topos::build_reverse_topo(const VariableSizeMeshContainer<int> &topo) {
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
      reverse_topo[topo[i][j]]
                  [count_mass[topo[i][j]] - count_i_mas[topo[i][j]]] = i;
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
VariableSizeMeshContainer<int>
topos::build_topoBNS(const VariableSizeMeshContainer<int> &topo) {
  vector<int> BlockSize;
  vector<int> temp;
  int Nx, Ny;
  VariableSizeMeshContainer<int> topoBNS(temp, BlockSize);
  for (Nx = 1; topo[Nx - 1][0] == 0; ++Nx) {
  }

  for (Ny = 1; topo[Nx + Ny - 2][0] == 1; ++Ny) {
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

// topoNN
VariableSizeMeshContainer<int>
topos::build_topoNN(const VariableSizeMeshContainer<int> &topoSN) {
  vector<int> BlockSize;
  vector<int> temp;
  vector<int> count_i_mas;
  vector<int> count_mass;
  int k;

  VariableSizeMeshContainer<int> topoNN(temp, BlockSize);
  int nS = topoSN.getBlockNumber();

  int nN = 0;
  for (int i = 0; i < nS; i++) {
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
