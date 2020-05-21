#include "parallel.hpp"

void parallel::build_list_of_neighbors(map<int, int> &list_of_neighbors, const vector<int> &part, int self_id) {
    set<int> temp;
    int cur_idx = 0;

    for (size_t i = 0; i < part.size(); ++i) {
        temp.insert(part[i]);
    }

    temp.erase(self_id);

    for (set<int>::iterator it = temp.begin(); it != temp.end(); ++it) {
        list_of_neighbors[*it] = cur_idx;
        cur_idx++;
    }

    return;
}

void parallel::build_list_send_recv(VariableSizeMeshContainer<int> &topoNN_2, map<int, int> &G2L, vector<int> &L2G,
                                    vector<int> &part, map<int, int> &list_of_neighbors, vector<set<int>> &send,
                                    vector<set<int>> &recv, unsigned n_own, unsigned self_id) {

    unsigned key, position;

    vector<set<int>> recv_temp(list_of_neighbors.size());
    vector<set<int>> send_temp(list_of_neighbors.size());

    for (unsigned i = 0; i < topoNN_2.getBlockNumber(); ++i) {
        for (unsigned j = 0; j < topoNN_2.getBlockSize(i); ++j) {
            key = G2L[topoNN_2[i][j]];

            if (part[key] != self_id) {
                position = list_of_neighbors[part[key]];    // Порядковый номер соседа
                send_temp[position].insert(L2G[i]);         // send
                recv_temp[position].insert(topoNN_2[i][j]); // recv
            }
        }
    }

    send = send_temp;
    recv = recv_temp;

    return;
}
