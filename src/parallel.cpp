#include "parallel.hpp"

parallel::Ne_scheme_bufs::~Ne_scheme_bufs() {
    delete[] recv_buf;
    delete[] send_buf;
}

void parallel::build_list_of_neighbors(map<int, int> &list_of_neighbors, const vector<int> &part, int self_id) {
    set<int> temp;
    list_of_neighbors.clear();

    int cur_idx = 0;

    for (size_t i = 0; i < part.size(); ++i) {
        temp.insert(part[i]);
    }

    temp.erase(self_id);

    for (set<int>::iterator it = temp.begin(); it != temp.end(); ++it) {
        list_of_neighbors[*it] = cur_idx;
        cur_idx++;
    }
}

void parallel::build_list_send_recv(VariableSizeMeshContainer<int> &topoNN, map<int, int> &G2L, vector<int> &L2G,
                                    vector<int> &part, map<int, int> &list_of_neighbors, vector<set<int>> &send,
                                    vector<set<int>> &recv, size_t n_own, int self_id) {

    size_t key, position;
    recv.clear();
    send.clear();

    recv.resize(list_of_neighbors.size());
    send.resize(list_of_neighbors.size());

    for (size_t i = 0; i < topoNN.getBlockNumber(); ++i) {
        for (size_t j = 0; j < topoNN.getBlockSize(i); ++j) {
            key = G2L[topoNN[i][j]];

            if (part[key] != self_id) {
                position = list_of_neighbors[part[key]]; // neighbourn number
                send[position].insert(L2G[i]);           // send
                recv[position].insert(topoNN[i][j]);     // recv
            }
        }
    }
}

// Considering 'send','recv' lists are in local numeration(if not -> G2L needed)
void parallel::update_halo(vector<double> &nodes_values, size_t n_own, map<int, int> &list_of_neighbors,
                           vector<set<int>> &send, vector<set<int>> &recv, int processor_id, MPI_Comm mpi_comm) {

    size_t scheme_size = list_of_neighbors.size();
    static vector<Ne_scheme_bufs> scheme(scheme_size);

    size_t j = 0, i = 0;
    size_t send_size, recv_size;
    size_t request_size = 2 * scheme_size;

    MPI_Request requests[request_size];
    MPI_Status statuses[request_size];

    int mpi_res;
    for (auto it = list_of_neighbors.cbegin(); it != list_of_neighbors.cend(); ++it, ++j) {
        send_size = send[j].size();
        recv_size = recv[j].size();

        scheme[j].neighbour_id = it->first; // for whom
        if (scheme[j].send_buf) {
            scheme[j].send_buf = new double(send_size); // interface for neighbour
        }
        if (scheme[j].recv_buf) {
            scheme[j].recv_buf = new double(recv_size); // halo for us
        }
        i = 0;
        for (auto send_it = send[j].cbegin(); send_it != send[j].cend(); ++send_it, ++i) {
            scheme[j].send_buf[i] = nodes_values[*send_it];
        }
        i = 0;
        for (auto recv_it = send[j].cbegin(); recv_it != recv[j].cend(); ++recv_it, ++i) {
            scheme[j].recv_buf[i] = nodes_values[*recv_it];
        }

        mpi_res =
            MPI_Isend(scheme[j].send_buf, send_size, MPI_DOUBLE, scheme[j].neighbour_id, j, mpi_comm, &requests[j]);

        if (mpi_res != MPI_SUCCESS) {
            cerr << "MPI Isend error: processor " << processor_id << " crashed! Exiting..." << endl;
            exit(1);
        }

        mpi_res =
            MPI_Irecv(scheme[j].recv_buf, recv_size, MPI_DOUBLE, scheme[j].neighbour_id, j, mpi_comm, &requests[j + 1]);

        if (mpi_res != MPI_SUCCESS) {
            cerr << "MPI Irevc error: processor " << processor_id << " crashed! Exiting..." << endl;
            exit(1);
        }
    }

    mpi_res = MPI_Waitall(request_size, requests, statuses); // blocking
    if (mpi_res != MPI_SUCCESS) {
        cerr << "MPI Waitall error: processor " << processor_id << " crashed! Exiting..." << endl;
        exit(1);
    }

    // considering halo nodes in x are sorted by owners order
    size_t offset = 0;
    for (auto it = list_of_neighbors.cbegin(); it != list_of_neighbors.cend(); ++it, ++j) {
        recv_size = recv[j].size();
        for (size_t k = 0; k < recv_size; ++k) {
            nodes_values[n_own + offset + k] = scheme[j].recv_buf[k];
        }
        offset += recv_size;
    }
}

int parallel::printf_master(int id, const char *fmt, ...) { // Write to stdout from Master process
    int r = 0;
    va_list ap;
    if (id == _MAIN_ID) {
        va_start(ap, fmt);
        r = vfprintf(stdout, fmt, ap);
        va_end(ap);
    }
    fflush(stdout);
    return r;
}
