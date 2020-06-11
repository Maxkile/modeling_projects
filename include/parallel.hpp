#ifndef PARALLEL_HPP
#define PARALLEL_HPP

#pragma once

#include "VariableSizeMeshContainer.hpp"
#include "stdafx.hpp"
#include "toposBuild.hpp"
#include <algorithm>
#include <cstdarg>
#include <map>
#include <mpi.h>
#include <set>

using namespace std;

namespace parallel {

struct Ne_scheme_bufs {
    size_t neighbour_id;

    double *send_buf;
    double *recv_buf;


    Ne_scheme_bufs() :  send_buf(nullptr), recv_buf(nullptr) {}

    ~Ne_scheme_bufs();
};

struct Decision {
    int id;
    double answer;
};

#define crash(...) exit(Crash(__VA_ARGS__)) // via exit define so static analyzer knows its an exit point
int Crash(const char *fmt, ...);

int pprintf(const char *fmt, ...);

int printf_master(int id, const char *log, ...);

// Internal interface to MPI barrier
static inline void barrier() {
    if (MPI_Barrier(MPI_COMM_WORLD) != MPI_SUCCESS)
        crash("Base lib barrier: MPI_Barrier failed! \n");
}

void build_list_of_neighbors(map<int, int> &list_of_neighbors, const vector<int> &part, int self_id);

void build_list_send_recv(VariableSizeMeshContainer<int> &topoNN, vector<int> &part, map<int, int> &list_of_neighbors,
                          vector<set<int>> &send, vector<set<int>> &recv, int self_id);

void gather_all(Decision *total, const vector<double> &local_solution, size_t n_own, const vector<int> &L2G);

void update_halo(vector<double> &x, map<int, int> &list_of_neighbors, vector<set<int>> &send, vector<set<int>> &recv,
                 int proccessor_id, MPI_Comm mpi_comm = MPI_COMM_WORLD);

} // namespace parallel

#endif // PARALLEL_HPP
