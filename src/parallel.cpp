#include "parallel.hpp"

int parallel::Crash(const char *fmt, ...) { // termination of program due to error
    va_list ap;
    fprintf(stderr, "\nEpic fail: MyID = %d\n", proc_id);

    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);

    fprintf(stderr, "\n");
    fflush(stderr);

    MPI_Abort(MCW, -1);

    return 0;
}

// Debug synchronous printf - Write to stdout + flush + barrier
int parallel::pprintf(const char *fmt, ...) {
    int r = 0;
    fflush(stdout);
    barrier();
    for (int p = 0; p < proc_number; p++) {
        parallel::barrier();
        if (proc_id != p)
            continue;
        //	fprintf(stdout, "%3d: ", proc_id);
        va_list ap;
        // stdout
        va_start(ap, fmt);
        r = vfprintf(stdout, fmt, ap);
        va_end(ap);
        fflush(stdout);
    }
    fflush(stdout);
    barrier();
    return (r);
}

int parallel::printf_master(int id, const char *fmt, ...) { // Write to stdout from Master process
    int r = 0;
    va_list ap;
    if (id == MASTER_ID) {
        va_start(ap, fmt);
        r = vfprintf(stdout, fmt, ap);
        va_end(ap);
    }
    fflush(stdout);
    return r;
}

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

    set<int>::iterator it_temp = temp.find(self_id);
    temp.erase(it_temp);

    for (set<int>::iterator it = temp.begin(); it != temp.end(); ++it) {
        list_of_neighbors[*it] = cur_idx;
        cur_idx++;
    }
}

void parallel::build_list_send_recv(VariableSizeMeshContainer<int> &topoNN, vector<int> &part,
                                    map<int, int> &list_of_neighbors, vector<set<int>> &send, vector<set<int>> &recv,
                                    int self_id) {

    int key, position;
    recv.clear();
    send.clear();

    recv.resize(list_of_neighbors.size());
    send.resize(list_of_neighbors.size());

    for (size_t i = 0; i < topoNN.getBlockNumber(); ++i) {
        for (size_t j = 0; j < topoNN.getBlockSize(i); ++j) {
            key = topoNN[i][j];

            if (part[i] == self_id && part[key] != self_id) {
                position = list_of_neighbors[part[key]]; // neighbourn number
                recv[position].insert(topoNN[i][j]);     // recv
                send[position].insert(i);                // send
            }
        }
    }
}

void parallel::gather_all(Decision *total, const vector<double> &local_solution, size_t n_own, const vector<int> &L2G) {
    int *recvcounts = nullptr;
    int *offsets = nullptr;
    size_t size = 0;

    Decision *buf_send;
    // Creating and registring custom data type
    MPI_Datatype builtin_types[2] = {MPI_INT, MPI_DOUBLE};
    int blocklens[2] = {1, 1};

    MPI_Aint displs[2];
    displs[0] = offsetof(Decision, id);
    displs[1] = offsetof(Decision, answer);

    MPI_Datatype Mpi_decision_datatype;
    MPI_Type_create_struct(2, blocklens, displs, builtin_types, &Mpi_decision_datatype);
    MPI_Type_commit(&Mpi_decision_datatype);
    ///////////////////////////////////// Registered

    int mpi_res;

    // Init
    if (proc_id == MASTER_ID) {
        recvcounts = new int[proc_number];
        offsets = new int[proc_number];
    }

    // Size
    size = n_own;

    // Gathering sizes
    mpi_res = MPI_Gather(&size, 1, MPI_INT, recvcounts, 1, MPI_INT, MASTER_ID, MCW);

    if (mpi_res != MPI_SUCCESS) {
        cout << "Sizes gathering error!" << endl;
        MPI_Finalize();

        if (proc_id == MASTER_ID) {
            delete[] recvcounts;
            delete[] offsets;
        }
        exit(1);
    }

    if (proc_id == MASTER_ID) {
        offsets[0] = 0;
        for (int i = 1; i < proc_number; ++i) {
            offsets[i] = offsets[i - 1] + recvcounts[i - 1];
        }
    }

    // Copying data...
    buf_send = new Decision[size];
    for (size_t i = 0; i < n_own; ++i) {
        buf_send[i].id = L2G[i];
        buf_send[i].answer = local_solution[i];
    }

    mpi_res = MPI_Gatherv(buf_send, size, Mpi_decision_datatype, total, recvcounts, offsets, Mpi_decision_datatype,
                          MASTER_ID, MCW);

    if (mpi_res != MPI_SUCCESS) {
        cout << "Primal gathering error!" << endl;
        MPI_Finalize();

        if (proc_id == MASTER_ID) {
            delete[] recvcounts;
            delete[] offsets;
            delete[] total;
        }

        delete[] buf_send;
        exit(1);
    }

    if (proc_id == MASTER_ID) {
        delete[] recvcounts;
        delete[] offsets;
    }
    delete[] buf_send;

    MPI_Type_free(&Mpi_decision_datatype);
}

void parallel::update_halo(vector<double> &x, map<int, int> &list_of_neighbors, vector<set<int>> &send,
                           vector<set<int>> &recv, int processor_id, MPI_Comm mpi_comm) {
    size_t scheme_size = list_of_neighbors.size();
    static vector<Ne_scheme_bufs> scheme(scheme_size);

    size_t j = 0, i = 0;
    size_t send_size, recv_size;
    size_t request_size = 2 * scheme_size;

    MPI_Request *requests = new MPI_Request[request_size];
    MPI_Status *statuses = new MPI_Status[request_size];

    int mpi_res;
    for (auto it = list_of_neighbors.cbegin(); it != list_of_neighbors.cend(); ++it, ++j) {

        send_size = send[j].size();
        recv_size = recv[j].size();

        scheme[j].neighbour_id = it->first; // for whom

        scheme[j].send_buf = new double[send[j].size()]; // interface for neighbour
        scheme[j].recv_buf = new double[recv[j].size()]; // halo for us

        i = 0;
        for (auto send_it = send[j].cbegin(); send_it != send[j].cend(); ++send_it, ++i) {
            scheme[j].send_buf[i] = x[*send_it];
        }

        mpi_res =
            MPI_Isend(scheme[j].send_buf, send_size, MPI_DOUBLE, scheme[j].neighbour_id, 0, mpi_comm, &requests[2 * j]);

        if (mpi_res != MPI_SUCCESS) {
            cerr << "MPI Isend error: processor " << processor_id << " crashed! Exiting..." << endl;
            exit(1);
        }

        mpi_res = MPI_Irecv(scheme[j].recv_buf, recv_size, MPI_DOUBLE, scheme[j].neighbour_id, 0, mpi_comm,
                            &requests[2 * j + 1]);

        if (mpi_res != MPI_SUCCESS) {
            cerr << "MPI Irevc error: processor " << processor_id << " crashed! Exiting..." << endl;
            exit(1);
        }
    }

    mpi_res = MPI_Waitall(request_size, requests, statuses); // blocking
    if (mpi_res != MPI_SUCCESS) {
        printf_master(processor_id, "MPI Waitall error: processor crashed! Exiting... %d", processor_id);
        MPI_Finalize();
        exit(1);
    }

    // considering halo nodes in x are sorted by owners order
    i = 0;
    for (auto it = list_of_neighbors.cbegin(); it != list_of_neighbors.cend(); ++it, ++i) {
        j = 0;
        for (auto recv_it = recv[i].cbegin(); recv_it != recv[i].cend(); ++recv_it, ++j) {
            x[*recv_it] = scheme[i].recv_buf[j];
        }
    }
}
