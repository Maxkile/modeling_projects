// clang-format off
#include "IO.hpp"

void IO::draw_mesh(int Nx, int Ny, int k3, int k4){
    int cur_number_type = k3 > 0 ? k3:k4;
    int type_elem = k3 > 0 ? 0:1;
    int prev_cur_number = cur_number_type;
    int prev_type = type_elem;
    char symb = ' ';

    cout << "\n";

    for (int i = 0; i < Nx; i++) {
        cout << " ";
        for (int j = 0; j < 8; j++)
            cout << "_";
    }
    cout << "\n";

    for (int cur_line = 0; cur_line < Ny; cur_line++) {
        symb = ' ';

        for (int cur_str = 0; cur_str < 4; cur_str++) {
            cout << "|";
            if (cur_str == 3)
                symb = '_';

            type_elem = prev_type;
            cur_number_type = type_elem == 0 ? k3:k4;
            if ((cur_line != 0) && (prev_cur_number != 0 )){
                cur_number_type = prev_cur_number;
            }

            for (int i = 0; i < Nx; i++) {
                for (int j = 0; j < 7 - 2 * cur_str; j++)
                    cout << symb;
                if (type_elem == 0) {
                    cout << "̸"; //U+0338
                }
                else {
                    cout << symb;
                }
                for (int j = 0; j < 2 * cur_str; j++)
                    cout << symb;
                cout << "|";

                cur_number_type--;
                if (cur_number_type == 0){
                    type_elem = (type_elem + 1) % 2;
                    cur_number_type = type_elem == 0 ? k3:k4;
                    if (cur_number_type == 0) {
                        type_elem = (type_elem + 1) % 2;
                        cur_number_type = type_elem == 0 ? k3:k4;
                    }
                }
            }
            cout << "\n";
        }
        prev_cur_number = cur_number_type;
        prev_type = type_elem;
    }

    cout << "\n";
}

int IO::num_elem(int Nx, int Ny, int k3, int k4){
    int nE;

    nE = (2 * k3 + k4) * int(((Nx - 1) * (Ny - 1)) / (k3 + k4));
    if ((((Nx - 1) * (Ny - 1)) % (k3 + k4)) >= k3)
        nE += 2 * k3 + (((Nx - 1) * (Ny - 1)) % (k3 + k4)) - k3;
    else
        nE += 2 * (((Nx - 1) * (Ny - 1)) % (k3 + k4));

    return nE;
}

void IO::MPI_gather_write(const std::string &filename, const vector<double> &solution, const size_t n_own,
                      const vector<int> &L2G, const size_t totalSize){
    std::ofstream fout;

    if (proc_id == MASTER_ID) {
        printf(LINE_SEPARATOR);
        cout << "\t\t Gather writing" << endl;
        printf(LINE_SEPARATOR);
        fout.flush();
    }

    parallel::Decision* total = nullptr;
    if (proc_id == MASTER_ID) {
        total = new parallel::Decision[totalSize];
    }

    parallel::gather_all(total, solution, n_own, L2G);

    if (proc_id == MASTER_ID) {
        fout.open(filename);
        fout << "Id: "
             << "    |   "
             << "Value:" << std::endl;
        for (size_t i = 0; i < totalSize; ++i) {
            fout << total[i].id << " " << total[i].answer << endl;
        }
        fout.flush();
        fout.close();
    }
    parallel::barrier();
    if (proc_id == MASTER_ID){
        cout << "Gather writing done" << endl;
        delete[] total;
    }
}

void IO::MPI_self_write(const std::string& filename,const vector<double>& solution,const size_t n_own,const vector<int>& L2G){
    std::ofstream fout;

    if (proc_id == MASTER_ID) {
        printf(LINE_SEPARATOR);
        cout << "\t\t Self process writing" << endl;
        printf(LINE_SEPARATOR);
        fout.flush();
    }
    parallel::barrier();
    for (int i = 0; i < proc_number; ++i) {
        if (i == proc_id) {

            std::cout << "Process: " << proc_id << " writing solution to file..." << endl;
            if (proc_id == 0) {
                fout.open(filename);
                fout << "Id: "
                     << "    |   "
                     << "Value:" << std::endl;
            } else {
                fout.open(filename, ios_base::app);
            }
            for (size_t i = 0; i < n_own; i++) {
                fout << L2G[i] << " " << solution[i] << endl;
            }
            cout << "Done" << endl;
            fout.flush();
            fout.close();
        }
        parallel::barrier();
    }
}

// clang-format on
