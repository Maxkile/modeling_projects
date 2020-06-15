// clang-format off
#include "IO.hpp"

int decision_comparator(const void* dec1,const void* dec2){
    const parallel::Decision *decision_1 = static_cast<const parallel::Decision*>(dec1);
    const parallel::Decision *decision_2 = static_cast<const parallel::Decision*>(dec2);

    int delta = decision_1->id > decision_2->id;
    if (delta < 0){
        return -1;
    } else if (delta > 0){
        return 1;
    }
    else{
        return 0;
    }

}

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

void IO::mpi_gather_write(OutputType solverOutType, const vector<double>& solution, const size_t n_own,
    const vector<int>& L2G, const size_t totalSize) {

    double time;

    parallel::printf_master(proc_id, LINE_SEPARATOR, "\t\t Gather writing",LINE_SEPARATOR);

    parallel::Decision* total = nullptr;
    if (proc_id == MASTER_ID) {
        total = new parallel::Decision[totalSize];
    }


    parallel::gather_all(total, solution, n_own, L2G);

    if (proc_id == MASTER_ID) {
        time = omp_get_wtime();
        qsort(total,totalSize,sizeof(parallel::Decision),decision_comparator);//
        time = omp_get_wtime() - time;
    }

    if (solverOutType == OutputType::STDOUT) {
        if (proc_id == MASTER_ID) {
            cout << "Id: "
                << "    |   "
                << "Value:" << std::endl;
            for (size_t i = 0; i < totalSize; ++i) {
                cout << total[i].id << " " << total[i].answer << endl;
            }
            cout.flush();
        }
    }
    else {
        std::ofstream fout;
        if (proc_id == MASTER_ID) {
            fout.open("solution_gather.txt");
            fout << "Id: "
                << "    |   "
                << "Value:" << std::endl;
            for (size_t i = 0; i < totalSize; ++i) {
                fout << total[i].id << " " << total[i].answer << endl;
            }
            fout.flush();
            fout.close();
        }
    }
        
    parallel::barrier();
    parallel::printf_master(proc_id, "Gather writing done.\nQuick sorting time: %.6g\n", time);
    if (proc_id == MASTER_ID){
        delete[] total;
    }
}

void IO::mpi_self_write(OutputType solverOutType, const vector<double>& solution, const size_t n_own, const vector<int>& L2G) {

    parallel::printf_master(proc_id, LINE_SEPARATOR, "\t\t Gather writing",LINE_SEPARATOR);
    parallel::barrier();

    if (solverOutType == OutputType::STDOUT) {
        for (int i = 0; i < proc_number; ++i) {
            if (i == proc_id) {
                printf(LINE_SEPARATOR);
                cout << "\t\t Process " << proc_id  << endl;
                printf(LINE_SEPARATOR);
                cout.flush();
                cout << "Id: "
                    << "    |   "
                    << "Value:" << std::endl;
                for (size_t i = 0; i < n_own; i++) {
                    cout << L2G[i] << " " << solution[i] << endl;
                }
                cout << "Done" << endl;
                cout.flush();
            }
            parallel::barrier();
        }
    }
    else {
        std::ofstream fout;
        for (int i = 0; i < proc_number; ++i) {
            if (i == proc_id) {

                std::cout << "Process: " << proc_id << " writing solution to file..." << endl;
                if (proc_id == 0) {
                    fout.open("solution.txt");
                    fout << "Id: "
                        << "    |   "
                        << "Value:" << std::endl;
                }
                else {
                    fout.open("solution.txt", ios_base::app);
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
}

// clang-format on
