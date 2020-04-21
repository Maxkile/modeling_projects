#include "IO.hpp"

void draw_mesh(int Nx, int Ny, int k3, int k4){
    int cur_number_type = k3 > 0 ? k3:k4;
    int type_elem = k3 > 0 ? 0:1;
    int prev_cur_number = cur_number_type;
    int prev_type = type_elem;
    char symb = ' ';

    cout << "\n";

    for (int i = 0; i < Ny; i++) {
        cout << " ";
        for (int j = 0; j < 8; j++)
            cout << "_";
    }
    cout << "\n";

    for (int cur_line = 0; cur_line < Nx; cur_line++) {
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

            for (int i = 0; i < Ny; i++) {
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

int num_elem(int Nx, int Ny, int k3, int k4){
    int nE;

    nE = (2 * k3 + k4) * int(((Nx - 1) * (Ny - 1)) / (k3 + k4));
    if ((((Nx - 1) * (Ny - 1)) % (k3 + k4)) >= k3)
        nE += 2 * k3 + (((Nx - 1) * (Ny - 1)) % (k3 + k4)) - k3;
    else
        nE += 2 * (((Nx - 1) * (Ny - 1)) % (k3 + k4));

    return nE;
}
