Simple MPI and OpenMP solver. Supports CSR, COO and ELLPACK sparse matrix formats. Vector "b" is set manually.

Example: "mpiexec -n 4 meshbuilder --gen 6 6 100 1000 3 3 2 2 --solver -coo 1"

Result is written to the "solution.txt" and "solution_gather.txt".
Vector-column b by default is [1,2,3,4,5,6, ... ,n]

To build:
  1) cmake CmakeLists.txt
  2) make -f Makefile

Note: make sure that composition of submeshes(Px,Py) must be equal to value "-n"(processes number) or decomposition will have no effect.


