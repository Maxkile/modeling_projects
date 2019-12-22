For Keldysh Institute.

Builds topologies in the beginning, sets sin(i + j) values(set_model_values method) and solves.

Example command: ./prod --gen 6 6 10 10 3 2 --solver csr 4
Result Ax=b write into the file "decision.txt"
Vector-column b by default is [1,2,3,4,5,6, ... ,n]


