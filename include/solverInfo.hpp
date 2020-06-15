#ifndef _SOLVER_INFO_HPP
#define _SOLVER_INFO_HPP

#pragma once

enum class SparseType : unsigned int { CSR, ELLPACK, COO };

enum class OutputStrategy : unsigned int { SEPARATE, GATHER };

enum class OutputType : unsigned int { STDOUT, FILE };

struct SolverInfo {
    unsigned threadsNumber;
    OutputStrategy outputStrategy;
    OutputType outputType;
    SparseType sparseMatrixType;
};

#endif
