#pragma once
#include "stdafx.hpp"

#include <iomanip>
#include <cmath>
#include <omp.h>

#include "platformDependencies.hpp"
#include "Sparse.hpp"


//vector matrix operations: consider we have same template arguments
namespace vmo
{
	template<typename VT, typename ST>
	std::vector<VT> axpby(const std::vector<VT>& x, const std::vector<VT>& y, const ST& a, const ST& b, unsigned short threadsNumber = 4)
	{
		std::vector<VT> result;
		size_t size = x.size();
		if (size != y.size())
		{
			std::cerr << "Incompatible sizes in axpby!" << std::endl;
			return result;
		}
		else
		{
			result.resize(size);

			#pragma omp parallel for num_threads(threadsNumber) 
			for (OPENMP_INDEX_TYPE i = 0; i < size; ++i)
			{
				result[i] = a * x[i] + b * y[i];
			}
			return result;
		}
	}


	template<typename VT>
	VT dot(const std::vector<VT>& x, const std::vector<VT>& y, unsigned short threadsNumber = 4)
	{
		VT result = 0;
		size_t size = x.size();

		if (size != y.size())
		{
			std::cerr << "Incompatible sizes in dot!" << std::endl;
			return result;//incorrent
		}
		else
		{
			#pragma omp parallel for num_threads(threadsNumber) reduction(+:result)
			for (OPENMP_INDEX_TYPE i = 0; i < size; ++i)
			{
				result += x[i] * y[i];
			}
			return result;
		}
	}

	template<typename VT>
	std::vector<VT> multiply(const std::vector<VT>& x, VT scalar, unsigned short threadsNumber = 4)
	{
		std::vector<VT> result;

		result.resize(x.size());

		#pragma omp parallel for num_threads(threadsNumber) 
		for (OPENMP_INDEX_TYPE i = 0; i < x.size(); ++i)
		{
			result[i] = scalar * x[i];
		}
		return result;
	}
	

	// Compute norm
	double norm(const std::vector<double>& vec){
		return sqrt(dot(vec,vec));
	}

#if 0
	// Compute determinant matrix T
	template <typename M>
	M det(M** T, int N)
	{
		M det__;
		int sub_j, s;
		M** subT;    // Субматрица как набор ссылок на исходную матрицу

		switch (N)
		{
		case 1:
			return T[0][0];
		case 2:
			return T[0][0] * T[1][1] - T[0][1] * T[1][0];
		default:
			if (N < 1)
			{
				throw 0;  // Некорректная матрица
			}
			subT = new M*[N - 1];  // Массив ссылок на столбцы субматрицы
			det__ = 0;
			s = 1;        // Знак минора
			for (int i = 0; i < N; i++)  // Разложение по первому столбцу
			{
				sub_j = 0;
				for (int j = 0; j < N; j++)// Заполнение субматрицы ссылками на исходные столбцы
				{
					if (i != j)      // исключить i строку
					{
						subT[sub_j++] = T[j] + 1;  // здесь + 1 исключает первый столбец
					}

				}
				det__ = det__ + s * T[i][0] * det(subT, N - 1);
				s = -s;
			};
			delete[] subT;
			return det__;
		};
	};


	// Check symmetric
	template <typename M>
	bool isSymmetric(M** mas, int n)
	{
		for (int i = 0; i < n - 1; ++i)
		{
			for (int j = i + 1; j < n; ++j)
			{
				if (mas[i][j] != mas[j][i])
				{
					return false;
				}
			}
		}
		return true;
	}

	// Check A > 0
	template <typename M>
	bool isPositive(M** mas, int n)
	{

		M** minor = nullptr;
		M value_minor;

		if (!isSymmetric(mas, n))
		{
			cout << "Matrix is not symmetric!" << endl;
			return false;
		}

		for (int i = 1; i <= n; i++)
		{
			if (i == 1)
			{
				value_minor = mas[0][0];
			}
			else
			{
				minor = new M*[i];
				for (int j = 1; j <= i; j++)
				{
					minor[j - 1] = new M[i];
				}
				for (int k = 0; k < i; k++)
				{
					for (int l = 0; l < i; l++)
					{
						minor[k][l] = mas[k][l];
					}
					value_minor = det(minor, i);
				}
				for (int j = 1; j <= i; j++)
				{
					delete[] minor[j - 1];
				}
				delete[] minor;
			}
			if (value_minor <= 0)
			{
				return false;
			}
		}
		return true;
	}
#endif


	//A = A^(T) > 0
	template<typename M, typename V>
	std::vector<double> conGradSolver(Sparse<M>& A, const std::vector<V>& b, size_t threadsNum = 4, size_t n_max = 100, double eps = 0.001, bool mode = 1)
	{
		std::vector<double> x_prev(A.getDenseColumns());
		if (A.getDenseColumns() != b.size())
		{
			std::cerr << "Incompatible sizes in Solver!" << std::endl;
		}
		else
		{
			

			size_t size = A.getDenseColumns();

			double** rev_M = new double*[size];
			for (size_t i = 0; i < size; ++i)
			{
				rev_M[i] = new double[size];
			}

			std::vector<M> diagonal = A.getDiagonal();

			//forming preconditioning matrix - M
			for (size_t i = 0; i < size; ++i)
			{
				for (size_t j = 0; j < size; ++j)
				{
					rev_M[i][j] = 0.0;
				}
				rev_M[i][i] = 1 / (diagonal[i] + 0.001);
				std::cout << diagonal[i] << std::endl;
			}

			SparseELL<double> reverse_M(rev_M, size, size);//preconditioner

			for (size_t i = 0; i < size; ++i)
			{
				delete[] rev_M[i];	
			}
			delete[] rev_M;

			//initializang parameters for algorithm

			std::vector<double> x_cur;

			std::vector<double> r_cur;
			std::vector<double> r_prev = b;//we think that x0 = (0,0,0,0,0,0 ... 0)^(T) -> b - Ax = b

			std::vector<double> p_cur;//p and p+1 are conjugated relative to A
			std::vector<double> p_prev;

			std::vector<double> q;
			std::vector<double> z;

			double delta_cur, delta_prev;
			double alpha, beta;

			size_t k = 1;

			do
			{
				z = reverse_M.spmv(r_prev, threadsNum);
				delta_cur = dot(r_prev, z, threadsNum);
				if (k == 1)
				{
					p_cur = z;
				}
				else
				{
					beta = delta_cur / delta_prev;
					p_cur = axpby(z, multiply(p_prev, beta, threadsNum), 1, 1, threadsNum);
				}
				q = A.spmv(p_cur, threadsNum);
				alpha = delta_cur / dot(p_cur, q);

				x_cur = axpby(x_prev, multiply(p_cur, alpha, threadsNum), 1, 1, threadsNum);

				r_cur = axpby(r_prev, multiply(q, alpha, threadsNum), 1, -1, threadsNum);
				
				if (mode)
				{
					if (k == 1)
					{
						std::cout << "| Iteration | norm value |" << endl;
					}
					std::cout << "|" << setw(11) << k << "|" << setw(12) << fixed << setprecision(3) << norm(r_cur) << "|" << endl;
				}
				
				if ((norm(p_cur) < eps) || (k > n_max))
				{	
					break;
				}

				r_prev = r_cur;
				x_prev = x_cur;
				p_prev = p_cur;
				delta_prev = delta_cur;
				
				k += 1;
				
			} while (true);

		}

		return x_prev;
	}
}

