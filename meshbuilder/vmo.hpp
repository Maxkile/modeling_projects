#pragma once
#include "stdafx.hpp"

#include <iomanip>
#include <cmath>
#include <omp.h>

#include "platformDependencies.hpp"
#include "Sparse.hpp"

extern double timeSpmv;

//vector matrix operations: consider we are having same template arguments
namespace vmo
{

	static double timeDot = 0;
	static double timeMultiply = 0;
	static double timeAxpby = 0;

	template<typename VT, typename ST>
	void axpby(std::vector<VT>& x, const std::vector<VT>& y, const ST& a, const ST& b, size_t threadsNumber = 4)
	{
		size_t size = x.size();
		if (size != y.size())
		{
			std::cerr << "Incompatible sizes in axpby!" << std::endl;
		}
		else
		{
			double start = omp_get_wtime();
			#pragma omp parallel for num_threads(threadsNumber) 
			for (OPENMP_INDEX_TYPE i = 0; i < size; ++i)
			{
				x[i] = a * x[i] + b * y[i];
			}
			timeAxpby+= omp_get_wtime() - start;
		}
	}


	template<typename VT>
	VT dot(const std::vector<VT>& x, const std::vector<VT>& y, size_t threadsNumber = 4)
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
			double start = omp_get_wtime();
			#pragma omp parallel for num_threads(threadsNumber) reduction(+:result)
			for (OPENMP_INDEX_TYPE i = 0; i < size; ++i)
			{
				result += x[i] * y[i];
			}
			timeDot+= omp_get_wtime() - start;
			return result;
		}
	}

	template<typename VT>
	void multiply(const std::vector<VT>& x, std::vector<VT>& result, VT scalar, size_t threadsNumber = 4)
	{
		double start = omp_get_wtime();
		#pragma omp parallel for num_threads(threadsNumber) 
		for (OPENMP_INDEX_TYPE i = 0; i < x.size(); ++i)
		{
			result[i] = scalar * x[i];
		}
		timeMultiply+= omp_get_wtime() - start;
	}
	

	// Compute eucledian norm
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

	// Compute eucledian norm
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
	std::vector<double> conGradSolver(Sparse<M>& A, const std::vector<V>& b, size_t threadsNum = 4, size_t n_max = 100, double eps = 0.0000001)
	{
		std::vector<double> x_prev(A.getDenseColumns());
		if (A.getDenseColumns() != b.size())
		{
			std::cerr << "Incompatible sizes in Solver!" << std::endl;
		}
		else
		{
			size_t size = A.getDenseColumns();
		
			std::vector<M> diagonal = A.getDiagonal();
			std::vector<M> rev_M(diagonal.size());
			
			//forming preconditioning matrix - M
			for (size_t i = 0; i < size; ++i)
			{
				rev_M[i] = 1 / (diagonal[i] + eps);
			}
				
			SparseELL<double> reverse_M(rev_M);//preconditioner from diagonal vector

			//initializang parameters for algorithm

			std::vector<double> x_cur;

			std::vector<double> r_prev = b;//we think that x0 = (0,0,0,0,0,0 ... 0)^(T) -> b - Ax = b

			std::vector<double> p_cur;//p and p+1 are conjugated relative to A
			std::vector<double> p_prev;

			std::vector<double> q(size);
			std::vector<double> z(size);
			std::vector<double> multiplyResult(size);

			double delta_cur, delta_prev;
			double alpha, beta;
			double b_norm = norm(b);

			size_t k = 1;

			std::cout << "| Iteration | Norm value |" << endl;

			double conGradSolverTime = omp_get_wtime();

 			do
			{
				reverse_M.spmv(r_prev, z, threadsNum);
				delta_cur = dot(r_prev, z, threadsNum);

				if (k == 1)
				{
					p_cur = z;
				}
				else
				{
					beta = delta_cur / delta_prev;
					multiply(p_prev, multiplyResult, beta ,threadsNum);
					axpby(multiplyResult, z, 1, 1, threadsNum);
					p_cur = multiplyResult;//don't know how to improve
				}

				A.spmv(p_cur, q, threadsNum);
				alpha = delta_cur / dot(p_cur, q);

				multiply(p_cur, multiplyResult, alpha, threadsNum);
				axpby(x_prev, multiplyResult, 1, 1, threadsNum);

				multiply(q, multiplyResult, alpha, threadsNum);
				axpby(r_prev, multiplyResult, 1, -1, threadsNum);

				std::cout << "|" << setw(11) << k << "|" << setw(12) << fixed << setprecision(7) << norm(r_prev)/b_norm << "|" << endl;
				
				p_prev = p_cur;
				delta_prev = delta_cur;
				
				k += 1;
				
			} while ((norm(r_prev)/b_norm >= eps) && (k <= n_max));
		
			conGradSolverTime = omp_get_wtime() - conGradSolverTime;

			std::cout << std::endl;
			std::cout << "Dot functions total time: " << timeDot << std::endl;
			std::cout << "Axpby functions total time: " << timeAxpby << std::endl;
			std::cout << "Multiply functions total time: " << timeMultiply << std::endl;
			std::cout << "Spmv functions total time: " << timeSpmv << std::endl;
			std::cout << std::endl;
			std::cout << "Conjugate gradient method solution total time: " << conGradSolverTime << std::endl;
			std::cout << std::endl;
		}
		
		return x_prev;
	}
}

