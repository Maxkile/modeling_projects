#pragma once
#include "stdafx.hpp"

#include <omp.h>
#include <algorithm>
#include <iomanip>
#include <cmath>
#include "platformDependencies.hpp"
#include "VariableSizeMeshContainer.hpp"

double timeSpmv = 0;

template<typename T>
class Sparse
{
protected:

	size_t denseRows, denseColumns;
	std::vector<size_t> JA;
	std::vector<T> A;

public:

	virtual void spmv(const std::vector<T>&, std::vector<T>& result, size_t threadsNumber = 4) const = 0;
	virtual void spmv(const T*, std::vector<T>& result, size_t, size_t threadsNumber = 4) const = 0;
	virtual T** getDenseMatrix() const = 0;
	virtual void setValues(const std::vector<T>&) = 0;
	virtual size_t getValuesSize() const = 0;
	virtual std::vector<T> getDiagonal() const = 0;
	
	
	inline void printA() const
	{
		std::cout << "A vector: ";
		for (size_t i = 0; i < A.size(); ++i)
		{
			std::cout << A[i] << " ";
		}
		std::cout << std::endl;
	}

	inline void printJa() const
	{
		std::cout << "JA vector: ";
		for (size_t i = 0; i < JA.size(); ++i)
		{
			std::cout << JA[i] << " ";
		}
		std::cout << std::endl;
	}

	size_t getDenseRows() const
	{
		return denseRows;
	}

	size_t getDenseColumns() const
	{
		return denseColumns;
	}

	void outputMatrix() const
	{
		T** matr = this->getDenseMatrix();

		std::cout << "Matrix: " << endl;
		for (size_t i = 0; i < denseRows; i++)
		{
			for (size_t j = 0; j < denseRows; j++)
			{
				std::cout << fixed << setprecision(3) << matr[i][j] << " ";
			}
			std::cout << endl;
		}
		std::cout << endl;

		for (size_t i = 0; i < this->denseRows; ++i)
		{
			delete[] matr[i];
		}
		delete[] matr;
	}

	virtual ~Sparse(){}
	
};

template<typename T>
class SparseELL : public Sparse<T>
{
	size_t rowOffset, valuesSize;


	size_t findRowOffset(const VariableSizeMeshContainer<int>& topoNN)
	{
		size_t maxRowOffset = 1;

		for (size_t i = 0; i < topoNN.getBlockNumber(); ++i)
		{
			if (topoNN.getBlockSize(i) > maxRowOffset)
			{
				maxRowOffset = topoNN.getBlockSize(i);
			}
		}
		return (maxRowOffset + 1);//"+1" - for diagonal
	}

	size_t findRowOffset(T** dense,size_t rows,size_t columns)
	{
		size_t maxRowOffset = 1;
		size_t rowElementCount;

		for (size_t i = 0; i < rows; ++i)
		{
			rowElementCount = 1;//"1" - for diagonal
			for (size_t j = 0; j < columns; ++j)
			{
				if (dense[i][j] != 0)
				{
					rowElementCount++;
				}
			}
			if (maxRowOffset < rowElementCount)
			{
				maxRowOffset = rowElementCount;
			}
		}
		return maxRowOffset;
	}

	void setPortraitValues()
	{
		size_t index_i, index_j, diagIndex;
		diagIndex = 0;
		T rowSum;
		for (size_t i = 0; i < this->denseRows; ++i)
		{
			rowSum = 0;
			for (size_t j = i * this->rowOffset; j < (i + 1) * this->rowOffset; ++j)
			{
				if (this->A[j] != 0)//padding
				{
					index_i = i;
					index_j = this->JA[j];
					if (index_i != index_j)
					{
						this->A[j] = sin(index_i + index_j);
						rowSum += abs(this->A[j]);
					}
					else
					{
						diagIndex = j;
					}
				}
			}
			this->A[diagIndex] = rowSum * 1.5;
		}
	}

public:

	SparseELL(const std::vector<T>& diagonal)//vector is put on primary diagonal,no padding
	{
		this->denseRows = this->denseColumns = this->valuesSize = diagonal.size();
		this->rowOffset = 1;//diagonal

		this->JA.reserve(this->rowOffset * diagonal.size());
		this->A.reserve(this->rowOffset * diagonal.size());

		for (size_t i = 0; i < this->denseRows; ++i)
		{
			this->JA.push_back(i);
			this->A.push_back(diagonal[i]);
		}
	}

	SparseELL(T** dense, size_t rows, size_t columns)//from dense(needed in solver)
	{
		this->valuesSize = 0;
		this->denseRows = rows;
		this->denseColumns = columns;
		this->rowOffset = findRowOffset(dense,rows,columns);

		this->JA.reserve(this->rowOffset * rows);
		this->A.reserve(this->rowOffset * rows);
		
		size_t rowElemCount;

		for (size_t i = 0; i < rows; ++i)
		{
			rowElemCount = 0;
			
			for (size_t j = 0; j < columns; ++j)
			{
				if (dense[i][j] != 0)
				{
					this->JA.push_back(j);
					this->A.push_back(dense[i][j]);
					this->valuesSize++;
					rowElemCount++;
				}
			}
			
			size_t paddingIndex = this->JA[this->JA.size() - 1];

			for (size_t j = rowElemCount; j < rowOffset; ++j)
			{
				this->JA.push_back(paddingIndex);
				this->A.push_back(0);//padding
			}
		}
	}

	SparseELL(const VariableSizeMeshContainer<int>& topoNN)//from topo
	{
		this->valuesSize = 0;
		this->denseRows = this->denseColumns = topoNN.getBlockNumber();
		this->rowOffset = findRowOffset(topoNN);

		this->JA.reserve(this->rowOffset * topoNN.getBlockNumber());
		this->A.reserve(this->rowOffset * topoNN.getBlockNumber());

		std::vector<int> temp;

		for (size_t i = 0; i < topoNN.getBlockNumber(); ++i)
		{
			temp.push_back(i);//diag
			for (size_t j = 0; j < topoNN.getBlockSize(i); ++j)
			{
				temp.push_back(topoNN[i][j]);
			}

			std::sort(temp.begin(), temp.end());

			for (size_t j = 0; j < temp.size(); ++j)
			{
				this->JA.push_back(temp[j]);
				this->A.push_back(1);
				this->valuesSize++;
			}

			size_t paddingIndex = this->JA[this->JA.size() - 1];

			for(size_t j = temp.size(); j < rowOffset;++j)
			{
				this->JA.push_back(paddingIndex);
				this->A.push_back(0);//padding
			}
			temp.clear();
		}
		setPortraitValues();
	}

	T** getDenseMatrix() const override
	{
		T** dense = new T*[this->denseRows];
		for (size_t i = 0; i < this->denseRows; ++i)
		{
			dense[i] = new T[this->denseColumns];
			for (size_t j = 0; j < this->denseColumns; ++j)
			{
				dense[i][j] = 0;
			}
		}


		for (size_t i = 0; i < this->denseRows; ++i)
		{
			for (size_t j = i * this->rowOffset; j < (i + 1) * this->rowOffset; ++j)
			{
				if (this->A[j] != 0)//padding
				{
					dense[i][this->JA[j]] = this->A[j];
				}
			}
		}
		return dense;
	}

	void spmv(const std::vector<T>& x, std::vector<T>& result, size_t threadsNumber = 4) const override
	{
		if (x.size() != this->denseColumns)
		{
			std::cerr << "Incompatible sizes in spmv!" << std::endl;
		}
		else
		{
			double start = omp_get_wtime();
			#pragma omp parallel for num_threads(threadsNumber)//TODO: ask about template reduction
			for (OPENMP_INDEX_TYPE i = 0; i < this->denseRows; ++i)
			{
				result[i] = 0;

				for (OPENMP_INDEX_TYPE j = i * this->rowOffset; j < (i + 1) * this->rowOffset; ++j)
				{
					result[i] += x[this->JA[j]] * this->A[j];
				}
			}
			timeSpmv += omp_get_wtime() - start;
		}
	}

	void spmv(const T* x, std::vector<T>& result, size_t xSize, size_t threadsNumber = 4) const override
	{

		if (xSize != this->denseColumns)
		{
			std::cerr << "Incompatible sizes in spmv!" << std::endl;
		}
		else
		{
			double start = omp_get_wtime();
			#pragma omp parallel for num_threads(threadsNumber) 
			for (OPENMP_INDEX_TYPE i = 0; i < this->denseRows; ++i)
			{
				result[i] = 0;

				for (OPENMP_INDEX_TYPE j = i * rowOffset; j < (i + 1) * rowOffset; ++j)
				{
					result[i] += x[this->JA[j]] * this->A[j];
				}
			}
			timeSpmv+= omp_get_wtime() - start;
		}
	}

	void setValues(const std::vector<T>& a) override
	{
		if (a.size() != this->getValuesSize())
		{
			std::cerr << "Incompatibe sizes of values vectors!" << std::endl;
		}
		else
		{
			for (size_t i = 0, j = 0; j < a.size(); i++)//i - for A, j - for vector "a"
			{
				if (this->A[i] != 0)
				{
					this->A[i] = a[j];
					j++;
				}
			}
		}
	}

	size_t getRowOffset() const
	{
		return rowOffset;
	}

	//must remember that have padding
	size_t getValuesSize() const override
	{
		return this->valuesSize;
	}

	std::vector<T> getDiagonal() const override
	{
		std::vector<T> diagonal;
		diagonal.reserve(this->denseColumns);
		for(size_t i = 1; i < this->denseRows + 1;++i)
		{			
			for(size_t j = (i - 1) * rowOffset; j < i * rowOffset;++j)
			{
				if ((this->JA[j] == i - 1) && (this->A[j] != 0))//avoiding padding : when padding index is index of diagonal
				{
					diagonal.push_back(this->A[j]);
				}
			}
		}
		return diagonal;		
	}
};

template <typename T>
class SparseCSR :public Sparse<T> {

	std::vector<size_t> IA;//sizes of A and JA are similar

	void setPortraitValues()
	{
		size_t index_i, index_j, diagIndex;
		diagIndex = 0;
		T rowSum;
		for (size_t i = 0; i < this->denseRows; ++i)
		{
			rowSum = 0;
			for (size_t j = IA[i]; j < IA[i + 1]; ++j)
			{
				index_i = i;
				index_j = this->JA[j];
				if (index_i != index_j)
				{
					this->A[j] = sin(index_i + index_j);
					rowSum += abs(this->A[j]);
				}
				else
				{
					diagIndex = j;
				}
			}
			this->A[diagIndex] = rowSum * 1.5;
		}
	}

public:

	SparseCSR(const VariableSizeMeshContainer<int>& topoNN)
	{
		std::vector<int> temp;
		this->denseColumns = this->denseRows = topoNN.getBlockNumber();
		int prev;

		this->JA.reserve(topoNN.getBlockNumber());
		this->A.reserve(topoNN.getBlockNumber());

		IA.push_back(0);
		for (size_t i = 0; i < topoNN.getBlockNumber(); i++)
		{
			temp.push_back(i);
			for (size_t j = 0; j < topoNN.getBlockSize(i); j++)
			{
				temp.push_back(topoNN[i][j]);
			}
			std::sort(temp.begin(), temp.end());
			for (size_t j = 0; j < temp.size(); j++)
			{
				this->JA.push_back(temp[j]);
			}
			prev = IA[IA.size() - 1];
			IA.push_back(temp.size() + prev);
			temp.clear();
		}
		IA.push_back(this->denseColumns);

		this->A.reserve(this->JA.size());
		for (size_t i = 0; i < this->JA.size(); ++i)
		{
			this->A.push_back(1);//portrait
		}
		setPortraitValues();
	}

	T** getDenseMatrix() const override
	{
		T** dense = new T*[this->denseRows];
		for (size_t i = 0; i < this->denseRows; ++i)
		{
			dense[i] = new T[this->denseColumns];
			for (size_t j = 0; j < this->denseColumns; ++j)
			{
				dense[i][j] = 0;
			}
		}

		for (size_t i = 0; i < this->denseRows; ++i)
		{
			for (size_t j = IA[i]; j < IA[i + 1]; ++j)
			{
				dense[i][this->JA[j]] = this->A[j];
			}
		}
		return dense;
	}

	void spmv(const std::vector<T>& x, std::vector<T>& result, size_t threadsNumber = 4) const override
	{
		if (x.size() != this->denseColumns)
		{
			std::cerr << "Incompatible sizes in spmv!" << std::endl;
		}
		else
		{
			double start = omp_get_wtime();
			#pragma omp parallel for num_threads(threadsNumber)
			for (OPENMP_INDEX_TYPE i = 0; i < this->denseRows; i++)
			{

				result[i] = 0;

				for (size_t j = IA[i]; j < IA[i + 1]; j++)
					result[i] += x[this->JA[j]] * (this->A[j]);
			}
			timeSpmv+= omp_get_wtime() - start;
		}
	}

	void spmv(const T* x, std::vector<T>& result, size_t xSize,  size_t threadsNumber = 4) const override
	{
		if (xSize != this->denseColumns)
		{
			std::cerr << "Incompatible sizes in spmv!" << std::endl;
		}
		else
		{
			double start = omp_get_wtime();
			#pragma omp parallel for num_threads(threadsNumber)
			for (OPENMP_INDEX_TYPE i = 0; i < this->denseRows; ++i) {

				result[i] = 0;

				for (size_t j = IA[i]; j < IA[i + 1]; j++)
					result[i] += x[this->JA[j]] * (this->A[j]);
			}
			timeSpmv+= omp_get_wtime() - start;
		}
	}

	void printIa() const
	{
		std::cout << "IA vector: ";
		for (size_t i = 0; i < IA.size(); ++i)
		{
			std::cout << IA[i] << " ";
		}
		std::cout << std::endl;
	}

	void setValues(const std::vector<T>& a) override
	{
		if (a.size() != this->getValuesSize())
		{
			std::cerr << "Incompatibe sizes of values vectors!" << std::endl;
		}
		else
		{
			for (size_t i = 0; i < a.size(); ++i)
			{
				this->A[i] = a[i];
			}
		}
	}

	size_t getValuesSize() const override
	{
		return this->A.size();
	}

	std::vector<T> getDiagonal() const override
	{
		std::vector<T> diagonal;
		diagonal.reserve(this->denseColumns);

		for(size_t i = 1; i < this->IA.size();++i)
		{			
			for(size_t j = this->IA[i - 1]; j < this->IA[i];++j)
			{
				if (this->JA[j] == i - 1)//diag elem
				{
					diagonal.push_back(this->A[j]);
				}
			}
		}
		return diagonal;
	}
	
};

template <typename T>
class SparseCOO :public Sparse<T> {

	std::vector<size_t> IA;//sizes of IA, A and JA are similar

	void setPortraitValues()
	{
		size_t index_i, index_j, curRow, diagIndex;
		T rowSum = 0;
		diagIndex = 0;
		curRow = 0;
		for (size_t i = 0; i < IA.size(); ++i)
		{
			index_i = IA[i];
			index_j = this->JA[i];
			if (curRow != index_i)//next row
			{
				this->A[diagIndex] = rowSum * 1.5;
				rowSum = 0;
				curRow = index_i;
			}
			if (index_i != index_j)
			{
				this->A[i] = sin(index_i + index_j);
				rowSum += abs(this->A[i]);
			}
			else
			{
				diagIndex = i;
			}	
		}
		this->A[diagIndex] = rowSum * 1.5;//last
	}

public:

	SparseCOO(const VariableSizeMeshContainer<int>& topoNN)
	{
		vector<int> temp;

		this->denseColumns = this->denseRows = topoNN.getBlockNumber();

		this->JA.reserve(topoNN.getBlockNumber());
		this->A.reserve(topoNN.getBlockNumber());

		for (size_t i = 0; i < topoNN.getBlockNumber(); ++i)
		{
			for (size_t j = 0; j < topoNN.getBlockSize(i); ++j)
			{
				temp.push_back(topoNN[i][j]);
			}
			temp.push_back(i);

			std::sort(temp.begin(), temp.end());

			for (size_t j = 0; j < topoNN.getBlockSize(i) + 1; ++j)
			{
				IA.push_back(i);
				this->JA.push_back(temp[j]);
				this->A.push_back(1);
			}
			temp.clear();
		}

		setPortraitValues();
	};

	T** getDenseMatrix() const override
	{
		T** dense = new T*[this->denseRows];
		for (size_t i = 0; i < this->denseRows; ++i)
		{
			dense[i] = new T[this->denseColumns];
			for (size_t j = 0; j < this->denseColumns; ++j)
			{
				dense[i][j] = 0;
			}
		}

		for (size_t i = 0; i < IA.size(); ++i)
		{
			dense[IA[i]][this->JA[i]] = this->A[i];
		}
		return dense;
	};

	void spmv(const std::vector<T>& x, std::vector<T>& result, size_t threadsNumber = 4) const override
	{
		if (x.size() != this->denseColumns)
		{
			std::cerr << "Incompatible sizes in spmv!" << std::endl;
		}
		else
		{
			for (size_t i = 0; i < result.size(); ++i) result[i] = 0;

			double start = omp_get_wtime();
			#pragma omp parallel for num_threads(threadsNumber)
			for (OPENMP_INDEX_TYPE i = 0; i < IA.size(); ++i)
			{
				result[IA[i]] += (this->A[i]) * x[this->JA[i]];
			}
			timeSpmv+= omp_get_wtime() - start;
		}
	}

	void spmv(const T* x, std::vector<T>& result, size_t xSize,  size_t threadsNumber = 4) const override
	{
		if (xSize != this->denseColumns)
		{
			std::cerr << "Incompatible sizes in spmv!" << std::endl;
		}
		else
		{
			for (size_t i = 0; i < xSize; ++i) result[i] = 0;

			double start = omp_get_wtime();
			#pragma omp parallel for num_threads(threadsNumber)
			for (OPENMP_INDEX_TYPE i = 0; i < IA.size(); ++i)
			{
				result[IA[i]] += (this->A[i]) * x[this->JA[i]];
			}
			timeSpmv+= omp_get_wtime() - start;
		}	
	};

	void printIa() const
	{
		std::cout << "IA vector: ";
		for (size_t i = 0; i < IA.size(); ++i)
		{
			std::cout << IA[i] << " ";
		}
		std::cout << std::endl;
	}

	void setValues(const std::vector<T>& a) override
	{
		if (a.size() != this->getValuesSize())
		{
			std::cerr << "Incompatibe sizes of values vectors!" << std::endl;
		}
		else
		{
			for (size_t i = 0; i < a.size(); i++)
			{
				this->A[i] = a[i];
			}
		}
	}

	size_t getValuesSize() const override
	{
		return this->A.size();
	}

	std::vector<T> getDiagonal() const override
	{
		std::vector<T> diagonal;
		diagonal.reserve(this->denseColumns);

		for(size_t i = 0; i < this->IA.size();++i)
		{
			if (this->IA[i] == this->JA[i])
			{
				diagonal.push_back(this->A[i]);
			}
		}
		return diagonal;
	}


};
