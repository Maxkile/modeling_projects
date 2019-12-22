#pragma once
#include "stdafx.hpp"

#include <omp.h>
#include <algorithm>
#include <iomanip>
#include <cmath>
#include "platformDependencies.hpp"
#include "VariableSizeMeshContainer.hpp"

template<typename T>
class Sparse
{
protected:

	size_t denseRows, denseColumns;
	std::vector<size_t> JA;
	std::vector<T> A;

public:

	virtual std::vector<T> spmv(const std::vector<T>&, size_t threadsNumber = 4) const = 0;
	virtual T* spmv(const T*, size_t, size_t threadsNumber = 4) const = 0;
	virtual T** getDenseMatrix() const = 0;
	virtual void setValues(const std::vector<T>&) = 0;
	virtual size_t getValuesSize() const = 0;
	virtual void set_model_values() = 0;
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
	
};

template<typename T>
class SparseELL : public Sparse<T>
{
public:

	SparseELL(T** dense, size_t rows, size_t columns)//from dense(needed in solver)
	{
		this->valuesSize = 0;
		this->denseRows = rows;
		this->denseColumns = columns;
		findRowOffset(dense, this->denseRows, this->denseRows);
		buildAJA(dense, this->denseRows, this->denseRows);
	}

	SparseELL(const VariableSizeMeshContainer<int>& topoNN)//from topo
	{
		this->valuesSize = 0;
		this->denseRows = this->denseColumns = topoNN.getBlockNumber();
		T** m_portrait = build_portrait(topoNN);

		findRowOffset(m_portrait, this->denseRows, this->denseRows);
		buildAJA(m_portrait, this->denseRows, this->denseRows);

		for (size_t i = 0; i < topoNN.getBlockNumber(); ++i)
		{
			delete[] m_portrait[i];
		}
		delete[] m_portrait;
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

	std::vector<T> spmv(const std::vector<T>& x, size_t threadsNumber = 4) const override
	{
		std::vector<T> y(x.size());
		if (x.size() != this->denseColumns)
		{
			std::cerr << "Incompatible sizes in spmv!" << std::endl;
			return y;
		}
		else
		{
			#pragma omp parallel for num_threads(threadsNumber)//TODO: ask about template reduction
			for (OPENMP_INDEX_TYPE i = 0; i < this->denseRows; ++i)
			{
				y[i] = 0;

				for (OPENMP_INDEX_TYPE j = i * this->rowOffset; j < (i + 1) * this->rowOffset; ++j)
				{
					y[i] += x[this->JA[j]] * this->A[j];
				}
			}

			return y;
		}
	}

	T* spmv(const T* x, size_t xSize, size_t threadsNumber = 4) const override
	{

		if (xSize != this->denseColumns)
		{
			std::cerr << "Incompatible sizes in spmv!" << std::endl;
			return nullptr;
		}
		else
		{
			T* y = new T[xSize];//need to deallocate then

			#pragma omp parallel for num_threads(threadsNumber) 
			for (OPENMP_INDEX_TYPE i = 0; i < this->denseRows; ++i)
			{
				y[i] = 0;

				for (OPENMP_INDEX_TYPE j = i * rowOffset; j < (i + 1) * rowOffset; ++j)
				{
					y[i] += x[this->JA[j]] * this->A[j];
				}
			}
			return y;
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
				if (this->JA[j] == i - 1)//diag elem
				{
					diagonal.push_back(this->A[j]);
				}
			}
		}
		return diagonal;		
	}

	void set_model_values() override {
		T** matr = getDenseMatrix();

		int k = 0; // number non zero element
		std::vector<int> diag;
		std::vector<T> values((*this).getValuesSize());
		double sum;
		
		// set non diag elements by sin(i+j)
		for (size_t i = 0; i < this->denseRows; i++)
			for (size_t j = 0; j < this->denseColumns; j++)
				if (matr[i][j] != 0){
					if (i != j)
						matr[i][j] = values[k] = sin(i + j);
					else
						diag.push_back(k);
					k++;
				}
					
		// set diag elements by (sum(abs(mat[i][j]))) * 1.1
		for (size_t i = 0; i < this->denseRows; i++){
			sum = 0;
			for (size_t j = 0; j < this->denseColumns; j++)
				if (j != i)
					sum += abs(matr[i][j]);
			sum *= 1.1;
			values[diag[i]] = sum;
		}

		for (size_t i = 0; i < this->denseRows; ++i)
		{
			delete[] matr[i];
		}
		delete[] matr;
		
		(*this).setValues(values);
	}
	
private:
	size_t rowOffset, valuesSize;

	void findRowOffset(T** dense, size_t rows, size_t columns)
	{
		this->rowOffset = 1;
		size_t currentRowOffset;

		for (size_t i = 0; i < rows; ++i)
		{
			currentRowOffset = 0;
			for (size_t j = 0; j < columns; ++j)
			{
				if (dense[i][j] != 0)
				{
					currentRowOffset++;
				}
			}
			if (currentRowOffset > this->rowOffset)
			{
				this->rowOffset = currentRowOffset;
			}
		}
	}

	void buildAJA(T** dense, size_t rows, size_t columns)
	{
		int rowElemCount;

		this->JA.reserve(this->rowOffset * rows);
		this->A.reserve(this->rowOffset * rows);

		for (size_t i = 0; i < rows; ++i)
		{
			rowElemCount = 0;
			for (size_t j = 0; j < columns; ++j)
			{
				if (dense[i][j] != 0)
				{
					this->valuesSize++;
					this->A.push_back(dense[i][j]);
					this->JA.push_back(j);
					rowElemCount++;
				}
			}
			size_t padding_ind = this->JA[this->JA.size() - 1];
			for (size_t count = rowElemCount; count < this->rowOffset; ++count)
			{
				this->A.push_back(0);
				this->JA.push_back(padding_ind);
			}
		}
	}

	T** build_portrait(const VariableSizeMeshContainer<int>& topoNN) const//not forget to free memory lately
	{
		size_t blockNumber = topoNN.getBlockNumber();

		T** portrait = new T*[blockNumber];
		for (size_t i = 0; i < blockNumber; ++i)
		{
			portrait[i] = new T[blockNumber];
			for (size_t j = 0; j < blockNumber; ++j)
			{
				portrait[i][j] = 0;
			}
		}

		for (size_t i = 0; i < blockNumber; ++i)
		{
			portrait[i][i] = 1;//diagonal
			for (size_t j = 0; j < topoNN.getBlockSize(i); ++j)
			{
				size_t _j = topoNN[i][j];
				portrait[i][_j] = 1;
				portrait[_j][i] = 1;
			}
		}
		return portrait;
	}
	
};

template <typename T>
class SparseCSR :public Sparse<T> {

	std::vector<size_t> IA;//sizes of A and JA are similar

public:

	SparseCSR(const VariableSizeMeshContainer<int>& topoNN)
	{
		std::vector<int> temp;
		this->denseColumns = this->denseRows = topoNN.getBlockNumber();
		int prev;

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
				Sparse<T>::JA.push_back(temp[j]);
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

		for (size_t i = 0; i < this->denseRows; ++i)
		{
			for (size_t j = IA[i]; j < IA[i + 1]; ++j)
			{
				dense[i][this->JA[j]] = this->A[j];
			}
		}
		return dense;
	};

	std::vector<T> spmv(const std::vector<T>& x, size_t threadsNumber = 4) const override
	{
		std::vector<T> result(x.size());

		if (x.size() != this->denseColumns)
		{
			std::cerr << "Incompatible sizes in spmv!" << std::endl;
			return result;
		}

		#pragma omp parallel for num_threads(threadsNumber)
		for (OPENMP_INDEX_TYPE i = 0; i < this->denseRows; i++)
		{

			result[i] = 0;

			for (size_t j = IA[i]; j < IA[i + 1]; j++)
				result[i] += x[this->JA[j]] * (this->A[j]);
		}

		return result;
	};

	T* spmv(const T* x, size_t xSize, size_t threadsNumber = 4) const override
	{
		T* result = new T[xSize];
		if (xSize != this->denseColumns)
		{
			std::cerr << "Incompatible sizes in spmv!" << std::endl;
			return nullptr;
		}

		#pragma omp parallel for num_threads(threadsNumber)
		for (OPENMP_INDEX_TYPE i = 0; i < this->denseRows; ++i) {

			result[i] = 0;

			for (size_t j = IA[i]; j < IA[i + 1]; j++)
				result[i] += x[this->JA[j]] * (this->A[j]);
		}

		return result;
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
	
	void set_model_values() override {
		std::cout << "Rows: " << this->denseRows << "Columns: " << this->denseColumns << std::endl;
		T** matr = getDenseMatrix();

		int k = 0; // number non zero element
		std::vector<int> diag;
		double sum;
		
		// set non diag elements by sin(i+j)
		for (size_t i = 0; i < this->denseRows; i++)
			for (size_t j = 0; j < this->denseColumns; j++)
				if (matr[i][j] != 0){
					if (i != j)
						matr[i][j] = this->A[k] = sin(i + j);
					else
						diag.push_back(k);
					k++;
				}
	
		// set diag elements by (sum(abs(mat[i][j]))) * 1.1
		for (size_t i = 0; i < this->denseRows; i++){
			sum = 0;
			for (size_t j = 0; j < this->denseColumns; j++)
				if (j != i)
					sum += abs(matr[i][j]);
			sum *= 1.1;
			this->A[diag[i]] = sum;
		}

		for(size_t i = 0;i < this->denseRows;++i)
		{
			delete[] matr[i];
		}
		delete[] matr;
	}
};

template <typename T>
class SparseCOO :public Sparse<T> {

	std::vector<size_t> IA;//sizes of IA, A and JA are similar

public:

	SparseCOO(const VariableSizeMeshContainer<int>& topoNN)
	{
		vector<int> temp;

		this->denseColumns = this->denseRows = topoNN.getBlockNumber();

		for (size_t i = 0; i < topoNN.getBlockNumber(); ++i)
		{
			for (size_t j = 0; j < topoNN.getBlockSize(i); j++)
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

	std::vector<T> spmv(const std::vector<T>& x, size_t threadsNumber = 4) const override
	{
		std::vector<T> result(x.size());

		if (x.size() != this->denseColumns)
		{
			std::cerr << "Incompatible sizes in spmv!" << std::endl;
			return result;
		}

		#pragma omp parallel for num_threads(threadsNumber)
		for (OPENMP_INDEX_TYPE i = 0; i < IA.size(); ++i)
		{
			result[IA[i]] += (this->A[i]) * x[this->JA[i]];
		}

		return result;
	};

	T* spmv(const T* x, size_t xSize, size_t threadsNumber = 4) const override
	{
		T* result = new T[xSize];

		if (xSize != this->denseColumns)
		{
			std::cerr << "Incompatible sizes in spmv!" << std::endl;
			return nullptr;
		}

		for (size_t i = 0; i < xSize; ++i) result[i] = 0;

		#pragma omp parallel for num_threads(threadsNumber)
		for (OPENMP_INDEX_TYPE i = 0; i < IA.size(); ++i)
		{
			result[IA[i]] += (this->A[i]) * x[this->JA[i]];
		}

		return result;
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

	void set_model_values() override {
		T** matr = getDenseMatrix();

		int k = 0; // number non zero element
		std::vector<int> diag;
		double sum;
		
		// set non diag elements by sin(i+j)
		for (size_t i = 0; i < this->denseRows; i++)
			for (size_t j = 0; j < this->denseColumns; j++)
				if (matr[i][j] != 0){
					if (i != j)
						matr[i][j] = this->A[k] = sin(i + j);
					else
						diag.push_back(k);
					k++;
				}
	
		// set diag elements by (sum(abs(mat[i][j]))) * 1.1
		for (size_t i = 0; i < this->denseRows; i++){
			sum = 0;
			for (size_t j = 0; j < this->denseColumns; j++)
				if (j != i)
					sum += abs(matr[i][j]);
			sum *= 1.1;
			this->A[diag[i]] = sum;
		}

		for (size_t i = 0; i < this->denseRows; ++i)
		{
			delete[] matr[i];
		}
		delete[] matr;
	}

};
