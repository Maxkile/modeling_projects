#pragma once

#include "VariableSizeMeshContainer.hpp"
#include <vector>
#include <omp.h>
#include <algorithm>
#include <iomanip>


template<typename T>
class Sparse
{
protected:

    size_t denseRows,denseColumns;
    std::vector<int> JA;
    std::vector<T> A;

public:

    virtual std::vector<T> spmv(const std::vector<T>&,size_t) const = 0;
    virtual T* spmv(const T*,size_t,size_t) const = 0;
    virtual T** getDenseMatrix() const = 0;

    inline void printA() const 
    {
        std::cout << "A vector: ";
        for(size_t i = 0;i < A.size();++i)
        {
            std::cout << A[i] << " ";
        }
        std::cout << std::endl;
    }

    inline void printJa() const 
    {
        std::cout << "JA vector: ";
        for(size_t i = 0;i < JA.size();++i)
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
        for (int i = 0; i < denseRows; i++) {

            for (int j = 0; j < denseColumns; j++)
                std::cout << fixed << setprecision(2) << matr[i][j] << " ";
            std::cout << endl;
        }
        std::cout << endl;
    }
	
};

// template<typename T, typename = std::enable_if_t<std::is_same<T,int>::value>>
template<typename T>
class SparseELL: public Sparse<T>
{
public:
    
    SparseELL(const VariableSizeMeshContainer<int>& topoNN,const std::vector<T>& nodeValues)
    {
        this->denseRows = this->denseColumns = topoNN.getBlockNumber();
        int** m_portrait = build_portrait(topoNN);

        for(size_t i = 0;i < topoNN.getBlockNumber();++i)
        {
            for(size_t j = 0;j < topoNN.getBlockSize(i);++j)
            {
                m_portrait[i][j] = 1;
                m_portrait[j][i] = 1;
            }
        }

        findRowOffset(m_portrait,this->denseRows,this->denseRows);
        buildAJA(m_portrait,this->denseRows,this->denseRows);
        
        for(size_t i = 0;i < topoNN.getBlockNumber();++i)
        {
            delete[] m_portrait[i];
        }
        delete[] m_portrait;
    }

    T** getDenseMatrix() const override
    {
        T** dense = new T*[this->denseRows];
        for(size_t i = 0;i < this->denseRows;++i)
        {
            dense[i] = new T[this->denseColumns];
            for(size_t j = 0;j < this->denseColumns;++j)
            {
                dense[i][j] = 0;
            }
        }

        for(size_t i = 0;i < this->denseRows;++i)
        {
            for(size_t j = i * this->rowOffset;j < (i + 1) * this->rowOffset;++j)
            {
                if (this->A[j] != 0)//padding
                {
                    dense[i][this->JA[j]] = this->A[j];
                }
            }
        }
        return dense;   
    }

    std::vector<T> spmv(const std::vector<T>& x,size_t threadsNumber = 4) const override
    {
        std::vector<T> y;
        if (x.size() != this->denseColumns)
        {
            std::cerr << "Incompatible sizes!" << std::endl;
            return y;
        }
        else
        {
            y.resize(x.size());

            #pragma omp parallel for num_threads(threadsNumber)//TODO: ask about template reduction
            for(size_t i = 0;i < this->denseRows;++i)
            {
                y[i] = 0;

                // #pragma omp parallel for reduction(+:result[i])
                for(size_t j = i * this->rowOffset;j < (i + 1) * this->rowOffset;++j)
                {
                    y[i] += this->A[j] * x[this->JA[j]];
                }
            }
            return y;
        }
    }

    T* spmv(const T* x,size_t xSize,size_t threadsNumber = 4) const override
    {
        
        if (xSize != this->denseColumns)
        {
            std::cerr << "Incompatible sizes!" << std::endl;
            return nullptr;
        }
        else
        {
            T* y = new T[xSize];//need to deallocate then
            for(size_t i = 0;i < this->denseRows;++i)
            {
                y[i] = 0;

                #pragma omp parallel for num_threads(threadsNumber) reduction(+:y[i])
                for(size_t j = i * rowOffset;j < (i + 1) * rowOffset;++j)
                {
                    y[i] += this->A[j] * x[this->JA[j]];
                }
            }
            return y;
        }
    }

    size_t getRowOffset() const
    {
        return rowOffset;
    }

private:
    size_t rowOffset;
    
    void findRowOffset(T** dense,size_t rows,size_t columns)
    {
        this->rowOffset = 1;
        size_t currentRowOffset;
        
        for(size_t i = 0; i < rows;++i)
        {
            currentRowOffset = 0;
            for(size_t j = 0; j < columns; ++j)
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

    void buildAJA(T** dense,size_t rows,size_t columns)
    {
        int rowElemCount;

        this->JA.reserve(this->rowOffset * rows);
        this->A.reserve(this->rowOffset * rows);

        for(size_t i = 0;i < rows;++i)
        {
            rowElemCount = 0;
            for(size_t j = 0;j < columns;++j)
            {
                if (dense[i][j] != 0)
                {
                    this->A.push_back(dense[i][j]);
                    this->JA.push_back(j);
                    rowElemCount++;
                }
            }
            size_t padding_ind = this->JA[this->JA.size() - 1];
            for(size_t count = rowElemCount;count < this->rowOffset;++count)
            {
                this->A.push_back(0);
                this->JA.push_back(padding_ind);
            }
        }
    }

    int** build_portrait(const VariableSizeMeshContainer<int>& topoNN) const//not forget to free memory lately
    {
        size_t blockNumber = topoNN.getBlockNumber();
        
        int** portrait = new int*[blockNumber];
        for(size_t i = 0; i < blockNumber; ++i)
        {
            portrait[i] = new int[blockNumber];
            for(size_t j = 0; j < blockNumber; ++j)
            {
                portrait[i][j] = 0;
            }
        }

        for(size_t i = 0;i < topoNN.getBlockNumber();++i)
        {
            portrait[i][i] = 1;//diagonal
            for(size_t j = 0;j < topoNN.getBlockSize(i);++j)
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
class SparseCSR:public Sparse<T>{

    std::vector<int> IA;

public:

    SparseCSR(const VariableSizeMeshContainer<int>& topoNN)
    {
        std::vector<int> temp;
        this->denseColumns = this->denseRows = topoNN.getBlockNumber();
        int prev;

        IA.push_back(0);
        for (int i = 0; i < topoNN.getBlockNumber(); i++){
            temp.push_back(i);
            for (int j = 0; j < topoNN.getBlockSize(i); j++){
                temp.push_back(topoNN[i][j]);
            }
            std::sort(temp.begin(),temp.end());
            for (int j = 0; j < temp.size(); j++)
                Sparse<T>::JA.push_back(temp[j]);
            prev = IA[IA.size() - 1];
            IA.push_back(temp.size()+prev);
            temp.clear();
        }
        IA.push_back(this->denseColumns);
    };

    T** getDenseMatrix() const override
    {
        T** dense = new T*[this->denseRows];
        for(size_t i = 0;i < this->denseRows;++i) {
            dense[i] = new T[this->denseColumns];
            for(size_t j = 0;j < this->denseColumns;++j)
            {
                dense[i][j] = 0;
            }
        }

        for(size_t i = 0; i < this->denseRows; ++i) {
            for (size_t j = IA[i]; j < IA[i+1]; ++j){
                dense[i][this->JA[j]] = this->A[j];
            }
        }
        return dense;
    };

    std::vector<T> spmv(const std::vector<T>& x,size_t threadsNumber = 4) const override
    {
        std::vector<T> result;

        if (x.size() != this->denseColumns)
        {
            std::cerr << "Incompatible sizes!" << std::endl;
            return result;
        }

        result.reserve(this->denseRows);
        
        
        omp_set_num_threads(threadsNumber);

        #pragma omp parallel for
        for (int i = 0; i < this->denseRows; i++){
            
            result[i] = 0;    

            // #pragma omp parallel for reduction(+:result[i])
            for (int j = IA[i]; j < IA[i+1]; j++)
                result[i] += x[this->JA[j]] * (this->A[j]);
        }

        return result;
    };

    T* spmv(const T* x,size_t xSize,size_t threadsNumber = 4) const override
    {
        T* result = new T[xSize];
        if (xSize != this->denseColumns)
        {
            std::cerr << "Incompatible sizes!" << std::endl;
            return nullptr;
        }

        for(int i = 0; i < this->denseRows;++i) {

            result[i] = 0;

			#pragma omp parallel for num_threads(threadsNumber) reduction(+:result[i])
            for (int j = IA[i]; j < IA[i+1]; j++)
                result[i] += x[this->JA[j]] * (this->A[j]);
        }

        return result;
    }

    void set_values(const std::vector<T>& a)
    {
        for (int i = 0; i < a.size(); i++)
            (this->A).push_back(a[i]);
    }

    void printIa() const
    {
        std::cout << "IA vector: ";
        for(size_t i = 0;i < IA.size();++i)
        {
            std::cout << IA[i] << " ";
        }
        std::cout << std::endl;
    }

};

template <typename T>
class SparseCOO:public Sparse<T>{

    std::vector<int> IA;

public:

    SparseCOO(const VariableSizeMeshContainer<int>& topoNN)
    {
        this->denseColumns = this->denseRows = topoNN.getBlockNumber();

        for (int i = 0; i < topoNN.getBlockNumber(); ++i){
            for (int j = 1; j < topoNN.getBlockSize(i); ++j){
                IA.push_back(topoNN[i][0]);
                (this->JA).push_back(topoNN[i][j]);
            }
        }
    };

    T** getDenseMatrix() const override
    {
        T** dense = new T*[this->denseRows];
        for(int i = 0;i < this->denseRows;++i)
        {
            dense[i] = new T[this->denseColumns];
            for(int j = 0;j < this->denseColumns;++j)
            {
                dense[i][j] = 0;
            }
        }

        for(int i = 0; i < IA.size(); ++i) {
            dense[IA[i]][this->JA[i]] = this->A[i];
        }
        return dense;
    };

    std::vector<T> spmv(const std::vector<T>& x,size_t threadsNumber = 4) const override
    {
        std::vector<T> result;

        if (x.size() != this->denseColumns)
        {
            std::cerr << "Incompatible sizes!" << std::endl;
            return result;
        }

        result.reserve(this->denseRows);
        for (int i = 0; i < result.size(); ++i) result[i] = 0;

        #pragma omp parallel for num_threads(threadsNumber) reduction(+:result[i])
        for (int i = 0; i < IA.size(); ++i)
        {
            result[IA[i]] += (this->A[i]) * x[this->JA[i]];
        }

        return result;
    };

    T* spmv(const T* x,size_t xSize,size_t threadsNumber = 4) const override
    {
        T* result = new T[xSize];

        if (xSize != this->denseColumns)
        {
            std::cerr << "Incompatible sizes!" << std::endl;
            return nullptr;
        }

        for(int i = 0; i < result.size(); ++i) result[i] = 0;

        #pragma omp parallel for num_threads(threadsNumber) reduction(+:result[i])
        for (int i = 0; i < IA.size(); ++i)
        {
            result[IA[i]] += (this->A[i]) * x[this->JA[i]];
        }

        return result;
    };

    void printIa() const
    {
        std::cout << "IA vector: ";
        for(size_t i = 0;i < IA.size();++i)
        {
            std::cout << IA[i] << " ";
        }
        std::cout << std::endl;
    }
    
    void set_values(const std::vector<T>& a)
    {
        for (int i = 0; i < a.size(); i++)
            (this->A).push_back(a[i]);
    }
};
