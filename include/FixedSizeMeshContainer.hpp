#pragma once
#include "stdafx.hpp"

#include <fstream>

using namespace std;

// Realizing one-dimensional effective-accessible and storing block vector with
// matrix interface.
template <typename T> class FixedSizeMeshContainer {

    size_t blockSize;
    vector<T> V; // main grid vector

public:
    FixedSizeMeshContainer();

    FixedSizeMeshContainer(size_t num);

    FixedSizeMeshContainer(const vector<T> &source, size_t blockSize);

    FixedSizeMeshContainer(const FixedSizeMeshContainer<T> &source);

    void setBlockSize(size_t newSize);

    // Can only add elements of vectors which sizes are divided without remainder
    // into {M}
    bool add(const vector<T> &extra);

    // Add another GridContainer elements to our container
    bool add(const FixedSizeMeshContainer<T> &extra);

    void reserve(size_t);

    void shrink_to_fit();

    FixedSizeMeshContainer &operator=(const FixedSizeMeshContainer<T> &target);

    FixedSizeMeshContainer operator+(const FixedSizeMeshContainer<T> &extra);

    FixedSizeMeshContainer &operator+=(const FixedSizeMeshContainer<T> &extra);

    T *operator[](size_t i);

    const T *operator[](size_t i) const;

    size_t getBlockNumber() const;

    size_t getBlockSize() const;

    size_t getTotalSize() const;

    void inline printContainer(ostream &out = cout) const {
        for (size_t i = 0; i < getBlockNumber(); ++i) {
            for (size_t j = 0; j < blockSize; ++j) {
                if (j % 2 == 0) {
                    out << operator[](i)[j];
                    out << " ";
                } else {
                    out << operator[](i)[j];
                }
            }
            out << endl;
        }
    }

    void inline printContainer(ofstream &out) const {

        for (size_t i = 0; i < getBlockNumber(); ++i) {
            for (size_t j = 0; j < blockSize; ++j) {
                out << operator[](i)[j];
                out << " ";
            }
            out << 0 << " "; // z - coordinate(for paraview)
            out << " ";
        }

        out << endl;
    }

    void clear();
};

// default constructor
template <typename T> FixedSizeMeshContainer<T>::FixedSizeMeshContainer() { blockSize = 0; }

template <typename T> FixedSizeMeshContainer<T>::FixedSizeMeshContainer(size_t num) : blockSize(num){};

template <typename T>
FixedSizeMeshContainer<T>::FixedSizeMeshContainer(const vector<T> &source, size_t blockSize) : blockSize(blockSize) {
    V.reserve(source.size());
    for (typename vector<T>::const_iterator it = source.begin(); it != source.end(); ++it) {
        V.push_back(*it);
    }
}

template <typename T>
FixedSizeMeshContainer<T>::FixedSizeMeshContainer(const FixedSizeMeshContainer<T> &source)
    : blockSize(source.blockSize) {
    V.resize(source.blockSize * source.getBlockNumber()); // allocating

    V = source.V;
}

template <typename T> void FixedSizeMeshContainer<T>::setBlockSize(size_t newSize) { blockSize = newSize; }

// Can only add elements of vectors which sizes are divided without remainder
// into {M}
template <typename T> bool FixedSizeMeshContainer<T>::add(const vector<T> &extra) {
    if (extra.size() % blockSize != 0) {
        cerr << "Error: Vector size must be a multiple of " << blockSize << endl;
        return false;
    } else {
        V.reserve(extra.size());
        for (typename vector<T>::const_iterator it = extra.begin(); it != extra.end(); ++it) {
            V.push_back(*it);
        }
        return true;
    }
}

// Add another GridContainer elements to our container
template <typename T> bool FixedSizeMeshContainer<T>::add(const FixedSizeMeshContainer<T> &extra) {
    if (extra.blockSize != blockSize) {
        cerr << "Error: Incompatible block sizes!" << endl;
        return false;
    } else {
        size_t odlBlocksNum = getBlockNumber();
        V.resize((getBlockNumber() + extra.getBlockNumber()) * blockSize);
        for (size_t i = odlBlocksNum; i < getBlockNumber(); ++i) {
            for (size_t j = 0; j < blockSize; ++j) {
                operator[](i)[j] = extra[i - odlBlocksNum][j];
            }
        }
        return true;
    }
}

template <typename T> void FixedSizeMeshContainer<T>::reserve(size_t size) { V.reserve(size); }

template <typename T> void FixedSizeMeshContainer<T>::shrink_to_fit() { V.shrink_to_fit(); }

template <typename T>
FixedSizeMeshContainer<T> &FixedSizeMeshContainer<T>::operator=(const FixedSizeMeshContainer<T> &extra) {
    blockSize = extra.blockSize;
    V.clear();
    V.resize(extra.blockSize * extra.BlocksNumber);
    V = extra.V;
    return (*this);
}

template <typename T>
FixedSizeMeshContainer<T> FixedSizeMeshContainer<T>::operator+(const FixedSizeMeshContainer<T> &extra) {
    FixedSizeMeshContainer temp_grid(*this);
    temp_grid.add(extra);
    return temp_grid;
}

template <typename T>
FixedSizeMeshContainer<T> &FixedSizeMeshContainer<T>::operator+=(const FixedSizeMeshContainer<T> &extra) {
    (*this).add(extra);
    return *this;
}

template <typename T> T *FixedSizeMeshContainer<T>::operator[](size_t i) {
    if (i >= getBlockNumber()) {
        cerr << "Error: Index out of bounds" << endl;
        return nullptr;
    } else {
        return &V[0] + blockSize * i;
    }
}

template <typename T> const T *FixedSizeMeshContainer<T>::operator[](size_t i) const {
    if (i >= getBlockNumber()) {
        cerr << "Error: Index out of bounds" << endl;
        return nullptr;
    } else {
        return &V[0] + blockSize * i;
    }
}

template <typename T> size_t FixedSizeMeshContainer<T>::getBlockNumber() const { return V.size() / blockSize; }

template <typename T> size_t FixedSizeMeshContainer<T>::getBlockSize() const { return blockSize; }

template <typename T> size_t FixedSizeMeshContainer<T>::getTotalSize() const { return V.size(); }

template <typename T> void FixedSizeMeshContainer<T>::clear() { V.clear(); }
