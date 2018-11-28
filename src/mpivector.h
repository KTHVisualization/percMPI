#pragma once

#include <mpi.h>
#include <string.h>

const unsigned int DEFAULT_CAPACITY = 2;
const float CAPACITY_MULTIPLIER = 2.0;

using ind = signed int;

namespace perc {

template <class T>
class MPIVector {
public:
    MPIVector();
    MPIVector(ind size);
    MPIVector(ind size, const T& value);
    MPIVector(const MPIVector<T>& vector);

    ~MPIVector();

    ind size() const;
    void resize(ind newCapacity);
    ind capacity() const;
    T* data();

    MPIVector<T>& operator=(MPIVector<T>&& other);
    T& operator[](ind index);

    void push_back(const T& value);
    void pop_back();

    int Send(int dest, int tag, MPI_Comm comm);

private:
    // first signed int for size, then number of size element of type T
    // allocated memory can hold capacity many elements
    void* SizeAndData;

    // Number of elements currently
    ind Capacity;
};

template <class T>
MPIVector<T>::MPIVector() : Capacity(DEFAULT_CAPACITY) {
    SizeAndData = malloc(sizeof(ind) + Capacity * sizeof(T));
    // Set size to 0
    static_cast<ind*>(SizeAndData)[0] = 0;
};

// This does not do any default construction of type T, memory is garbage after construction
template <class T>
MPIVector<T>::MPIVector(ind size) : Capacity(size) {
    // Allocate memory for size and size number of T´s
    SizeAndData = malloc(sizeof(ind) + Capacity * sizeof(T));
    static_cast<ind*>(SizeAndData)[0] = Capacity;
};

template <class T>
MPIVector<T>::MPIVector(ind size, const T& value) : Capacity(size) {
    // Allocate memory for size and size number of T´s
    SizeAndData = malloc(sizeof(ind) + Capacity * sizeof(T));
    static_cast<ind*>(SizeAndData)[0] = Capacity;

    // Set data
    T* thisData = this.data();
    for (ind i = 0; i < Capacity; i++) {
        thisData[i] = value;
    }
};

template <class T>
MPIVector<T>::MPIVector(const MPIVector<T>& vector) {
    // Allocate memory for size and size number of T´s
    Capacity = vector.capacity();
    SizeAndData = malloc(sizeof(ind) + Capacity * sizeof(T));
    ind vectorSize = vector.size();
    static_cast<ind*>(SizeAndData)[0] = vectorSize;

    // Set data
    T* thisData = this.data();
    T* vectorData = vector.data();
    memcpy(thisData, vectorData, vectorSize * sizeof(T));
    /* Equivalent to
    for (ind i = 0; i < vectorSize; i++) {
        thisData[i] = vectorData[i];
    } */
};

template <class T>
MPIVector<T>::~MPIVector() {
    free(SizeAndData);
};

template <class T>
ind MPIVector<T>::size() const {
    return static_cast<ind*>(SizeAndData)[0];
};

template <class T>
void MPIVector<T>::resize(ind newCapacity) {
    // We already have the desired capacity
    if (newCapacity == Capacity) return;

    ind currentSize = size();
    assert(currentSize <= newCapacity && "Attempting to make vector smaller than its contents.");

    void* tempSizeAndData = malloc(sizeof(ind) + newCapacity * sizeof(T));
    // Number of actual elements stays the same
    static_cast<ind*>(tempSizeAndData)[0] = currentSize;

    // Copy data into new vector
    T* thisData = data();
    T* tempData = reinterpret_cast<T*>(static_cast<ind*>(tempSizeAndData) + 1);
    memcpy(tempData, thisData, currentSize * sizeof(T));

    free(SizeAndData);
    Capacity = newCapacity;
    SizeAndData = tempSizeAndData;
};

template <class T>
ind MPIVector<T>::capacity() const {
    return Capacity;
};

template <class T>
T* MPIVector<T>::data() {
    // Cast to ind-pointer so we can increment by one ind (the size)
    // then cast result to T pointer
    return reinterpret_cast<T*>(static_cast<ind*>(SizeAndData) + 1);
};

template <class T>
MPIVector<T>& MPIVector<T>::operator=(MPIVector&& other) {
    // Allocate memory for size and size number of T´s in other
    Capacity = other.capacity();
    void* tempSizeAndData = malloc(sizeof(ind) + Capacity * sizeof(T));
    ind otherSize = other.size();
    static_cast<ind*>(tempSizeAndData)[0] = otherSize;

    // Set data
    T* tempData = reinterpret_cast<T*>(static_cast<ind*>(tempSizeAndData) + 1);
    T* otherData = other.data();
    memcpy(tempData, otherData, otherSize * sizeof(T));

    // Out with the old, in with the new (capacity already set)
    free(SizeAndData);
    SizeAndData = tempSizeAndData;
    return *this;
};

template <class T>
T& MPIVector<T>::operator[](ind index) {
    assert(index < size() && "index out of bounds.");
    return data()[index];
};

template <class T>
void MPIVector<T>::push_back(const T& value) {
    ind currSize = static_cast<ind*>(SizeAndData)[0];
    if (currSize == Capacity) {
        resize(Capacity * CAPACITY_MULTIPLIER);
    }
    data()[currSize] = value;
    static_cast<ind*>(SizeAndData)[0] = currSize + 1;
};

template <class T>
void MPIVector<T>::pop_back() {
    // Just decrement size
    ind& currSize = static_cast<ind*>(SizeAndData)[0];
    assert(currSize > 0 && "pop on empty array.");
    currSize--;
    // For now omitted: Making the array smaller if size falls below 1/CAPACITY_MULTIPLIER *
    // Capacity
};

template <class T>
int MPIVector<T>::Send(int dest, int tag, MPI_Comm comm){};

}  // namespace perc