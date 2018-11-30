#pragma once

#include <mpi.h>
#include <string.h>
#include <cassert>

const unsigned int DEFAULT_CAPACITY = 2;
const float CAPACITY_MULTIPLIER = 2.0;

using ind = signed int;

namespace perc {

template <class T>
class MPIVector {
public:
    /* Constructors and Destructor */
    MPIVector();
    MPIVector(ind capacity);
    MPIVector(ind capacity, const T& value);
    MPIVector(const MPIVector<T>& vector);
    MPIVector(MPIVector<T>&& vector);
    ~MPIVector();

    /* Assignment */
    MPIVector<T>& operator=(const MPIVector<T>& other);
    MPIVector<T>& operator=(MPIVector<T>&& other);

    /* Capacity, Element and Data Access */
    ind size() const;
    ind capacity() const;
    T* data();
    T& operator[](ind index);

    /* Modifiers */
    void push_back(const T& value);
    void pop_back();
    void resize(ind newCapacity);

    int Send(int dest, int tag, MPI_Comm comm);

    template <class U>
    friend std::ostream& operator<<(std::ostream& os, const MPIVector<U>& vector) {
        os << "vector with " << vector.size() << " elements and " << vector.capacity()
           << " capacity " << std::endl;
        return os;
    }

    /* ----- Helper functions ---- */

    static T* extractDataFromPointer(void* pointer) {
        // Cast to ind-pointer so we can increment by one ind (that saves the number of elements
        // currently in the vector), then cast result to T pointer
        reinterpret_cast<T*>(static_cast<ind*>(pointer) + 1);
    }

    static ind extractSizeFromPointer(void* pointer) {
        // Size of the vector is saved in the first position
        return static_cast<ind*>(pointer)[0];
    }

    static ind calculatePointerSize(ind capacity) {
        // Memory to be allocated is space for desired capacity elements and one ind encoding the
        // fillstatus of the vector
        return sizeof(ind) + capacity * sizeof(T);
    }

    static void setSizeinPointer(void* pointer, ind newSize) {
        static_cast<ind*>(pointer)[0] = newSize;
    }

private:
    // first signed int for size, then number of size element of type T
    // allocated memory can hold capacity many elements
    void* SizeAndData;

    // Number of elements currently
    ind Capacity;
};

template <class T>
MPIVector<T>::MPIVector() : Capacity(DEFAULT_CAPACITY) {
    SizeAndData = malloc(calculatePointerSize(Capacity));
    setSizeinPointer(SizeAndData, 0);
};

// This does not do any default construction of type T, memory is garbage after construction
template <class T>
MPIVector<T>::MPIVector(ind capacity) : Capacity(capacity) {
    SizeAndData = malloc(calculatePointerSize(Capacity));
    setSizeinPointer(SizeAndData, Capacity);
};

template <class T>
MPIVector<T>::MPIVector(ind capacity, const T& value) : Capacity(capacity) {
    SizeAndData = malloc(calculatePointerSize(Capacity));
    setSizeinPointer(SizeAndData, Capacity);

    // Set every element to the given value
    T* vectorData = data();
    for (ind i = 0; i < Capacity; i++) {
        vectorData[i] = value;
    }
};

// Copy constructor
template <class T>
MPIVector<T>::MPIVector(const MPIVector<T>& vector) {
    Capacity = vector.capacity();
    SizeAndData = malloc(calculatePointerSize(Capacity));
    ind vectorSize = vector.size();
    setSizeinPointer(SizeAndData, vectorSize);

    // Overwrite the current data content
    T* currentData = data();
    T* vectorData = vector.data();
    memcpy(currentData, vectorData, vectorSize * sizeof(T));
    /* Equivalent to
    for (ind i = 0; i < vectorSize; i++) {
        currentData[i] = vectorData[i];
    } */
};

// Move copy constructor
template <class T>
MPIVector<T>::MPIVector(MPIVector<T>&& vector) : SizeAndData(nullptr), Capacity(0) {
    SizeAndData = std::move(vector.SizeAndData);
    Capacity = vector.capacity;
    vector.Capacity = 0;
};

template <class T>
MPIVector<T>::~MPIVector() {
    free(SizeAndData);
};

// Copy Assignment
template <class T>
MPIVector<T>& MPIVector<T>::operator=(const MPIVector& other) {
    // Detect self-assignment
    if (this != &other) {
        // Allocate memory for size and size number of T´s in other
        Capacity = other.capacity();
        void* tempSizeAndData = malloc(calculatePointerSize(Capacity));
        ind otherSize = other.size();
        setSizeinPointer(tempSizeAndData, otherSize);

        // Set data
        T* tempData = extractDataFromPointer(tempSizeAndData);
        T* otherData = other.data();
        memcpy(tempData, otherData, otherSize * sizeof(T));

        // Out with the old, in with the new (capacity already set)
        free(SizeAndData);
        SizeAndData = tempSizeAndData;
    }
    return *this;
};

// Move Assignment
template <class T>
MPIVector<T>& MPIVector<T>::operator=(MPIVector&& other) {
    if (this != &other) {
        SizeAndData = std::move(other.SizeAndData);
        Capacity = other.Capacity;
    }
    return *this;
};

template <class T>
ind MPIVector<T>::size() const {
    return extractSizeFromPointer(SizeAndData);
};

template <class T>
ind MPIVector<T>::capacity() const {
    return Capacity;
};

template <class T>
T* MPIVector<T>::data() {
    return extractDataFromPointer(SizeAndData);
};

template <class T>
T& MPIVector<T>::operator[](ind index) {
    assert(index < size() && "Index out of bounds.");
    return data()[index];
};

template <class T>
void MPIVector<T>::push_back(const T& value) {
    ind currSize = size();
    if (currSize == Capacity) {
        resize(Capacity * CAPACITY_MULTIPLIER);
    }
    data()[currSize] = value;
    setSizeinPointer(SizeAndData, currSize + 1);
};

template <class T>
void MPIVector<T>::pop_back() {
    ind currSize = size();
    assert(currSize > 0 && "Pop on empty array.");
    // Just decrement size
    setSizeinPointer(SizeAndData, currSize - 1);
    // For now omitted: Making the array smaller if number of elements in the vector falls below
    // 1/CAPACITY_MULTIPLIER * Capacity
};

template <class T>
void MPIVector<T>::resize(ind newCapacity) {
    // We already have the desired capacity
    if (newCapacity == Capacity) return;

    ind currentSize = size();
    assert(currentSize <= newCapacity && "Attempting to make vector smaller than its contents.");

    void* tempSizeAndData = malloc(calculatePointerSize(newCapacity));
    // Number of actual elements stays the same
    setSizeinPointer(tempSizeAndData, currentSize);

    // Copy data into the new vector
    T* currentData = data();
    T* tempData = extractDataFromPointer(tempSizeAndData);
    memcpy(tempData, currentData, currentSize * sizeof(T));

    free(SizeAndData);
    Capacity = newCapacity;
    SizeAndData = tempSizeAndData;
};

template <class T>
int MPIVector<T>::Send(int dest, int tag, MPI_Comm comm){};

}  // namespace perc