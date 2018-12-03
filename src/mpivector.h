#pragma once

#include <mpi.h>
#include <string.h>
#include <cassert>

const unsigned int DEFAULT_CAPACITY = 16;
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
    const T* data() const;
    T& operator[](ind index);

    /* Modifiers */
    void push_back(const T& value);
    void pop_back();
    void resize(ind newCapacity);

    int Send(int dest, int tag, MPI_Comm comm);
    int Recv(int src, int tag, MPI_Comm comm, MPI_Status* status);

    friend std::ostream& operator<<(std::ostream& os, const MPIVector& vector) {
        os << "vector with " << vector.size() << " elements and " << vector.capacity()
           << " capacity " << std::endl;
        return os;
    }

    void printContents() const;

    /* ----- Helper functions ---- */

    static T* extractDataFromPointer(void* pointer) {
        // Cast to ind-pointer so we can increment by one ind (that saves the number of elements
        // currently in the vector), then cast result to T pointer
        return reinterpret_cast<T*>(reinterpret_cast<ind*>(pointer) + 1);
    }

    static ind extractSizeFromPointer(void* pointer) {
        // Size of the vector is saved in the first position
        return reinterpret_cast<ind*>(pointer)[0];
    }

    static ind calculatePointerSize(ind capacity) {
        // Memory to be allocated is space for desired capacity elements and one ind encoding the
        // fillstatus of the vector
        return sizeof(ind) + capacity * sizeof(T);
    }

    static void setSizeinPointer(void* pointer, ind newSize) {
        reinterpret_cast<ind*>(pointer)[0] = newSize;
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
}

// This does not do any default construction of type T, memory is garbage after construction
template <class T>
MPIVector<T>::MPIVector(ind capacity) : Capacity(capacity) {
    SizeAndData = malloc(calculatePointerSize(Capacity));
    setSizeinPointer(SizeAndData, Capacity);
}

template <class T>
MPIVector<T>::MPIVector(ind capacity, const T& value) : Capacity(capacity) {
    SizeAndData = malloc(calculatePointerSize(Capacity));
    setSizeinPointer(SizeAndData, Capacity);

    // Set every element to the given value
    T* vectorData = data();
    for (ind i = 0; i < Capacity; i++) {
        vectorData[i] = value;
    }
}

// Copy constructor
template <class T>
MPIVector<T>::MPIVector(const MPIVector<T>& vector) {
#ifdef TEST
    std::cout << "Called copy construction." << std::endl;
#endif
    Capacity = vector.capacity();
    SizeAndData = malloc(calculatePointerSize(Capacity));
    ind vectorSize = vector.size();
    setSizeinPointer(SizeAndData, vectorSize);

    // Overwrite the current data content
    T* currentData = data();
    const T* vectorData = vector.data();
    memcpy(currentData, vectorData, vectorSize * sizeof(T));
    /* Equivalent to
    for (ind i = 0; i < vectorSize; i++) {
        currentData[i] = vectorData[i];
    } */
}

// Move copy constructor
template <class T>
MPIVector<T>::MPIVector(MPIVector<T>&& vector) : SizeAndData(nullptr), Capacity(0) {
#ifdef TEST
    std::cout << "Called move constructor." << std::endl;
#endif
    SizeAndData = vector.SizeAndData;
    vector.SizeAndData = nullptr;

    Capacity = vector.Capacity;
    vector.Capacity = 0;
}

template <class T>
MPIVector<T>::~MPIVector() {
#ifdef TEST
    std::cout << "Destructing";
    if (SizeAndData)
        std::cout << " a " << *this;
    else
        std::cout << " a vector that holds no data anymore." << std::endl;
    // Could be a null pointer due to move construction
    if (SizeAndData) free(SizeAndData);
#endif
}

// Copy Assignment
template <class T>
MPIVector<T>& MPIVector<T>::operator=(const MPIVector& other) {
#ifdef TEST
    std::cout << "Called copy assignment." << std::endl;
#endif
    // Detect self-assignment
    if (this != &other) {
        // Allocate memory for size and size number of T´s in other
        Capacity = other.capacity();
        void* tempSizeAndData = malloc(calculatePointerSize(Capacity));
        ind otherSize = other.size();
        setSizeinPointer(tempSizeAndData, otherSize);

        // Set data
        T* tempData = extractDataFromPointer(tempSizeAndData);
        const T* otherData = other.data();
        memcpy(tempData, otherData, otherSize * sizeof(T));

        // Out with the old, in with the new (capacity already set)
        free(SizeAndData);
        SizeAndData = tempSizeAndData;
    }
    return *this;
}

// Move Assignment
template <class T>
MPIVector<T>& MPIVector<T>::operator=(MPIVector&& other) {
#ifdef TEST
    std::cout << "Called move assignment." << std::endl;
#endif
    if (this != &other) {
        SizeAndData = other.SizeAndData;
        other.SizeAndData = nullptr;

        Capacity = other.Capacity;
        other.Capacity = 0;
    }
    return *this;
}

template <class T>
ind MPIVector<T>::size() const {
    return extractSizeFromPointer(SizeAndData);
}

template <class T>
ind MPIVector<T>::capacity() const {
    return Capacity;
}

template <class T>
T* MPIVector<T>::data() {
    return extractDataFromPointer(SizeAndData);
}

template <class T>
const T* MPIVector<T>::data() const {
    return extractDataFromPointer(SizeAndData);
}

template <class T>
T& MPIVector<T>::operator[](ind index) {
    assert(index < size() && "Index out of bounds.");
    return data()[index];
}

template <class T>
void MPIVector<T>::push_back(const T& value) {
    ind currSize = size();
    if (currSize == Capacity) {
        resize(Capacity * CAPACITY_MULTIPLIER);
    }
    data()[currSize] = value;
    setSizeinPointer(SizeAndData, currSize + 1);
}

template <class T>
void MPIVector<T>::pop_back() {
    ind currSize = size();
    assert(currSize > 0 && "Pop on empty array.");
    // Just decrement size
    setSizeinPointer(SizeAndData, currSize - 1);
    // For now omitted: Making the array smaller if number of elements in the vector falls below
    // 1/CAPACITY_MULTIPLIER * Capacity
}

template <class T>
void MPIVector<T>::resize(ind newCapacity) {
#ifdef TEST
    std::cout << "Resizing from " << Capacity << " to " << newCapacity << "." << std::endl;
#endif
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
}

template <class T>
int MPIVector<T>::Send(int dest, int tag, MPI_Comm comm) {
    // Use MPI_BYTE to send data, capacity is not sent, will be set to size by receiver
    std::cout << size() << " " << Capacity << " " << calculatePointerSize(size()) << std::endl;
    return MPI_Send(SizeAndData, calculatePointerSize(size()), MPI_BYTE, dest, tag, comm);
}

template <class T>
int MPIVector<T>::Recv(int src, int tag, MPI_Comm comm, MPI_Status* status) {
    // Resize if we cannot hold the data to be received
    int err = MPI_Probe(src, tag, comm, status);
    int messageSize;
    MPI_Get_count(status, MPI_BYTE, &messageSize);
    ind vectorSize = (messageSize - sizeof(ind)) / sizeof(T);
    if (Capacity < vectorSize) resize(vectorSize);
    MPI_Recv(SizeAndData, messageSize, MPI_BYTE, src, tag, comm, MPI_STATUS_IGNORE);
}

template <class T>
void MPIVector<T>::printContents() const {
    std::cout << "[";
    ind vectorSize = size();
    const T* vectorData = data();
    for (ind i = 0; i < vectorSize - 1; i++) {
        std::cout << vectorData[i] << ",";
    }
    std::cout << vectorData[vectorSize - 1] << "]" << std::endl;
}

}  // namespace perc
