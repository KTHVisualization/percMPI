#pragma once

#include <mpi.h>

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
    void resize(ind size);
    ind capacity() const;
    void reserve(ind size);

    MPIVector& operator=(MPIVector&& other);
    MPIVector& operator[](ind index);

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

}  // namespace perc