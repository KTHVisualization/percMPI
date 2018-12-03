#pragma once

#include <mpi.h>
#include <vector>
#include <string.h>
#include <cassert>

const unsigned int DEFAULT_CAPACITY = 16;
const float CAPACITY_MULTIPLIER = 2.0;

using ind = signed int;

namespace perc {

class MPICommunication {
public:
    template <typename T>
    static int SendVector(std::vector<T>& vector, int dest, int tag, MPI_Comm comm);

    template <typename T>
    static int IsendVector(std::vector<T>& vector, int dest, int tag, MPI_Comm comm,
                           MPI_Request* request);

    template <typename T>
    static int RecvVectorUknownSize(std::vector<T>& vector, int src, int tag, MPI_Comm comm,
                                    MPI_Status* status);
};

template <typename T>
int MPICommunication::SendVector(std::vector<T>& vector, int dest, int tag, MPI_Comm comm) {
    return MPI_Send(&vector[0], vector.size() * sizeof(T), MPI_BYTE, dest, tag, comm);
}

template <typename T>
int MPICommunication::IsendVector(std::vector<T>& vector, int dest, int tag, MPI_Comm comm,
                                  MPI_Request* request) {
    return MPI_Isend(&vector[0], vector.size() * sizeof(T), MPI_BYTE, dest, tag, comm, request);
}

template <typename T>
int MPICommunication::RecvVectorUknownSize(std::vector<T>& vector, int src, int tag, MPI_Comm comm,
                                           MPI_Status* status) {
    // Resize if we cannot hold the data to be received
    int err = MPI_Probe(src, tag, comm, status);
    if (err != MPI_SUCCESS) return err;
    int messageSize;
    MPI_Get_count(status, MPI_BYTE, &messageSize);
    ind vectorSize = messageSize / sizeof(T);
    if (vector.capacity() < vectorSize) vector.resize(vectorSize);
    return MPI_Recv(&vector[0], messageSize, MPI_BYTE, src, tag, comm, MPI_STATUS_IGNORE);
}

}  // namespace perc
