#pragma once

#include "vec.h"
#include <mpi.h>
#include <vector>
#include <string.h>
#include <cassert>
#include <cmath>

namespace perc {

#define MPI_IND MPI_LONG_LONG

class MPICommunication {
public:
    const static ind RANK_SHIFT = 7;

    enum Tags {
        VECTORSIZES,
        NUMCLUSTERS,
        NUMNEWCLUSTERS,
        MERGES,
        TOTALVOLUME,
        MAXVOLUME,
        VOLUMES,
        PLOGS,
        STARTOFPLOG,
        REDPOINTERS,
        GREENPOINTERS = 1 << (RANK_SHIFT - 1),  // Up to 27 green blocks -> ids take up 5 bits
        REDINDICES,
        ERRORFLAG,
        LOADTIME,
        COMMUNCATIONTIME,
        WATERSHEDTIME,
        MEMESTIMATE
    };

    template <typename T>
    static int SendVector(const std::vector<T>& vector, int dest, int tag, MPI_Comm comm);

    template <typename T>
    static int IsendVector(const std::vector<T>& vector, int dest, int tag, MPI_Comm comm,
                           MPI_Request* request);

    template <typename T>
    static int RecvVectorUknownSize(std::vector<T>& vector, int src, int tag, MPI_Comm comm,
                                    MPI_Status* status);

    static void handleError(int errorCode);

    static int computeProcess(const vec3i& idx, const vec3i& blockSize, const vec3i& numNodes);
};

template <typename T>
int MPICommunication::SendVector(const std::vector<T>& vector, int dest, int tag, MPI_Comm comm) {
    return MPI_Send(vector.data(), vector.size() * sizeof(T), MPI_BYTE, dest, tag, comm);
}

template <typename T>
int MPICommunication::IsendVector(const std::vector<T>& vector, int dest, int tag, MPI_Comm comm,
                                  MPI_Request* request) {
    // Sends an empty message for size equal to zero, bugger does not matter send
    return MPI_Isend(vector.data(), vector.size() * sizeof(T), MPI_BYTE, dest, tag, comm, request);
}

template <typename T>
int MPICommunication::RecvVectorUknownSize(std::vector<T>& vector, int src, int tag, MPI_Comm comm,
                                           MPI_Status* status) {

    int err = MPI_Probe(src, tag, comm, status);
    if (err != MPI_SUCCESS) return err;
    int messageSize;
    MPI_Get_count(status, MPI_BYTE, &messageSize);
    ind vectorSize = messageSize / sizeof(T);
    // Resize if we cannot hold the data to be received
    if (vector.size() != vectorSize) vector.resize(vectorSize);
    return MPI_Recv(vector.data(), messageSize, MPI_BYTE, src, tag, comm, MPI_STATUS_IGNORE);
}

inline void MPICommunication::handleError(int errorCode) {
    int currProcess;
    MPI_Comm_rank(MPI_COMM_WORLD, &currProcess);
    if (errorCode != MPI_SUCCESS) {
        char errorString[MPI_MAX_ERROR_STRING];
        int errorStringLength, errorClass;
        MPI_Error_class(errorCode, &errorClass);
        MPI_Error_string(errorClass, errorString, &errorStringLength);
        std::cerr << "Error in rank: " << currProcess
                  << "(class: " << std::string(errorString, errorString + errorStringLength);
        MPI_Error_string(errorCode, errorString, &errorStringLength);
        std::cerr << ", description: " << std::string(errorString, errorString + errorStringLength)
                  << ")" << std::endl;
    }
}

// This works only for the original blockSize settings: Do not use this from the localNode,
inline int MPICommunication::computeProcess(const vec3i& idx, const vec3i& blockSize,
                                            const vec3i& numNodes) {
    vec3i idxNode;
    for (int dim = 0; dim < 3; dim++) {
        idxNode[dim] = static_cast<int>(floor(static_cast<double>(idx[dim]) / blockSize[dim]));
    }
    return idxNode.toIndexOfTotal(numNodes);
}

}  // namespace perc
