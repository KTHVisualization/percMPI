// Basic MPI setup adapted from www.mpitutorial.com

#include <mpi.h>
#include <stdio.h>
#include <iostream>
#include <cassert>
#include <cmath>
#include "vec.h"
#include "datablock.h"
#include "mpivector.h"
#include "performancetimer.h"

using namespace perc;

void testingMPIVectors() {
    // MPI init
    MPI_Init(NULL, NULL);

    int numProcesses;
    MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);

    int currProcess;
    MPI_Comm_rank(MPI_COMM_WORLD, &currProcess);

    // Processor 0 is just testing some basic things before sending anything
    if (currProcess == 0) {
        // Initially empty vector, push some elements
        MPIVector<double> doubleVector;
        std::cout << doubleVector;
        doubleVector.push_back(5.0);
        doubleVector.push_back(5.0);
        doubleVector.push_back(6.0);
        std::cout << doubleVector[2] << std::endl;
        doubleVector.push_back(4.5);
        std::cout << doubleVector[3] << std::endl;
        std::cout << doubleVector;
        doubleVector.printContents();
        doubleVector.pop_back();
        std::cout << doubleVector;

        MPIVector<int> intVector(7, 9);
        std::cout << intVector;
        intVector.printContents();
        intVector.push_back(4);
        intVector.pop_back();
        intVector.pop_back();
        intVector.printContents();
        std::cout << intVector;

        // Move constructor test
        MPIVector<double> a;
        MPIVector<double> b(std::move(a));

        // Move assignment
        MPIVector<double> c;
        c = std::move(b);

        // Copy constructor
        MPIVector<int> d;
        MPIVector<int> e(d);

        // Copy assignment
        MPIVector<int> f;
        f = d;

        // Send a int vector
        if (numProcesses > 1) {
            intVector.Send(1, 0, MPI_COMM_WORLD);
        }

    } else if (currProcess == 1) {
        MPIVector<int> received;
        MPI_Status status;
        received.Recv(sizeof(ind) + 6 * sizeof(int), 0, 0, MPI_COMM_WORLD, &status);
        std::cout << received;
        received.printContents();
    }

    MPI_Finalize();
}

void percolationSequential() { std::cout << "Running percolation in a single node" << std::endl; }

void percolation() {}

void percolationDistributed() {}

// Input args:
// - path:        folder with data
// - rms:         name of rms file
// - xT, yT, zT:  total size
// - xB, yB, zB:  block size
int main(int argc, char** argv) {

#ifdef TEST

    testingMPIVectors();

#else

    // Directory path and sizes from args
    if (argc < 9) {
        std::cerr << "Not enough arguments.\n";
        return 1;
    }

    // First argument is the executable name.
    // for (int a = 0; a < argc; ++a) std::cout << a << ": " << argv[a] << std::endl;
    char* baseFolder = argv[1];
    char* rmsFilename = argv[2];
    vec3i totalSize(atoi(argv[3]), atoi(argv[4]), atoi(argv[5]));
    vec3i blockSize(atoi(argv[6]), atoi(argv[7]), atoi(argv[8]));
    vec3i numNodes;

    for (int n = 0; n < 3; ++n) {
        double s = static_cast<double>(totalSize[n]) / blockSize[n];
        numNodes[n] = static_cast<int>(ceil(s));
    }

    // MPI init
    MPI_Init(NULL, NULL);

    int numProcesses;
    MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);

    if (numNodes.prod() > numProcesses) {
        std::cerr << "Too few nodes. Needing " << numNodes.x << " x " << numNodes.y << " x "
                  << numNodes.z << '\n';
        return 1;
    }

    int currProcess;
    MPI_Comm_rank(MPI_COMM_WORLD, &currProcess);

    vec3i idxNode = vec3i::fromIndexOfTotal(currProcess, numNodes);
    vec3i blockOffset = blockSize * idxNode;
    blockSize = vec3i::min(totalSize, blockSize * (idxNode + 1)) - blockOffset;

    DataBlock data(blockSize, blockOffset, totalSize);
    bool result = data.loadData(1, baseFolder, rmsFilename);
    // Print status
    std::cout << "Processor " << currProcess << ", index " << idxNode << ", size " << blockSize
              << (result ? " - Worked!\n" : " - Nope :(\n") << std::endl;

    // Finalize the MPI environment. No more MPI calls can be made after this
    MPI_Finalize();

#endif
}