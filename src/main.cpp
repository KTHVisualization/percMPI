// Basic MPI setup adapted from www.mpitutorial.com

#include <mpi.h>
#include <stdio.h>
#include <iostream>
#include <cassert>
#include <cmath>
#include "vec.h"
#include "datablock.h"

using namespace perc;

// Input args:
// - path:        folder with data
// - rms:         name of rms file
// - xT, yT, zT:  total size
// - xB, yB, zB:  block size
int main(int argc, char** argv) {

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
}