// Basic MPI setup adapted from www.mpitutorial.com

#include <mpi.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath>
#include "vec.h"
#include "datablock.h"
#include "mpicommuncation.h"
#include "performancetimer.h"
#include "whiteblock.h"
#include "greenblock.h"
#include "clusterlistrecording.h"

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
        std::vector<double> doubleVector;
        doubleVector.push_back(5.0);
        doubleVector.push_back(5.0);
        doubleVector.push_back(6.0);
        doubleVector.push_back(4.5);

        std::vector<int> intVector(7, 9);
        intVector.push_back(4);

        // Send a int vector
        if (numProcesses > 2) {
            MPI_Request* requests = new MPI_Request[2];
            MPICommunication::IsendVector(intVector, 1, 0, MPI_COMM_WORLD, requests);
            MPICommunication::IsendVector(doubleVector, 2, 1, MPI_COMM_WORLD, &requests[1]);

            // Wait for all messages to be finished
            MPI_Waitall(1, requests, MPI_STATUSES_IGNORE);
        }

    } else if (currProcess == 1) {
        std::vector<int> received;
        MPI_Status status;
        MPICommunication::RecvVectorUknownSize(received, 0, 0, MPI_COMM_WORLD, &status);
    } else if (currProcess == 2) {
        std::vector<double> received;
        MPI_Status status;
        MPICommunication::RecvVectorUknownSize(received, 0, 1, MPI_COMM_WORLD, &status);
    }

    MPI_Finalize();
}

void percolationSequential() { std::cout << "Running percolation in a single node" << std::endl; }

void percolation() {}

void percolationDistributed() {}

void testClusterMerges() {
    ind numMerges = 15;
    ind numClusters = 20;

    std::vector<std::vector<ClusterMerge>> merges(2);
    for (auto& g : merges)
        for (ind merge = 0; merge < numMerges / 2; ++merge) {
            ind f = rand() % numClusters;
            ind o = rand() % numClusters;
            g.push_back(ClusterMerge(f, o));
            std::cout << f << " -> " << o << '\n';
        }

    // Graph from the edge list.
    auto graphs = ClusterMerge::mergeClustersFromLists(merges);
    std::cout << "\nGraphs:\n";
    for (auto& g : graphs) {
        for (ind n : g) std::cout << ID(n).baseID() << " - ";
        std::cout << '\n';
    }

    // Into one single vector for transfer.
    std::cout << "\nAs one vector\n";
    auto graphAsOne = ClusterMerge::mergeClusterAsList(graphs);
    for (ind n : graphAsOne) std::cout << ID(n).baseID() << ", ";
}

// Input args:
// - path:        folder with data
// - rms:         name of rms file
// - xT, yT, zT:  total size
// - xB, yB, zB:  block size
int main(int argc, char** argv) {

#ifdef TEST

    testClusterMerges();

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

    // blockSize = {193, 194, 1000};

    // Print status
    std::cout << "Processor " << currProcess << ", index " << idxNode << ", size " << blockSize
              << std::endl;

    PerformanceTimer timer;
    timer.Reset();
    float timeElapsed;

    WhiteBlock localBlock(blockSize, blockOffset, totalSize);
    GreenBlock globalBlock(blockSize, blockOffset, blockSize, numNodes);

    // To compare against
    WhiteBlock* groundtruth = WhiteBlock::makeGroundtruth(blockSize, blockOffset, totalSize);

    timeElapsed = timer.ElapsedTimeAndReset();
    std::cout << "Loaded and sorted data in " << timeElapsed << " seconds." << std::endl;

    // TODO, put settings in command line arguments
    float hMin = 0.0;
    float hMax = 2;
    assert(hMax > hMin && "HMax needs to be larger than hMin.");
    int hSamples = 100;
    float hStep = (hMax - hMin) / (hSamples - 1);

    // Keep track of threshold h, number of components, volume largest component, volume total
    std::vector<float> h;
    h.reserve(hSamples);
    std::vector<ind> numClusters;
    numClusters.reserve(hSamples);
    std::vector<float> maxVolumes;
    maxVolumes.reserve(hSamples);
    std::vector<float> totalVolumes;
    totalVolumes.reserve(hSamples);

    for (float currentH = hMax; currentH >= hMin - 1e-5; currentH -= hStep) {

        // TODO: Make this dependent on process if (currProcess != numProcesses - 1) {}

        // Local Block does watershed and sends result over to the global block
        localBlock.doWatershed(currentH);
        localBlock.sendData();
        // Global Block
        globalBlock.receiveData();
        globalBlock.doWatershed(currentH);
        globalBlock.sendData();

        // Prepate local block for next step
        localBlock.receiveData();

        groundtruth->doWatershed(currentH);

        std::cout << currentH << "/ " << hStep << "\t - " << globalBlock.numClusters() << std::endl;
        h.push_back(currentH);
#ifndef NDEBUG
        if (globalBlock.totalNumClusters() != groundtruth->totalNumClusters() ||
            globalBlock.totalTotalVolume() != groundtruth->totalTotalVolume()) {

            if (globalBlock.totalNumClusters() != groundtruth->totalNumClusters())
                std::cout << "Number of clusters is " << globalBlock.totalNumClusters() << '\\'
                          << groundtruth->totalNumClusters() << std::endl;
            if (globalBlock.totalTotalVolume() != groundtruth->totalTotalVolume())
                std::cout << "Total volume is " << ind(globalBlock.totalTotalVolume()) << '\\'
                          << ind(groundtruth->totalTotalVolume()) << std::endl;

            // Compare newet additions by related cluster volume.
            auto groundStats = groundtruth->getVoluminaForAddedVertices(currentH + hStep);
            for (int i = 0; i < groundStats.size(); ++i) {
                vec3i idx = groundStats[i].first;
                vec3i dummy(-1, -1, -1);
                double combinedVolume = 0;
                // TODO: Hard coded uglyness
                // White
                if (idx.liesWithin(vec3i(0, 0, 10), vec3i(193, 194, 980))) {
                    ClusterID cluster = *localBlock.findClusterID(idx, dummy);
                    combinedVolume = localBlock.getClusterVolume(cluster);
                } else {
                    ClusterID cluster = *localBlock.findClusterID(idx, dummy);
                    combinedVolume += localBlock.getClusterVolume(cluster);
                    cluster = *globalBlock.findClusterID(idx, dummy);
                    combinedVolume += globalBlock.getClusterVolume(cluster);
                }
                if (groundStats[i].second != combinedVolume) {
                    std::cout << '\t' << groundStats[i].first << ":\tcorrect "
                              << groundStats[i].second << " != " << combinedVolume << std::endl;
                }
            }

            /*auto whiteStats = localBlock.getVoluminaForAddedVertices(currentH + hStep);
            auto greenStats = globalBlock.getVoluminaForAddedVertices(currentH + hStep);

            // Combine green and white status into one
            std::vector<std::pair<vec3i, double>> combinedStats;
            auto itWhite = whiteStats.begin();
            auto itGreen = greenStats.begin();
            while (itWhite != whiteStats.end() && itGreen != greenStats.end()) {
                auto whiteIdx = itWhite->first.toIndexOfTotal(localBlock.TotalSize);
                auto greenIdx = itGreen->first.toIndexOfTotal(globalBlock.TotalSize);
                if (whiteIdx < greenIdx) {
                    combinedStats.push_back(*itWhite);
                    itWhite++;
                } else if (whiteIdx > greenIdx) {
                    combinedStats.push_back(*itGreen);
                    itGreen++;
                } else {
                    combinedStats.emplace_back(itGreen->first, itGreen->second + itWhite->second);
                    itWhite++;
                    itGreen++;
                }
            }
            while (itWhite != whiteStats.end()) {
                combinedStats.push_back(*itWhite);
                itWhite++;
            }
            while (itGreen != greenStats.end()) {
                combinedStats.push_back(*itGreen);
                itGreen++;
            }*/

            /*assert(groundStats.size() == combinedStats.size() && "Test is incorrect.");
            for (int i = 0; i < groundStats.size(); ++i) {
                assert(groundStats[i].first == combinedStats[i].first && "Test is incorrect.");
                if (groundStats[i].second != combinedStats[i].second) {
                    std::cout << '\t' << groundStats[i].first << ":\tcorrect "
                              << groundStats[i].second << " != " << combinedStats[i].second
                              << std::endl;
                }
            }*/

            assert(globalBlock.totalNumClusters() == groundtruth->totalNumClusters() &&
                   "Groundtruth has other num clusters.");
            assert(globalBlock.totalTotalVolume() == groundtruth->totalTotalVolume() &&
                   "Groundtruth has other total volume.");
        }
#endif
        numClusters.push_back(globalBlock.totalNumClusters());
        maxVolumes.push_back(globalBlock.totalMaxVolume());
        totalVolumes.push_back(globalBlock.totalTotalVolume());
    }

    delete groundtruth;

    timeElapsed = timer.ElapsedTimeAndReset();
    std::cout << "Watershedding took " << timeElapsed << " seconds." << std::endl;

    const int maxClusters = *(std::max_element(numClusters.cbegin(), numClusters.cend()));

    std::ofstream percFile;
    std::string fileName = "percolation.csv";

    percFile.open(fileName, std::ios::out);

    if (percFile.is_open()) {
        percFile << "H; Number of connected components; Maximum number of connected components; "
                    "Number of connected components / Maximum number of connected components;  "
                    "Largest Volume ; Total Volume; Largest Volume / Total Volume;"
                 << std::endl;
        for (int line = 0; line < h.size(); line++) {
            percFile << h[line] << ";" << float(numClusters[line]) << ";" << float(maxClusters)
                     << ";" << float(numClusters[line]) / float(maxClusters) << ";"
                     << maxVolumes[line] << ";" << totalVolumes[line] << ";"
                     << maxVolumes[line] / totalVolumes[line] << ";" << std::endl;
        }
    }

    timeElapsed = timer.ElapsedTime();
    std::cout << "Writing statistics took " << timeElapsed << " seconds." << std::endl;

    // Finalize the MPI environment. No more MPI calls can be made after this
    MPI_Finalize();

#endif
}