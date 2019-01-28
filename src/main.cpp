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
#include "localblock.h"
#include "globalblock.h"
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

    // testingMPIVectors();
    // testClusterMerges();

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

    // One master node and one node for each block
#ifndef SINGLENODE
    if (numNodes.prod() + 1 > numProcesses) {
        std::cerr << "Too few nodes. Needing " << numNodes.x << " x " << numNodes.y << " x "
                  << numNodes.z << " and one additional master" << std::endl;
        return 1;
    }
#endif

    int currProcess;
    MPI_Comm_rank(MPI_COMM_WORLD, &currProcess);

    vec3i idxNode = vec3i::fromIndexOfTotal(currProcess, numNodes);
    vec3i blockOffset = blockSize * idxNode;
    blockSize = vec3i::min(totalSize, blockSize * (idxNode + 1)) - blockOffset;

    // blockSize = {193, 194, 100};
    // totalSize = {193, 194, 100};

    // Print status
    std::cout << "Processor " << currProcess << ", index " << idxNode << ", size " << blockSize
              << std::endl;

    // TODO, put settings in command line arguments
    float hMin = 0.0;
    float hMax = 2;
    assert(hMax > hMin && "HMax needs to be larger than hMin.");
    int hSamples = 101;
    float hStep = (hMax - hMin) / (hSamples - 1);

    // Keep track of threshold h, number of components, volume largest component, volume total
    std::vector<float> h;
    std::vector<ind> numClusters;
    std::vector<float> maxVolumes;
    std::vector<float> totalVolumes;

    PerformanceTimer timer;
    timer.Reset();
    float timeElapsed;

#ifdef SINGLENODE

    if (numProcesses > 1) {
        std::cerr << "Using single node mode with " << numProcesses << " processes. " << std::endl;
        return 1;
    }

    h.reserve(hSamples);
    numClusters.reserve(hSamples);
    maxVolumes.reserve(hSamples);
    totalVolumes.reserve(hSamples);

#ifdef COMMUNICATION

    std::cout << "Running single node communication test: WhiteRed and WhiteRedGreen." << std::endl;

#ifndef NDEBUG
    // Keeps groundtruth block to compare against
    LocalBlock* groundtruth = LocalBlock::makeGroundtruth(blockSize, blockOffset, totalSize);
#endif  // !NDEBUG

    // WhiteRed
    LocalBlock* localBlockWhiteRed =
        LocalBlock::makeWhiteRedTest(blockSize, blockOffset, totalSize);
    GlobalBlock* globalBlockRed = GlobalBlock::makeWhiteRedTest(blockSize, blockOffset, totalSize);

#ifdef GREEN
    // WhiteRedGreen
    LocalBlock* localBlockWhiteRedGreen =
        LocalBlock::makeWhiteRedGreenTest(blockSize, blockOffset, totalSize);
    GlobalBlock* globalBlockRedGreen =
        GlobalBlock::makeWhiteRedGreenTest(blockSize, blockOffset, totalSize);
#endif  // GREEN

    timeElapsed = timer.ElapsedTimeAndReset();
    std::cout << "Loaded and sorted data in " << timeElapsed << " seconds." << std::endl;

    for (float currentH = hMax; currentH >= hMin - 1e-5; currentH -= hStep) {
        localBlockWhiteRed->doWatershed(currentH);
        localBlockWhiteRed->sendData();
        globalBlockRed->receiveData();
        globalBlockRed->doWatershed(currentH);
        globalBlockRed->sendData();
        localBlockWhiteRed->receiveData();
        std::cout << currentH << "/ " << hStep << "\t - "
                  << "(WhiteRed) " << globalBlockRed->numClustersCombined() << std::endl;

#ifdef GREEN
        MPI_Barrier(MPI_COMM_WORLD);

        localBlockWhiteRedGreen->doWatershed(currentH);
        localBlockWhiteRedGreen->sendData();
        globalBlockRedGreen->receiveData();
        globalBlockRedGreen->doWatershed(currentH);
        globalBlockRedGreen->sendData();
        localBlockWhiteRedGreen->receiveData();
        std::cout << currentH << "/ " << hStep << "\t - "
                  << "(WhiteRedGreen) " << globalBlockRedGreen->numClustersCombined() << std::endl;
#endif  // GREEN
        h.push_back(currentH);

#ifndef NDEBUG
        groundtruth->doWatershed(currentH);
        // WhiteRed
        if (globalBlockRed->numClustersCombined() != groundtruth->numClustersCombined() ||
            globalBlockRed->totalVolumeCombined() != groundtruth->totalVolumeCombined()) {

            if (globalBlockRed->numClustersCombined() != groundtruth->numClustersCombined())
                std::cout << "(WhiteRed) Number of clusters is "
                          << globalBlockRed->numClustersCombined() << '\\'
                          << groundtruth->numClustersCombined() << std::endl;
            if (globalBlockRed->totalVolumeCombined() != groundtruth->totalVolumeCombined())
                std::cout << "(WhiteRed) Total volume is "
                          << ind(globalBlockRed->totalVolumeCombined()) << '\\'
                          << ind(groundtruth->totalVolumeCombined()) << std::endl;

            // Compare newet additions by related cluster volume.
            auto groundStats = groundtruth->getVoluminaForAddedVertices(currentH + hStep);
            auto localStats = localBlockWhiteRed->getVoluminaForAddedVertices(currentH + hStep);
            auto globalStats = globalBlockRed->getVoluminaForAddedVertices(currentH + hStep);

            assert(localStats.size() + globalStats.size() == groundStats.size() &&
                   "(WhiteRed) Test is incorrect.");
            for (int i = 0; i < groundStats.size(); ++i) {
                auto& currentStat = groundStats[i];
                // Find in both local and global
                auto localIt =
                    std::find_if(localStats.begin(), localStats.end(),
                                 [&currentStat](const std::pair<vec3i, double>& element) {
                                     return currentStat.first == element.first;
                                 });

                auto globalIt =
                    std::find_if(globalStats.begin(), globalStats.end(),
                                 [&currentStat](const std::pair<vec3i, double>& element) {
                                     return currentStat.first == element.first;
                                 });

                // Green
                if (globalIt != globalStats.end()) {
                    if (currentStat.second != globalIt->second) {
                        std::cout << "(WhiteRed in Green part) \t" << currentStat.first
                                  << ":\tcorrect " << currentStat.second
                                  << " != " << globalIt->second << std::endl;
                    }
                }
                // Red or white
                else if (localIt != localStats.end()) {
                    vec3i dummy(-1, -1, -1);
                    ClusterID* id = localBlockWhiteRed->findClusterID(currentStat.first, dummy);
                    // Red
                    if (id->isGlobal()) {
                        if (localIt->second != 0) {
                            std::cout
                                << "(WhiteRed in Red part) \t correct 0 != " << localIt->second
                                << std::endl;
                        }
                        double volRed = globalBlockRed->getClusterVolume(*id);
                        if (currentStat.second != volRed) {
                            std::cout << "(WhiteRed in Red part) \t" << currentStat.first
                                      << ":\tcorrect " << currentStat.second << " != " << volRed
                                      << std::endl;
                        }
                    } else {
                        if (currentStat.second != localIt->second) {
                            std::cout << "(WhiteRed in White part) \t" << currentStat.first
                                      << ":\tcorrect " << currentStat.second
                                      << " != " << localIt->second << std::endl;
                        }
                    }
                }
            }

            assert(globalBlockRed->numClustersCombined() == groundtruth->numClusters() &&
                   "(WhiteRed) Groundtruth has other num clusters.");
            assert(globalBlockRed->totalVolumeCombined() == groundtruth->totalVolume() &&
                   "(WhiteRed) Groundtruth has other total volume.");
        }
#ifdef GREEN
        // WhiteRedGreen
        if (globalBlockRedGreen->numClustersCombined() != groundtruth->numClustersCombined() ||
            globalBlockRedGreen->totalVolumeCombined() != groundtruth->totalVolumeCombined()) {

            if (globalBlockRedGreen->numClustersCombined() != groundtruth->numClustersCombined())
                std::cout << "(WhiteRedGreen) Number of clusters is "
                          << globalBlockRedGreen->numClustersCombined() << '\\'
                          << groundtruth->numClustersCombined() << std::endl;
            if (globalBlockRedGreen->totalVolumeCombined() != groundtruth->totalVolumeCombined())
                std::cout << "(WhiteRedGreen) Total volume is "
                          << ind(globalBlockRedGreen->totalVolumeCombined()) << '\\'
                          << ind(groundtruth->totalVolumeCombined()) << std::endl;

            // Compare newet additions by related cluster volume.
            auto groundStats = groundtruth->getVoluminaForAddedVertices(currentH + hStep);
            auto localStats =
                localBlockWhiteRedGreen->getVoluminaForAddedVertices(currentH + hStep);
            auto globalStats = globalBlockRedGreen->getVoluminaForAddedVertices(currentH + hStep);

            assert(localStats.size() + globalStats.size() == groundStats.size() &&
                   "(WhiteRedGreen) Test is incorrect.");
            for (int i = 0; i < groundStats.size(); ++i) {
                auto& currentStat = groundStats[i];
                // Find in both local and global
                auto localIt =
                    std::find_if(localStats.begin(), localStats.end(),
                                 [&currentStat](const std::pair<vec3i, double>& element) {
                                     return currentStat.first == element.first;
                                 });

                auto globalIt =
                    std::find_if(globalStats.begin(), globalStats.end(),
                                 [&currentStat](const std::pair<vec3i, double>& element) {
                                     return currentStat.first == element.first;
                                 });

                // Green
                if (globalIt != globalStats.end()) {
                    if (currentStat.second != globalIt->second) {
                        std::cout << "(WhiteRedGreen in Green part) \t" << currentStat.first
                                  << ":\tcorrect " << currentStat.second
                                  << " != " << globalIt->second << std::endl;
                    }
                }
                // Red or white
                else if (localIt != localStats.end()) {
                    vec3i dummy(-1, -1, -1);
                    ClusterID* id =
                        localBlockWhiteRedGreen->findClusterID(currentStat.first, dummy);
                    // Red
                    if (id->isGlobal()) {
                        if (localIt->second != 0) {
                            std::cout
                                << "(WhiteRedGreen in Red part) \t correct 0 != " << localIt->second
                                << std::endl;
                        }
                        double volRed = globalBlockRedGreen->getClusterVolume(*id);
                        if (currentStat.second != volRed) {
                            std::cout << "(WhiteRedGreen in Red part) \t" << currentStat.first
                                      << ":\tcorrect " << currentStat.second << " != " << volRed
                                      << std::endl;
                        }
                    } else {
                        if (currentStat.second != localIt->second) {
                            std::cout << "(WhiteRedGreen in White part) \t" << currentStat.first
                                      << ":\tcorrect " << currentStat.second
                                      << " != " << localIt->second << std::endl;
                        }
                    }
                }

                assert(globalBlockRedGreen->numClustersCombined() == groundtruth->numClusters() &&
                       "(WhiteRedGreen) Groundtruth has other num clusters.");
                assert(globalBlockRedGreen->totalVolumeCombined() == groundtruth->totalVolume() &&
                       "(WhiteRedGreen) Groundtruth has other total volume.");
            }
        }
#endif  // GREEN
#endif  // !NDEBUG
        numClusters.push_back(globalBlockRed->numClustersCombined());
        maxVolumes.push_back(globalBlockRed->maxVolumeCombined());
        totalVolumes.push_back(globalBlockRed->totalVolumeCombined());
    }

    delete localBlockWhiteRed;
    delete globalBlockRed;
#ifdef GREEN
    delete localBlockWhiteRedGreen;
    delete globalBlockRedGreen;
#endif  // GREEN
#ifndef NDEBUG
    delete groundtruth;
#endif  // !NDEBUG

#else  // !COMMUNICATION

#ifndef TEST

    // Just localNode White (no comparison to groundtruth (This is groundtruth essentially))

    LocalBlock localBlock(blockSize, blockOffset, totalSize);

    timeElapsed = timer.ElapsedTimeAndReset();
    std::cout << "Loaded and sorted data in " << timeElapsed << " seconds." << std::endl;

    for (float currentH = hMax; currentH >= hMin - 1e-5; currentH -= hStep) {
        localBlock.doWatershed(currentH);
        h.push_back(currentH);
        std::cout << currentH << "/ " << hStep << "\t - " << localBlock.numClusters() << std::endl;
        numClusters.push_back(localBlock.numClusters());
        maxVolumes.push_back(localBlock.maxVolume());
        totalVolumes.push_back(localBlock.totalVolume());
    }

#else  // TEST

    std::cout << "Running single node tests: Local WhiteRed and Global Green." << std::endl;

#ifndef NDEBUG
    // Keeps groundtruth block to compare against
    LocalBlock* groundtruth = LocalBlock::makeGroundtruth(blockSize, blockOffset, totalSize);
#endif  // !NDEBUG

    LocalBlock* localBlockWhiteRed =
        LocalBlock::makeWhiteRedTest(blockSize, blockOffset, totalSize);
#ifdef GREEN
    GlobalBlock* globalBlockGreen = GlobalBlock::makeGreenTest(blockSize, blockOffset, totalSize);
#endif  // GREEN

    timeElapsed = timer.ElapsedTimeAndReset();
    std::cout << "Loaded and sorted data in " << timeElapsed << " seconds." << std::endl;

    for (float currentH = hMax; currentH >= hMin - 1e-5; currentH -= hStep) {
        localBlockWhiteRed->doWatershed(currentH);
        // Send and receive models the receiving and sending process, but actually just cleans
        // up
        localBlockWhiteRed->sendData();
        localBlockWhiteRed->receiveData();
        std::cout << currentH << "/ " << hStep << "\t - "
                  << "(WhiteRed) " << localBlockWhiteRed->numClustersCombined() << std::endl;

#ifdef GREEN
        globalBlockGreen->doWatershed(currentH);
        // No send/receive
        std::cout << currentH << "/ " << hStep << "\t - "
                  << "(Green) " << globalBlockGreen->numClustersCombined() << std::endl;
#endif  //  GREEN

        h.push_back(currentH);
#ifndef NDEBUG
        groundtruth->doWatershed(currentH);
        // WhiteRed
        if (localBlockWhiteRed->numClustersCombined() != groundtruth->numClustersCombined() ||
            localBlockWhiteRed->totalVolumeCombined() != groundtruth->totalVolumeCombined()) {

            if (localBlockWhiteRed->numClustersCombined() != groundtruth->numClustersCombined())
                std::cout << "(WhiteRed) Number of clusters is "
                          << localBlockWhiteRed->numClustersCombined() << '\\'
                          << groundtruth->numClustersCombined() << std::endl;
            if (localBlockWhiteRed->totalVolumeCombined() != groundtruth->totalVolumeCombined())
                std::cout << "(WhiteRed) Total volume is "
                          << ind(localBlockWhiteRed->totalVolumeCombined()) << '\\'
                          << ind(groundtruth->totalVolumeCombined()) << std::endl;

            // Compare newet additions by related cluster volume.
            auto groundStats = groundtruth->getVoluminaForAddedVertices(currentH + hStep);
            auto whiteRedStats = localBlockWhiteRed->getVoluminaForAddedVertices(currentH + hStep);

            assert(groundStats.size() == whiteRedStats.size() && "(WhiteRed) Test is incorrect.");
            for (int i = 0; i < groundStats.size(); ++i) {
                assert(groundStats[i].first == whiteRedStats[i].first &&
                       "(WhiteRed) Test is incorrect.");
                if (groundStats[i].second != whiteRedStats[i].second) {
                    std::cout << "(WhiteRed) \t" << groundStats[i].first << ":\tcorrect "
                              << groundStats[i].second << " != " << whiteRedStats[i].second
                              << std::endl;
                }
            }

            assert(localBlockWhiteRed->numClustersCombined() == groundtruth->numClusters() &&
                   "(WhiteRed) Groundtruth has other num clusters.");
            assert(localBlockWhiteRed->totalVolumeCombined() == groundtruth->totalVolume() &&
                   "(WhiteRed) Groundtruth has other total volume.");
        }
#ifdef GREEN
        // Green
        if (globalBlockGreen->numClustersCombined() != groundtruth->numClustersCombined() ||
            globalBlockGreen->totalVolumeCombined() != groundtruth->totalVolumeCombined()) {

            if (globalBlockGreen->numClustersCombined() != groundtruth->numClustersCombined())
                std::cout << "(Green) Number of clusters is "
                          << globalBlockGreen->numClustersCombined() << '\\'
                          << groundtruth->numClustersCombined() << std::endl;
            if (globalBlockGreen->totalVolumeCombined() != groundtruth->totalVolumeCombined())
                std::cout << "(Green) Total volume is "
                          << ind(globalBlockGreen->totalVolumeCombined()) << '\\'
                          << ind(groundtruth->totalVolumeCombined()) << std::endl;

            // Compare newest additions by related cluster volume.
            auto groundStats = groundtruth->getVoluminaForAddedVertices(currentH + hStep);
            auto greenStats = globalBlockGreen->getVoluminaForAddedVertices(currentH + hStep);

            assert(groundStats.size() == greenStats.size() && "(Green) Test is incorrect.");
            for (int i = 0; i < groundStats.size(); ++i) {
                assert(groundStats[i].first == greenStats[i].first && "(Green) Test is incorrect.");
                if (groundStats[i].second != greenStats[i].second) {
                    std::cout << "(Green) \t" << groundStats[i].first << ":\tcorrect "
                              << groundStats[i].second << " != " << greenStats[i].second
                              << std::endl;
                }
            }

            assert(globalBlockGreen->numClustersCombined() == groundtruth->numClusters() &&
                   "(Green) Groundtruth has other num clusters.");
            assert(globalBlockGreen->totalVolumeCombined() == groundtruth->totalVolume() &&
                   "(Green) Groundtruth has other total volume.");
        }
#endif  // GREEN
#endif  // !NDEBUG
        numClusters.push_back(localBlockWhiteRed->numClustersCombined());
        maxVolumes.push_back(localBlockWhiteRed->maxVolumeCombined());
        totalVolumes.push_back(localBlockWhiteRed->totalVolumeCombined());
    }

    delete localBlockWhiteRed;
#ifdef GREEN
    delete globalBlockGreen;
#endif  // GREEN
#ifndef NDEBUG
    delete groundtruth;
#endif  // !NDEBUG

#endif  // TEST
#endif  // COMMUNICATION

// Not single node use -> The real parallel distributed thing
#else  // !SINGLENODE

    if (numProcesses == 1) {
        std::cerr << "Using multiple node mode with a single processes. " << std::endl;
        return 1;
    }

#ifndef NDEBUG
    // Keeps groundtruth block to compare against
    LocalBlock* groundtruth = LocalBlock::makeGroundtruth(blockSize, blockOffset, totalSize);
#endif  // NDEBUG

    // Master process
    if (currProcess == 0) {

        std::cout << "Running in parallel." << std::endl;

        // Prepare output
        h.reserve(hSamples);
        numClusters.reserve(hSamples);
        maxVolumes.reserve(hSamples);
        totalVolumes.reserve(hSamples);

        GlobalBlock globalBlock(blockSize, blockOffset, totalSize, numNodes);

        timeElapsed = timer.ElapsedTimeAndReset();
        std::cout << "Processor " << currProcess << ": Loaded and sorted data in " << timeElapsed
                  << " seconds." << std::endl;

        for (float currentH = hMax; currentH >= hMin - 1e-5; currentH -= hStep) {
            globalBlock.receiveData();
            globalBlock.doWatershed(currentH);
            globalBlock.sendData();
            std::cout << currentH << "/ " << hStep << "\t - " << globalBlock.numClustersCombined()
                      << std::endl;
            h.push_back(currentH);
#ifndef NDEBUG
            groundtruth->doWatershed(currentH);
            if (globalBlock.numClustersCombined() != groundtruth->numClustersCombined() ||
                globalBlock.totalVolumeCombined() != groundtruth->totalVolumeCombined()) {

                if (globalBlock.numClustersCombined() != groundtruth->numClustersCombined())
                    std::cout << "Number of clusters is " << globalBlock.numClustersCombined()
                              << '\\' << groundtruth->numClustersCombined() << std::endl;
                if (globalBlock.totalVolumeCombined() != groundtruth->totalVolumeCombined())
                    std::cout << "Total volume is " << ind(globalBlock.totalVolumeCombined())
                              << '\\' << ind(groundtruth->totalVolumeCombined()) << std::endl;

                // Let the other processes know something is wrong
                bool correct;
#ifndef COLLECTIVES
                MPI_Request* requests = new MPI_Request[numNodes.prod()];
                for (ind p = 0; p < numNodes.prod(); p++) {
                    ind processindex = p + 1;
                    MPI_Isend(&correct, 1, MPI_CXX_BOOL, processindex, MPICommunication::ERRORFLAG,
                              MPI_COMM_WORLD, &requests[p]);
                }
#else
                // TODO: Broadcast
#endif

                // We can only check if newely processed green indices have the correct volume
                // (White and red data is held by other nodes)
                auto groundStats = groundtruth->getVoluminaForAddedVertices(currentH + hStep);
                auto greenStats = globalBlock.getVoluminaForAddedVertices(currentH + hStep);

                for (int i = 0; i < greenStats.size(); ++i) {
                    auto& currentStat = greenStats[i];
                    // Find corresponding value in ground truth
                    auto truthIt =
                        std::find_if(groundStats.begin(), groundStats.end(),
                                     [&currentStat](const std::pair<vec3i, double>& element) {
                                         return currentStat.first == element.first;
                                     });
                    assert(truthIt != groundStats.end() && "Test is incorrect.");
                    if (currentStat.second != truthIt->second) {
                        std::cout << "Green part " << currProcess << "\t" << currentStat.first
                                  << ":\tcorrect " << truthIt->second
                                  << " != " << currentStat.second << std::endl;
                    }
                }

                std::vector<vec3i> redIndices;
                MPI_Status* statuses = new MPI_Status[numNodes.prod()];

                // Receive red indices to check from each node
                for (ind p = 0; p < numNodes.prod(); p++) {
                    ind processindex = p + 1;
                    // redIndices.clear();
                    MPICommunication::RecvVectorUknownSize(redIndices, processindex,
                                                           MPICommunication::REDINDICES,
                                                           MPI_COMM_WORLD, &statuses[p]);

                    for (vec3i& redIndex : redIndices) {
                        auto truthIt =
                            std::find_if(groundStats.begin(), groundStats.end(),
                                         [&redIndex](const std::pair<vec3i, double>& element) {
                                             return redIndex == element.first;
                                         });

                        vec3i dummy(-1, -1, -1);
                        ClusterID* id = globalBlock.findClusterID(redIndex, dummy);
                        double volRed = globalBlock.getClusterVolume(*id);

                        if (volRed != truthIt->second) {
                            std::cout << "Red part " << currProcess << "\t" << redIndex
                                      << ":\tcorrect " << truthIt->second << " != " << volRed
                                      << std::endl;
                        }
                    }
                }

                assert(globalBlock.numClustersCombined() == groundtruth->numClustersCombined() &&
                       "Groundtruth has other num clusters.");
                assert(globalBlock.totalVolumeCombined() == groundtruth->totalVolumeCombined() &&
                       "Groundtruth has other total volume.");
            }

#endif  // NDEBUG
            numClusters.push_back(globalBlock.numClustersCombined());
            maxVolumes.push_back(globalBlock.maxVolumeCombined());
            totalVolumes.push_back(globalBlock.totalVolumeCombined());
        }
    }

    // All other processes: currProcess != 0
    else {
        LocalBlock localBlock(blockSize, blockOffset, totalSize);
        timeElapsed = timer.ElapsedTimeAndReset();
        std::cout << "Processor " << currProcess << ": Loaded and sorted data in " << timeElapsed
                  << " seconds." << std::endl;

        for (float currentH = hMax; currentH >= hMin - 1e-5; currentH -= hStep) {
            localBlock.doWatershed(currentH);
            localBlock.sendData();
            std::cout << currentH << "/ " << hStep << "\t - " << localBlock.numClusters()
                      << std::endl;
#ifndef NDEBUG
            groundtruth->doWatershed(currentH);

            bool correct;
            MPI_Recv(&correct, 1, MPI_CXX_BOOL, 0, MPICommunication::ERRORFLAG, MPI_COMM_WORLD,
                     MPI_STATUS_IGNORE);

            if (!correct) {
                std::vector<vec3i> redIndices;

                // Check white indices and collect red indices so that global node can check it
                auto groundStats = groundtruth->getVoluminaForAddedVertices(currentH + hStep);
                auto WhiteAndRedStats = localBlock.getVoluminaForAddedVertices(currentH + hStep);

                for (int i = 0; i < WhiteAndRedStats.size(); ++i) {
                    auto& currentStat = WhiteAndRedStats[i];
                    // Find in both local and global
                    auto truthIt =
                        std::find_if(groundStats.begin(), groundStats.end(),
                                     [&currentStat](const std::pair<vec3i, double>& element) {
                                         return currentStat.first == element.first;
                                     });
                    assert(truthIt != groundStats.end() && "Test is incorrect.");
                    if (currentStat.second != truthIt->second) {
                        if (currentStat.second == 0) {
                            // Volume has been nulled, we assume it to be part of red
                            redIndices.push_back(currentStat.first);
                        } else {
                            // Not nulled (so white?), but incorrect.
                            std::cout << "White part in process " << currProcess << "\t"
                                      << currentStat.first << ":\tcorrect " << truthIt->second
                                      << " != " << currentStat.second << std::endl;
                        }
                    }
                }

                MPI_Request request;
                int err = MPICommunication::IsendVector(redIndices, 0, MPICommunication::ERRORFLAG,
                                                        MPI_COMM_WORLD, &request);
                MPI_Wait(&request, MPI_STATUS_IGNORE);
            }

#endif  // NDEBUG

            // Prepate local block for next step
            localBlock.receiveData();
        }
    }

#ifndef NDEBUG
    delete groundtruth;
#endif

#endif  // SINGLENODE

    timeElapsed = timer.ElapsedTimeAndReset();
    std::cout << "Processor " << currProcess << ": Watershedding took " << timeElapsed
              << " seconds." << std::endl;

    // Only master process writes out statistics
    if (currProcess == 0) {

        const int maxClusters = *(std::max_element(numClusters.cbegin(), numClusters.cend()));

        std::ofstream percFile;
        std::string fileName = "percolation.csv";

        percFile.open(fileName, std::ios::out);

        if (percFile.is_open()) {
            percFile
                << "H; Number of connected components; Maximum number of connected components; "
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
    }

    // Finalize the MPI environment. No more MPI calls can be made after this
    MPI_Finalize();
}