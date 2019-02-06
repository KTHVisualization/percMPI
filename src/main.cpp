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

enum Mode {
    REAL = 0,
    TEST_LOCALONLY_WHITERED = 1,
    TEST_GLOBALONLY_GREEN = 2,
    TEST_WHITERED = 3,
    TEST_WHITEREDGREEN = 4,
    TEST_COMMUNICATION = 10,
    TEST_MERGING = 11,
};

void testingMPIVectors() {
    // MPI init
    MPI_Init(NULL, NULL);

    int numProcesses;
    MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);

    int currProcess;
    MPI_Comm_rank(MPI_COMM_WORLD, &currProcess);

    // Rank  0 is just testing some basic things before sending anything
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

void watershedSequential(vec3i blockSize, vec3i blockOffset, vec3i totalSize, float hMin,
                         float hMax, float hStep, std::vector<float>& h,
                         std::vector<ind>& numClusters, std::vector<float>& maxVolumes,
                         std::vector<float>& totalVolumes, PerformanceTimer& timer) {
    // Just localNode White: Groundtruth for all other cases
    LocalBlock localBlock(blockSize, blockOffset, totalSize);

    float timeElapsed = timer.ElapsedTimeAndReset();
    std::cout << "Loaded and sorted data in " << timeElapsed << " seconds." << std::endl;

    for (float currentH = hMax; currentH >= hMin - 1e-5; currentH -= hStep) {
        localBlock.doWatershed(currentH);
        h.push_back(currentH);
        std::cout << currentH << "/ " << hStep << "\t - " << localBlock.numClusters() << std::endl;
        numClusters.push_back(localBlock.numClusters());
        maxVolumes.push_back(localBlock.maxVolume());
        totalVolumes.push_back(localBlock.totalVolume());
    }
}

void watershedParallelSingleRank(vec3i numNodes, vec3i blockSize, vec3i blockOffset,
                                 vec3i totalSize, float hMin, float hMax, float hStep,
                                 std::vector<float>& h, std::vector<ind>& numClusters,
                                 std::vector<float>& maxVolumes, std::vector<float>& totalVolumes,
                                 PerformanceTimer& timer) {
#ifndef NDEBUG
    // Keeps groundtruth block to compare against
    LocalBlock* groundtruth = LocalBlock::makeGroundtruth(totalSize, vec3i(0), totalSize);
    // The first slice of blocks for debugging.
    std::vector<char> debugSlice(totalSize.x * totalSize.y, ' ');
#endif  // NDEBUG
    std::vector<LocalBlock> localBlocks;

    vec3i idxNode;
    for (ind nodeIdx = 0; nodeIdx < numNodes.prod(); ++nodeIdx) {
        // std::cout << "Node " << nodeIdx << std::endl;
        idxNode = vec3i::fromIndexOfTotal(nodeIdx, numNodes);
        blockOffset = blockSize * idxNode;
        blockSize = vec3i::min(totalSize, blockSize * (idxNode + 1)) - blockOffset;

        localBlocks.emplace_back(blockSize, blockOffset, totalSize);
#ifndef NDEBUG
        localBlocks.back().outputFrontBlocks(debugSlice, 0, 0);
#endif
    }
    // Master block
    GlobalBlock globalBlock(blockSize, totalSize, numNodes);

#ifndef NDEBUG
    for (ind y = 0; y < totalSize.y; ++y) {
        for (ind x = 0; x < totalSize.x; ++x) {
            std::cout << debugSlice[x + y * totalSize.x];
        }
        std::cout << '\n';
    }
#endif

    // Check that the global block has the same information for each node that the node also has
    for (ind nodeIdx = 0; nodeIdx < numNodes.prod(); ++nodeIdx) {
        // Red in Local, Gray in Global
        // Green in Global, Gray in Local
    }

    float timeElapsed = timer.ElapsedTimeAndReset();
    std::cout << "Loaded and sorted data in " << timeElapsed << " seconds." << std::endl;

    for (float currentH = hMax; currentH >= hMin - 1e-5; currentH -= hStep) {
        for (auto itLocalBlock = localBlocks.begin(); itLocalBlock != localBlocks.end();
             itLocalBlock++) {
            itLocalBlock->doWatershed(currentH);
            itLocalBlock->sendData();
        }

        globalBlock.receiveData();
        globalBlock.doWatershed(currentH);
        globalBlock.sendData();

        for (auto itLocalBlock = localBlocks.begin(); itLocalBlock != localBlocks.end();
             itLocalBlock++) {
            itLocalBlock->receiveData();
        }

        std::cout << currentH << "/ " << hStep << "\t - " << globalBlock.numClustersCombined()
                  << std::endl;
        h.push_back(currentH);
#ifndef NDEBUG
        groundtruth->doWatershed(currentH);

        if (globalBlock.numClustersCombined() != groundtruth->numClustersCombined() ||
            globalBlock.totalVolumeCombined() != groundtruth->totalVolumeCombined()) {

            if (globalBlock.numClustersCombined() != groundtruth->numClustersCombined())
                std::cout << "Number of clusters is " << globalBlock.numClustersCombined() << '\\'
                          << groundtruth->numClustersCombined() << std::endl;
            if (globalBlock.totalVolumeCombined() != groundtruth->totalVolumeCombined())
                std::cout << "Total volume is " << ind(globalBlock.totalVolumeCombined()) << '\\'
                          << ind(groundtruth->totalVolumeCombined()) << std::endl;

            // Compate against green
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
                    std::cout << "Green part "
                              << "\t" << currentStat.first << ":\tcorrect " << truthIt->second
                              << " != " << currentStat.second << std::endl;
                }
            }

            for (auto itLocalBlock = localBlocks.begin(); itLocalBlock != localBlocks.end();
                 itLocalBlock++) {
                // TODO: Check red and white for each node
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

#ifndef NDEBUG
    delete groundtruth;
#endif  // NDEBUG
}

void whatershedMultipleRanks(int currProcess, vec3i numNodes, vec3i blockSize, vec3i blockOffset,
                             vec3i totalSize, float hMin, float hMax, float hStep,
                             std::vector<float>& h, std::vector<ind>& numClusters,
                             std::vector<float>& maxVolumes, std::vector<float>& totalVolumes,
                             PerformanceTimer& timer) {
#ifndef NDEBUG
    // Keeps groundtruth block to compare against
    LocalBlock* groundtruth = LocalBlock::makeGroundtruth(totalSize, vec3i(0), totalSize);
#endif  // NDEBUG

    // Master process
    if (currProcess == 0) {

        GlobalBlock globalBlock(blockSize, totalSize, numNodes);

        float timeElapsed = timer.ElapsedTimeAndReset();
        std::cout << "Rank " << currProcess << ": Loaded and sorted data in " << timeElapsed
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
            bool correct =
                globalBlock.numClustersCombined() == groundtruth->numClustersCombined() &&
                globalBlock.totalVolumeCombined() == groundtruth->totalVolumeCombined();

            // Let the other processes know something is wrong
            MPI_Request* requests = new MPI_Request[numNodes.prod()];
            for (ind p = 0; p < numNodes.prod(); p++) {
                ind processIndex = p + 1;
                MPI_Isend(&correct, 1, MPI_CXX_BOOL, processIndex, MPICommunication::ERRORFLAG,
                          MPI_COMM_WORLD, &requests[p]);
            }
            MPI_Waitall(numNodes.prod(), requests, MPI_STATUSES_IGNORE);
            if (!correct) {

                if (globalBlock.numClustersCombined() != groundtruth->numClustersCombined())
                    std::cout << "Number of clusters is " << globalBlock.numClustersCombined()
                              << '\\' << groundtruth->numClustersCombined() << std::endl;
                if (globalBlock.totalVolumeCombined() != groundtruth->totalVolumeCombined())
                    std::cout << "Total volume is " << ind(globalBlock.totalVolumeCombined())
                              << '\\' << ind(groundtruth->totalVolumeCombined()) << std::endl;

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
                    ind processIndex = p + 1;
                    std::cout << "0: Receiving red indices from process  " << processIndex
                              << std::endl;
                    redIndices.clear();
                    MPICommunication::RecvVectorUknownSize(redIndices, processIndex,
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
        float timeElapsed = timer.ElapsedTimeAndReset();
        std::cout << "Rank " << currProcess << ": Loaded and sorted data in " << timeElapsed
                  << " seconds." << std::endl;

        for (float currentH = hMax; currentH >= hMin - 1e-5; currentH -= hStep) {
            localBlock.doWatershed(currentH);
            localBlock.sendData();
            localBlock.receiveData();
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
        }
    }

#ifndef NDEBUG
    delete groundtruth;
#endif  // NDEBUG
}

void watershedLocalWhiteRed(vec3i blockSize, vec3i blockOffset, vec3i totalSize, float hMin,
                            float hMax, float hStep, std::vector<float>& h,
                            std::vector<ind>& numClusters, std::vector<float>& maxVolumes,
                            std::vector<float>& totalVolumes, PerformanceTimer& timer) {

#ifndef NDEBUG
    // Keeps groundtruth block to compare against
    LocalBlock* groundtruth = LocalBlock::makeGroundtruth(blockSize, blockOffset, totalSize);
#endif  // !NDEBUG

    LocalBlock* localBlockWhiteRed =
        LocalBlock::makeWhiteRedTest(blockSize, blockOffset, totalSize);

    float timeElapsed = timer.ElapsedTimeAndReset();
    std::cout << "Loaded and sorted data in " << timeElapsed << " seconds." << std::endl;

    for (float currentH = hMax; currentH >= hMin - 1e-5; currentH -= hStep) {
        localBlockWhiteRed->doWatershed(currentH);
        // Send and receive models the receiving and sending process, but actually just cleans
        // up
        localBlockWhiteRed->sendData();
        localBlockWhiteRed->receiveData();
        std::cout << currentH << "/ " << hStep << "\t - "
                  << localBlockWhiteRed->numClustersCombined() << std::endl;

        h.push_back(currentH);
#ifndef NDEBUG
        groundtruth->doWatershed(currentH);
        // WhiteRed
        if (localBlockWhiteRed->numClustersCombined() != groundtruth->numClustersCombined() ||
            localBlockWhiteRed->totalVolumeCombined() != groundtruth->totalVolumeCombined()) {

            if (localBlockWhiteRed->numClustersCombined() != groundtruth->numClustersCombined())
                std::cout << "Number of clusters is " << localBlockWhiteRed->numClustersCombined()
                          << '\\' << groundtruth->numClustersCombined() << std::endl;
            if (localBlockWhiteRed->totalVolumeCombined() != groundtruth->totalVolumeCombined())
                std::cout << "Total volume is " << ind(localBlockWhiteRed->totalVolumeCombined())
                          << '\\' << ind(groundtruth->totalVolumeCombined()) << std::endl;

            // Compare newet additions by related cluster volume.
            auto groundStats = groundtruth->getVoluminaForAddedVertices(currentH + hStep);
            auto whiteRedStats = localBlockWhiteRed->getVoluminaForAddedVertices(currentH + hStep);

            assert(groundStats.size() == whiteRedStats.size() && "Test is incorrect.");
            for (int i = 0; i < groundStats.size(); ++i) {
                assert(groundStats[i].first == whiteRedStats[i].first && "Test is incorrect.");
                if (groundStats[i].second != whiteRedStats[i].second) {
                    std::cout << groundStats[i].first << ":\tcorrect " << groundStats[i].second
                              << " != " << whiteRedStats[i].second << std::endl;
                }
            }

            assert(localBlockWhiteRed->numClustersCombined() == groundtruth->numClusters() &&
                   "Groundtruth has other num clusters.");
            assert(localBlockWhiteRed->totalVolumeCombined() == groundtruth->totalVolume() &&
                   "Groundtruth has other total volume.");
        }

#endif  // !NDEBUG
        numClusters.push_back(localBlockWhiteRed->numClustersCombined());
        maxVolumes.push_back(localBlockWhiteRed->maxVolumeCombined());
        totalVolumes.push_back(localBlockWhiteRed->totalVolumeCombined());
    }

    delete localBlockWhiteRed;
#ifndef NDEBUG
    delete groundtruth;
#endif  // !NDEBUG
}

void watershedGlobalGreen(vec3i blockSize, vec3i blockOffset, vec3i totalSize, float hMin,
                          float hMax, float hStep, std::vector<float>& h,
                          std::vector<ind>& numClusters, std::vector<float>& maxVolumes,
                          std::vector<float>& totalVolumes, PerformanceTimer& timer) {

#ifndef NDEBUG
    // Keeps groundtruth block to compare against
    LocalBlock* groundtruth = LocalBlock::makeGroundtruth(blockSize, blockOffset, totalSize);
#endif  // !NDEBUG

    GlobalBlock* globalBlockGreen = GlobalBlock::makeGreenTest(blockSize, blockOffset, totalSize);

    float timeElapsed = timer.ElapsedTimeAndReset();
    std::cout << "Loaded and sorted data in " << timeElapsed << " seconds." << std::endl;

    for (float currentH = hMax; currentH >= hMin - 1e-5; currentH -= hStep) {
        globalBlockGreen->doWatershed(currentH);
        // No send/receive
        std::cout << currentH << "/ " << hStep << "\t - " << globalBlockGreen->numClustersCombined()
                  << std::endl;
        h.push_back(currentH);
#ifndef NDEBUG
        groundtruth->doWatershed(currentH);
        if (globalBlockGreen->numClustersCombined() != groundtruth->numClustersCombined() ||
            globalBlockGreen->totalVolumeCombined() != groundtruth->totalVolumeCombined()) {

            if (globalBlockGreen->numClustersCombined() != groundtruth->numClustersCombined())
                std::cout << "Number of clusters is " << globalBlockGreen->numClustersCombined()
                          << '\\' << groundtruth->numClustersCombined() << std::endl;
            if (globalBlockGreen->totalVolumeCombined() != groundtruth->totalVolumeCombined())
                std::cout << "Total volume is " << ind(globalBlockGreen->totalVolumeCombined())
                          << '\\' << ind(groundtruth->totalVolumeCombined()) << std::endl;

            // Compare newest additions by related cluster volume.
            auto groundStats = groundtruth->getVoluminaForAddedVertices(currentH + hStep);
            auto greenStats = globalBlockGreen->getVoluminaForAddedVertices(currentH + hStep);

            assert(groundStats.size() == greenStats.size() && "Test is incorrect.");
            for (int i = 0; i < groundStats.size(); ++i) {
                assert(groundStats[i].first == greenStats[i].first && "Test is incorrect.");
                if (groundStats[i].second != greenStats[i].second) {
                    std::cout << groundStats[i].first << ":\tcorrect " << groundStats[i].second
                              << " != " << greenStats[i].second << std::endl;
                }
            }

            assert(globalBlockGreen->numClustersCombined() == groundtruth->numClusters() &&
                   "Groundtruth has other num clusters.");
            assert(globalBlockGreen->totalVolumeCombined() == groundtruth->totalVolume() &&
                   "Groundtruth has other total volume.");
        }
#endif  // !NDEBUG
        numClusters.push_back(globalBlockGreen->numClustersCombined());
        maxVolumes.push_back(globalBlockGreen->maxVolumeCombined());
        totalVolumes.push_back(globalBlockGreen->totalVolumeCombined());
    }

    delete globalBlockGreen;
#ifndef NDEBUG
    delete groundtruth;
#endif  // !NDEBUG
}

void watershedWhiteRed(vec3i blockSize, vec3i blockOffset, vec3i totalSize, float hMin, float hMax,
                       float hStep, std::vector<float>& h, std::vector<ind>& numClusters,
                       std::vector<float>& maxVolumes, std::vector<float>& totalVolumes,
                       PerformanceTimer& timer) {
#ifndef NDEBUG
    // Keeps groundtruth block to compare against
    LocalBlock* groundtruth = LocalBlock::makeGroundtruth(blockSize, blockOffset, totalSize);
#endif  // !NDEBUG

    // WhiteRed
    LocalBlock* localBlockWhiteRed =
        LocalBlock::makeWhiteRedTest(blockSize, blockOffset, totalSize);
    GlobalBlock* globalBlockRed = GlobalBlock::makeWhiteRedTest(blockSize, blockOffset, totalSize);

    float timeElapsed = timer.ElapsedTimeAndReset();
    std::cout << "Loaded and sorted data in " << timeElapsed << " seconds." << std::endl;

    for (float currentH = hMax; currentH >= hMin - 1e-5; currentH -= hStep) {
        localBlockWhiteRed->doWatershed(currentH);
        localBlockWhiteRed->sendData();
        globalBlockRed->receiveData();
        globalBlockRed->doWatershed(currentH);
        globalBlockRed->sendData();
        localBlockWhiteRed->receiveData();
        std::cout << currentH << "/ " << hStep << "\t - " << globalBlockRed->numClustersCombined()
                  << std::endl;
        h.push_back(currentH);

#ifndef NDEBUG
        groundtruth->doWatershed(currentH);
        // WhiteRed
        if (globalBlockRed->numClustersCombined() != groundtruth->numClustersCombined() ||
            globalBlockRed->totalVolumeCombined() != groundtruth->totalVolumeCombined()) {

            if (globalBlockRed->numClustersCombined() != groundtruth->numClustersCombined())
                std::cout << "Number of clusters is " << globalBlockRed->numClustersCombined()
                          << '\\' << groundtruth->numClustersCombined() << std::endl;
            if (globalBlockRed->totalVolumeCombined() != groundtruth->totalVolumeCombined())
                std::cout << "Total volume is " << ind(globalBlockRed->totalVolumeCombined())
                          << '\\' << ind(groundtruth->totalVolumeCombined()) << std::endl;

            // Compare newet additions by related cluster volume.
            auto groundStats = groundtruth->getVoluminaForAddedVertices(currentH + hStep);
            auto localStats = localBlockWhiteRed->getVoluminaForAddedVertices(currentH + hStep);
            auto globalStats = globalBlockRed->getVoluminaForAddedVertices(currentH + hStep);

            assert(localStats.size() + globalStats.size() == groundStats.size() &&
                   "Test is incorrect.");
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
                        std::cout << "(Green part) \t" << currentStat.first << ":\tcorrect "
                                  << currentStat.second << " != " << globalIt->second << std::endl;
                    }
                }
                // Red or white
                else if (localIt != localStats.end()) {
                    vec3i dummy(-1, -1, -1);
                    ClusterID* id = localBlockWhiteRed->findClusterID(currentStat.first, dummy);
                    // Red
                    if (id->isGlobal()) {
                        if (localIt->second != 0) {
                            std::cout << "(Red part) \t correct 0 != " << localIt->second
                                      << std::endl;
                        }
                        double volRed = globalBlockRed->getClusterVolume(*id);
                        if (currentStat.second != volRed) {
                            std::cout << "(Red part) \t" << currentStat.first << ":\tcorrect "
                                      << currentStat.second << " != " << volRed << std::endl;
                        }
                    } else {
                        if (currentStat.second != localIt->second) {
                            std::cout << "(White part) \t" << currentStat.first << ":\tcorrect "
                                      << currentStat.second << " != " << localIt->second
                                      << std::endl;
                        }
                    }
                }
            }

            assert(globalBlockRed->numClustersCombined() == groundtruth->numClusters() &&
                   "Groundtruth has other num clusters.");
            assert(globalBlockRed->totalVolumeCombined() == groundtruth->totalVolume() &&
                   "Groundtruth has other total volume.");
        }
#endif  // !NDEBUG
        numClusters.push_back(globalBlockRed->numClustersCombined());
        maxVolumes.push_back(globalBlockRed->maxVolumeCombined());
        totalVolumes.push_back(globalBlockRed->totalVolumeCombined());
    }

    delete localBlockWhiteRed;
    delete globalBlockRed;
#ifndef NDEBUG
    delete groundtruth;
#endif  // !NDEBUG
}

void watershedWhiteRedGreen(vec3i blockSize, vec3i blockOffset, vec3i totalSize, float hMin,
                            float hMax, float hStep, std::vector<float>& h,
                            std::vector<ind>& numClusters, std::vector<float>& maxVolumes,
                            std::vector<float>& totalVolumes, PerformanceTimer& timer) {

#ifndef NDEBUG
    // Keeps groundtruth block to compare against
    LocalBlock* groundtruth = LocalBlock::makeGroundtruth(blockSize, blockOffset, totalSize);
#endif  // !NDEBUG

    LocalBlock* localBlockWhiteRedGreen =
        LocalBlock::makeWhiteRedGreenTest(blockSize, blockOffset, totalSize);
    GlobalBlock* globalBlockRedGreen =
        GlobalBlock::makeWhiteRedGreenTest(blockSize, blockOffset, totalSize);

    float timeElapsed = timer.ElapsedTimeAndReset();
    std::cout << "Loaded and sorted data in " << timeElapsed << " seconds." << std::endl;

    for (float currentH = hMax; currentH >= hMin - 1e-5; currentH -= hStep) {
        localBlockWhiteRedGreen->doWatershed(currentH);
        localBlockWhiteRedGreen->sendData();
        globalBlockRedGreen->receiveData();
        globalBlockRedGreen->doWatershed(currentH);
        globalBlockRedGreen->sendData();
        localBlockWhiteRedGreen->receiveData();
        std::cout << currentH << "/ " << hStep << "\t - "
                  << globalBlockRedGreen->numClustersCombined() << std::endl;
        h.push_back(currentH);

#ifndef NDEBUG
        groundtruth->doWatershed(currentH);
        if (globalBlockRedGreen->numClustersCombined() != groundtruth->numClustersCombined() ||
            globalBlockRedGreen->totalVolumeCombined() != groundtruth->totalVolumeCombined()) {

            if (globalBlockRedGreen->numClustersCombined() != groundtruth->numClustersCombined())
                std::cout << "Number of clusters is " << globalBlockRedGreen->numClustersCombined()
                          << '\\' << groundtruth->numClustersCombined() << std::endl;
            if (globalBlockRedGreen->totalVolumeCombined() != groundtruth->totalVolumeCombined())
                std::cout << "Total volume is " << ind(globalBlockRedGreen->totalVolumeCombined())
                          << '\\' << ind(groundtruth->totalVolumeCombined()) << std::endl;

            // Compare newet additions by related cluster volume.
            auto groundStats = groundtruth->getVoluminaForAddedVertices(currentH + hStep);
            auto localStats =
                localBlockWhiteRedGreen->getVoluminaForAddedVertices(currentH + hStep);
            auto globalStats = globalBlockRedGreen->getVoluminaForAddedVertices(currentH + hStep);

            assert(localStats.size() + globalStats.size() == groundStats.size() &&
                   "Test is incorrect.");
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
                        std::cout << "(Green part) \t" << currentStat.first << ":\tcorrect "
                                  << currentStat.second << " != " << globalIt->second << std::endl;
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
                            std::cout << "(Red part) \t correct 0 != " << localIt->second
                                      << std::endl;
                        }
                        double volRed = globalBlockRedGreen->getClusterVolume(*id);
                        if (currentStat.second != volRed) {
                            std::cout << "(Red part) \t" << currentStat.first << ":\tcorrect "
                                      << currentStat.second << " != " << volRed << std::endl;
                        }
                    } else {
                        if (currentStat.second != localIt->second) {
                            std::cout << "(White part) \t" << currentStat.first << ":\tcorrect "
                                      << currentStat.second << " != " << localIt->second
                                      << std::endl;
                        }
                    }
                }

                assert(globalBlockRedGreen->numClustersCombined() == groundtruth->numClusters() &&
                       "Groundtruth has other num clusters.");
                assert(globalBlockRedGreen->totalVolumeCombined() == groundtruth->totalVolume() &&
                       "Groundtruth has other total volume.");
            }
        }
#endif  // !NDEBUG
        numClusters.push_back(globalBlockRedGreen->numClustersCombined());
        maxVolumes.push_back(globalBlockRedGreen->maxVolumeCombined());
        totalVolumes.push_back(globalBlockRedGreen->totalVolumeCombined());
    }

    delete localBlockWhiteRedGreen;
    delete globalBlockRedGreen;
#ifndef NDEBUG
    delete groundtruth;
#endif  // !NDEBUG
}

// Input args:
// - path:        folder with data
// - rms:         name of rms file
// - xT, yT, zT:  total size
// - xB, yB, zB:  block size
// - s:           mode
int main(int argc, char** argv) {

    // Directory path and sizes from args
    if (argc < 10) {
        std::cerr << "Not enough arguments.\n";
        return 1;
    }

    ind mode = atoi(argv[9]);

    // Test cases -> no need to parse other info
    if (mode == TEST_COMMUNICATION) {
        std::cout << "Testing MPI Communication." << std::endl;
        testingMPIVectors();
        return 0;
    } else if (mode == TEST_MERGING) {
        std::cout << "Testing Cluster Merges." << std::endl;
        testClusterMerges();
        return 0;
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

#ifndef SINGLENODE
    // One master node and one node for each block
    if (numNodes.prod() + 1 > numProcesses) {
        std::cerr << "Too few nodes. Needing " << numNodes.x << " x " << numNodes.y << " x "
                  << numNodes.z << " and one additional master" << std::endl;
        MPI_Finalize();
        return 1;
    }
#else
    if (numProcesses > 1) {
        std::cerr << "Using single node mode with " << numProcesses << " processes. " << std::endl;
        MPI_Finalize();
        return 2;
    }
#endif  // SINGLENODE

    int currProcess;
    MPI_Comm_rank(MPI_COMM_WORLD, &currProcess);

    // Process 0 is master rank -> -1
    vec3i idxNode;
    if (numProcesses == 1) {
        idxNode = vec3i::fromIndexOfTotal(0, numNodes);
    } else {
        idxNode = vec3i::fromIndexOfTotal(currProcess - 1, numNodes);
    }
    vec3i blockOffset = blockSize * idxNode;
    blockSize = vec3i::min(totalSize, blockSize * (idxNode + 1)) - blockOffset;

#ifndef SINGLENODE
    // Print process info
    if (currProcess != 0 && mode == REAL)
        std::cout << "Rank " << currProcess << ", index " << idxNode << ", offset " << blockOffset
                  << ", size " << blockSize << std::endl;
#endif

    // TODO, put settings in command line arguments
    float hMin = 0.0;
    float hMax = 2;
    assert(hMax > hMin && "HMax needs to be larger than hMin.");
    int hSamples = 101;
    float hStep = (hMax - hMin) / (hSamples - 1);

    // Keep track of threshold h, number of components, volume largest component, volume
    // total
    std::vector<float> h;
    std::vector<ind> numClusters;
    std::vector<float> maxVolumes;
    std::vector<float> totalVolumes;

    if (currProcess == 0) {
        h.reserve(hSamples);
        numClusters.reserve(hSamples);
        maxVolumes.reserve(hSamples);
        totalVolumes.reserve(hSamples);
    }

    PerformanceTimer timer;
    timer.Reset();
    float timeElapsed;

    switch (mode) {
        case REAL:
#ifdef SINGLENODE
#ifndef COMMUNICATION
            std::cout << "Watershedding sequentially." << std::endl;
            watershedSequential(totalSize, vec3i(0), totalSize, hMin, hMax, hStep, h, numClusters,
                                maxVolumes, totalVolumes, timer);
#else
            std::cout << "Watershedding sequentially, but distributed." << std::endl;
            watershedParallelSingleRank(numNodes, blockSize, blockOffset, totalSize, hMin, hMax,
                                        hStep, h, numClusters, maxVolumes, totalVolumes, timer);
#endif  // COMMUNCATION
#else   // !SINGLENODE
#ifdef COMMUNICATION
            if (currProcess == 0)
                std::cout << "Watershedding parallel and distributed." << std::endl;
            whatershedMultipleRanks(currProcess, numNodes, blockSize, blockOffset, totalSize, hMin,
                                    hMax, hStep, h, numClusters, maxVolumes, totalVolumes, timer);
#else   // !COMMUNCATION
            std::cerr << "Cannot watershed on multiple nodes without communcation "
                         "enabled."
                      << std::endl;
            MPI_Finalize();
            return 1;
#endif  // COMMUNCATION
#endif  // SINGLENODE
            break;

        case TEST_LOCALONLY_WHITERED:
#ifndef SINGLENODE
            if (currProcess == 0)
                std::cerr << "Testing single node mode (Local, WhiteRed) is not intended for "
                             "multiple nodes."
                          << std::endl;
            MPI_Finalize();
            return 1;
#else  // SINGLENODE
#ifdef COMMUNICATION
            std::cerr << "Testing without communication locally (WhiteRed) cannot be used with "
                         "communication."
                      << std::endl;
            MPI_Finalize();
            return 1;
#else
            std::cout << "Watershedding locally (WhiteRed) without communication," << std::endl;
            watershedLocalWhiteRed(totalSize, vec3i(0), totalSize, hMin, hMax, hStep, h,
                                   numClusters, maxVolumes, totalVolumes, timer);
#endif  //! COMMUNICATION
#endif  // !SINGLENODE
            break;

        case TEST_GLOBALONLY_GREEN:
#ifndef SINGLENODE
            if (currProcess == 0)
                std::cerr << "Testing single node mode (Global, Green) is not intended for "
                             "multiple nodes."
                          << std::endl;
            MPI_Finalize();
            return 1;
#else   // SINGLENODE
            std::cout << "Watershedding locally (Green) without communication." << std::endl;
            watershedGlobalGreen(totalSize, vec3i(0), totalSize, hMin, hMax, hStep, h, numClusters,
                                 maxVolumes, totalVolumes, timer);
#endif  // !SINGLENODE
            break;

        case TEST_WHITERED:
#ifndef SINGLENODE
            if (currProcess == 0)
                std::cerr << "Testing single node mode (Global, Green) is not intended for "
                             "multiple nodes."
                          << std::endl;
            MPI_Finalize();
            return 1;
#else  // SINGLENODE
#ifndef COMMUNICATION
            std::cerr
                << "Testing communication locally (WhiteRed) cannot be used without communication."
                << std::endl;
            MPI_Finalize();
            return 1;
#else
            std::cout
                << "Watershedding locally (WhiteRed) with communication (Global without Green)."
                << std::endl;
            watershedWhiteRed(totalSize, vec3i(0), totalSize, hMin, hMax, hStep, h, numClusters,
                              maxVolumes, totalVolumes, timer);
#endif  // COMMUNICATION
#endif  // !SINGLENODE
            break;

        case TEST_WHITEREDGREEN:
#ifndef SINGLENODE
            if (currProcess == 0)
                std::cerr << "Testing single node mode (Global, Green) is not intended for "
                             "multiple nodes."
                          << std::endl;
            MPI_Finalize();
            return 1;
#else  // SINGLENODE
#ifndef COMMUNICATION
            std::cerr << "Testing communication locally (WhiteRedGreen) cannot be used without "
                         "communication."
                      << std::endl;
            MPI_Finalize();
            return 1;
#else
            std::cout
                << "Watershedding locally (WhiteRedGreen) with communication (Global with Green)."
                << std::endl;
            watershedWhiteRedGreen(totalSize, vec3i(0), totalSize, hMin, hMax, hStep, h,
                                   numClusters, maxVolumes, totalVolumes, timer);
#endif  // COMMUNICATION
#endif  // !SINGLENODE
            break;
        default:
            if (currProcess == 0)
                std::cerr << "Mode " << mode << " is not a valid mode." << std::endl;
            return 1;
            break;
    }

    // Only master process writes out statistics
    if (currProcess == 0) {
        timeElapsed = timer.ElapsedTimeAndReset();
        std::cout << "Watershedding took " << timeElapsed << " seconds." << std::endl;

        const int maxClusters = *(std::max_element(numClusters.cbegin(), numClusters.cend()));

        std::ofstream percFile;
        std::string fileName = "percolation.csv";

        percFile.open(fileName, std::ios::out);

        if (percFile.is_open()) {
            percFile << "H; Number of connected components; Maximum number of connected "
                        "components; "
                        "Number of connected components / Maximum number of connected "
                        "components;  "
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