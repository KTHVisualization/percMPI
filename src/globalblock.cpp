#include <cmath>
#include "globalblock.h"
#include "performancetimer.h"
#include "mpicommuncation.h"
#include "interfaceblockbuilder.h"

namespace perc {

GlobalBlock::GlobalBlock(const vec3i& blockSize, const vec3i& totalSize, const vec3i& numNodes)
    : BlockSize(blockSize), NumNodes(numNodes), UnionFindBlock(totalSize), GOGs(GLOBAL_LIST) {

#ifndef NDEBUG
    for (int n = 0; n < 3; ++n) {
        double s = static_cast<double>(totalSize[n]) / blockSize[n];
        assert(numNodes[n] == static_cast<int>(ceil(s)) && "Block size and num nodes disagree.");
    }
#endif

    PerformanceTimer timer;
    timer.Reset();

    ind totalBlockSize = 0;
    ind numBlocks = 0;

    // Number of little 3x3x3 blocks.
    // Disregard size: The slivers cutting through them will amount to exactly 27 cells.
    numBlocks += (numNodes - 1).prod();

    // The 3x3xS edges between them.
    for (ind dim = 0; dim < 3; ++dim) {
        vec3i numSticks = numNodes - 1;
        numSticks[dim]++;
        numBlocks += numSticks.prod();

        // Length of sticks without little blocks.
        vec3i sizeSticks = numSticks;
        sizeSticks[dim] = totalSize[dim] - 3 * (numNodes[dim] - 1);
        // Two slivers cut through each edge. This makes 6 intersections per perpendicualr plane,
        // hence 9-6 = 3 muliplier.
        totalBlockSize += sizeSticks.prod() * 3;

        // The slivers.
        vec3i sizeSliver = totalSize;
        // The number of slivers.
        sizeSliver[dim] = numNodes[dim] - 1;
        totalBlockSize += sizeSliver.prod();

        vec3i numSlivers = numNodes;
        numSlivers[dim]--;
        numBlocks += numSlivers.prod();
    }

    std::cout << "Loading " << numBlocks << " blocks, with memory: " << totalBlockSize*4.0/(1024*1024) << "Mb." << std::endl;
    GOGSubBlocks.reserve(numBlocks);
    MemoryGreen = new ID[totalBlockSize];
    MemoryGreenSize = totalBlockSize;

    ID* memOngoing = MemoryGreen;
    ind fancyblockid = 0;

    // Reusable variables.
    vec3i whiteOffset, whiteSize;
    ID* memoryRed;
    ind sizeRed;

    // Directions.
    std::vector<vec3i> directions;
    directions.resize(7);

    // Results for all nodes.
    std::vector<std::vector<ind>> neighbors(numNodes.prod());
    PerProcessData.resize(NumNodes.prod());
    MemoryRedSize = 0;

    // Assemble.
    std::vector<vec3i> neighBlocks;
    neighBlocks.reserve(7);
    for (ind z = 0; z < numNodes.z; ++z)
        for (ind y = 0; y < numNodes.y; ++y)
            for (ind x = 0; x < numNodes.x; ++x) {
                directions.clear();
                vec3i node(x, y, z);
                vec3i min = node * blockSize;
                vec3i potentialMax = (node + 1) * blockSize;
                vec3i max = vec3i::min(potentialMax, totalSize);

                for (ind dim = 0; dim < 3; ++dim)
                    if (potentialMax[dim] < totalSize[dim]) {
                        // Direction the green side lies at.
                        vec3i dir(0);
                        dir[dim] = 1;

                        // Doing this first ensures a sorted vector
                        directions.push_back(dir);

                        // Combine the new direction with all previous ones (-> exclude the one just
                        // pushed)
                        ind size = directions.size();
                        for (ind d = 0; d < size - 1; ++d)
                            directions.push_back(directions[d] + dir);
                    }

                memoryRed = nullptr;
                // Buildng red adjacent blocks.
                PerProcessData[node.toIndexOfTotal(numNodes)].LOGstartIndex = LOGSubBlocks.size();
                interfaceblockbuilder::buildRedBlocks<GrayProcessor>(
                    max - min, min, totalSize, LOGSubBlocks, memoryRed, sizeRed, whiteSize,
                    whiteOffset, *this, []() { return GrayProcessor(); });
                PerProcessData[node.toIndexOfTotal(numNodes)].MemoryLOG = memoryRed;
                PerProcessData[node.toIndexOfTotal(numNodes)].MemoryLOGSize = sizeRed;
                MemoryRedSize += sizeRed;

                for (vec3i& dir : directions) {
                    // Building green adjacent blocks.
                    GOGSubBlocks.push_back(interfaceblockbuilder::buildGreenBlock<GreenProcessor>(
                        dir, whiteSize, whiteOffset, totalSize, memOngoing, *this,
                        [this]() { return GreenProcessor(GOGs); }));

                    timer.Reset();
                    GOGSubBlocks.back().loadData();
		    std::cout << " Loaded #" << fancyblockid++ << " of size " << GOGSubBlocks.back().blockSize().prod() * 4.0 / (1024.0*1024.0) << " in " << timer.ElapsedTimeAndReset() 
<< std::endl; 
                    // Add to neighborhood lists.
                    // Current node.
                    neighbors[node.toIndexOfTotal(numNodes)].push_back(GOGSubBlocks.size() - 1);

                    // Respective neighboring nodes.
                    neighBlocks.clear();
                    for (ind dim = 0; dim < 3; ++dim)
                        if (dir[dim]) {
                            // Direction the green side lies at.
                            vec3i neighDir(0);
                            neighDir[dim] = 1;

                            // Doing this first ensures a sorted vector
                            neighBlocks.push_back(neighDir);

                            // Combine the new direction with all previous ones (-> exclude the one
                            // just pushed)
                            ind size = neighBlocks.size();
                            for (ind d = 0; d < size - 1; ++d)
                                neighBlocks.push_back(neighBlocks[d] + neighDir);
                        }

                    for (auto& n : neighBlocks)
                        neighbors[(node + n).toIndexOfTotal(numNodes)].push_back(
                            GOGSubBlocks.size() - 1);
                }
            }

    ReceivedMerges.resize(NumNodes.prod() + 1);

    for (ind nodeIdx = 0; nodeIdx < neighbors.size(); ++nodeIdx) {
        PerProcessData[nodeIdx].GreenAdjacent = neighbors[nodeIdx];
        /*std::cout << nodeIdx + 1 << ": ";
        for (auto gogIdx : PerProcessData[nodeIdx].GreenAdjacent) {
            std::cout << GOGSubBlocks[gogIdx].blockOffset() << "\t";
        }
        std::cout << std::endl;*/
        ReceivedMerges[nodeIdx] = &PerProcessData[nodeIdx].Merges;
    }

    ReceivedMerges[NumNodes.prod()] = &GOGs.Merges;

    assert(memOngoing == MemoryGreen + totalBlockSize &&
           "Allocated and filled memory do not match.");
}

GlobalBlock::GlobalBlock(const vec3i& totalSize)
    : UnionFindBlock(totalSize)
    , GOGs(GLOBAL_LIST)
    , NumNodes(vec3i(1, 1, 1))
    , BlockSize(totalSize)
    , NumNewClusters(0)
    , NumClustersLocal(0)
    , MaxVolumeLocal(0)
    , TotalVolumeLocal(0) {}

GlobalBlock* GlobalBlock::makeGreenTest(const vec3i& blockSize, const vec3i& blockOffset,
                                        const vec3i& totalSize) {
    GlobalBlock* block = new GlobalBlock(totalSize);
    block->GOGSubBlocks.reserve(1);
    block->GOGSubBlocks.emplace_back(blockSize, blockOffset, totalSize, *block,
                                     GreenProcessor(block->GOGs));
    block->GOGSubBlocks.begin()->loadData();

    return block;
}

GlobalBlock* GlobalBlock::makeWhiteRedTest(const vec3i& blockSize, const vec3i& blockOffset,
                                           const vec3i& totalSize) {
    GlobalBlock* block = new GlobalBlock(totalSize);

    vec3i min = blockOffset;
    vec3i max = blockOffset + blockSize;

    // Create 3 red slices: one left, two right.
    vec3i sliceSize = vec3i(blockSize.x, blockSize.y, 10);
    min[2] += sliceSize.z;
    max[2] -= sliceSize.z * 2;

    // Set up data for sending and receiving as well as red blocks
    block->ReceivedMerges.resize(1);

    InfoPerProcess dataPerProcess;
    block->ReceivedMerges[0] = &dataPerProcess.Merges;
    block->ReceivedMerges[1] = &(block->GOGs).Merges;

    // No green indices
    std::vector<ind> greenIndices;
    dataPerProcess.GreenAdjacent = greenIndices;

    dataPerProcess.MemoryLOGSize = sliceSize.prod() * 3;
    dataPerProcess.MemoryLOG = new ID[dataPerProcess.MemoryLOGSize];
    dataPerProcess.LOGstartIndex = 0;

    block->PerProcessData.push_back(dataPerProcess);

    block->LOGSubBlocks.reserve(3);
    // Left slice.
    block->LOGSubBlocks.emplace_back(sliceSize, blockOffset, totalSize, *block, GrayProcessor(),
                                     dataPerProcess.MemoryLOG);

    // Two slices right.
    block->LOGSubBlocks.emplace_back(sliceSize, vec3i(blockOffset.x, blockOffset.y, max.z),
                                     totalSize, *block, GrayProcessor(),
                                     dataPerProcess.MemoryLOG + sliceSize.prod());
    block->LOGSubBlocks.emplace_back(
        sliceSize, vec3i(blockOffset.x, blockOffset.y, max.z + sliceSize.z), totalSize, *block,
        GrayProcessor(), dataPerProcess.MemoryLOG + 2 * sliceSize.prod());

    return block;
}  // namespace perc

GlobalBlock* GlobalBlock::makeWhiteRedGreenTest(const vec3i& blockSize, const vec3i& blockOffset,
                                                const vec3i& totalSize) {
    GlobalBlock* block = new GlobalBlock(totalSize);

    vec3i min = blockOffset;
    vec3i max = blockOffset + blockSize;

    // Create 3 green slices: two left, one right.
    // Create 2 red slices: one left, one right,
    vec3i sliceSize = vec3i(blockSize.x, blockSize.y, 10);
    min[2] += sliceSize.z * 3;
    max[2] -= sliceSize.z * 2;

    // Set up data for sending and receiving as well as red blocks
    block->ReceivedMerges.resize(2);

    InfoPerProcess dataPerProcess;
    block->ReceivedMerges[0] = &dataPerProcess.Merges;
    block->ReceivedMerges[1] = &(block->GOGs).Merges;

    std::vector<ind> greenIndices = {1, 2};
    dataPerProcess.GreenAdjacent = greenIndices;

    dataPerProcess.MemoryLOGSize = sliceSize.prod() * 2;
    dataPerProcess.MemoryLOG = new ID[dataPerProcess.MemoryLOGSize];
    dataPerProcess.LOGstartIndex = 0;
    block->PerProcessData.push_back(dataPerProcess);

    block->LOGSubBlocks.reserve(2);
    // Left slice.
    block->LOGSubBlocks.emplace_back(sliceSize,
                                     vec3i(blockOffset.x, blockOffset.y, min.z - sliceSize.z),
                                     totalSize, *block, GrayProcessor(), dataPerProcess.MemoryLOG);

    // Right slice.
    block->LOGSubBlocks.emplace_back(sliceSize, vec3i(blockOffset.x, blockOffset.y, max.z),
                                     totalSize, *block, GrayProcessor(),
                                     dataPerProcess.MemoryLOG + sliceSize.prod());

    // Green data.
    block->GOGSubBlocks.reserve(3);

    // Two left slices.
    block->GOGSubBlocks.emplace_back(sliceSize, vec3i(blockOffset.x, blockOffset.y, blockOffset.z),
                                     totalSize, *block, GreenProcessor(block->GOGs));

    block->GOGSubBlocks.emplace_back(
        sliceSize, vec3i(blockOffset.x, blockOffset.y, blockOffset.z + sliceSize.z), totalSize,
        *block, GreenProcessor(block->GOGs));

    // Right slice.
    block->GOGSubBlocks.emplace_back(sliceSize,
                                     vec3i(blockOffset.x, blockOffset.y, max.z + sliceSize.z),
                                     totalSize, *block, GreenProcessor(block->GOGs));

    for (auto& green : block->GOGSubBlocks) green.loadData();

    return block;
}

void GlobalBlock::doWatershed(const double minVal) {

    ind numClusters = GOGs.numClusters();
    for (auto& gogBlock : GOGSubBlocks) {
        gogBlock.doWatershed(minVal);
    }
    // Number of clusters can only have increased (Merges only recorded)
    NumNewClusters += GOGs.numClusters() - numClusters;

    // Merges have only been recorded -> Do merges now
    // PerformanceTimer timer;
    // timer.Reset();
    // std::cout << "Merge ";
    Merges = ClusterMerge::mergeClusterAsList(ClusterMerge::mergeClustersFromLists(ReceivedMerges));
    // std::cout << "took " << timer.ElapsedTimeAndReset() << " seconds." << std::endl;

    // Merges are processed, clear for next step
    GOGs.Merges.clear();

    // Do the representatives merge and repointering first
    repointerMultipleMerges(Merges);
    // Merge Clusters and Representatives in the list (this process will be recomputed by all other
    // blocks that receive the merges)
    GOGs.mergeClusterFromList(Merges);

    checkConsistency();
}

ClusterID* GlobalBlock::findClusterID(const vec3i& idx, vec3i& lastClusterID) {
    // Figure out in which process the vertex lies
    int process = MPICommunication::computeProcess(idx, BlockSize, NumNodes);
    if (process < PerProcessData.size()) {
        InfoPerProcess& processData = PerProcessData[process];

        for (auto& gogIndex : processData.GreenAdjacent) {
            auto& gogBlock = GOGSubBlocks[gogIndex];
            if (GOGSubBlocks[gogIndex].contains(idx))
                return gogBlock.findClusterID(idx, lastClusterID);
        }
        for (ind logIndex = processData.LOGstartIndex;
             logIndex < processData.LOGstartIndex + 6 && logIndex < LOGSubBlocks.size();
             logIndex++) {
            auto& logBlock = LOGSubBlocks[logIndex];
            if (logBlock.contains(idx)) return logBlock.findClusterID(idx, lastClusterID);
        }
    } else {
        // Legacy: Go through all blocks
        for (auto& gogBlock : GOGSubBlocks)
            if (gogBlock.contains(idx)) return gogBlock.findClusterID(idx, lastClusterID);
        for (auto& logBlock : LOGSubBlocks)
            if (logBlock.contains(idx)) return logBlock.findClusterID(idx, lastClusterID);
    }

    assert(false && "Can not find block containing this idx.");
    return nullptr;
}

ID* GlobalBlock::setID(const vec3i& idx, const ID& id) {
    // One might want to do this more cleverly, especialy in the sheet tree.
    ID* ptr = nullptr;

    int process = MPICommunication::computeProcess(idx, BlockSize, NumNodes);
    if (process < PerProcessData.size()) {
        InfoPerProcess& processData = PerProcessData[process];

        for (auto& gogIndex : processData.GreenAdjacent) {
            auto& gogBlock = GOGSubBlocks[gogIndex];
            if (GOGSubBlocks[gogIndex].contains(idx)) {
                ptr = gogBlock.PointerBlock.getPointer(idx);
                break;
            }
        }
        if (!ptr) {
            for (ind logIndex = processData.LOGstartIndex;
                 logIndex < processData.LOGstartIndex + 6 && logIndex < LOGSubBlocks.size();
                 logIndex++) {
                auto& logBlock = LOGSubBlocks[logIndex];
                if (logBlock.contains(idx)) {
                    ptr = logBlock.PointerBlock.getPointer(idx);
                    break;
                }
            }
        }

    } else {
        // Legacy: Go through all blocks
        for (auto& gogBlock : GOGSubBlocks)
            if (gogBlock.contains(idx)) {
                ptr = gogBlock.PointerBlock.getPointer(idx);
                break;
            }

        if (!ptr) {
            for (auto& logBlock : LOGSubBlocks)
                if (logBlock.contains(idx)) {
                    ptr = logBlock.PointerBlock.getPointer(idx);
                    break;
                }
        }
    }

    assert(ptr && "Can not find block containing this idx.");
    *ptr = id;

    return ptr;
}

double GlobalBlock::getClusterVolume(ClusterID cluster) { return GOGs.getClusterVolume(cluster); }

void GlobalBlock::repointerMultipleMerges(const std::vector<ind>& connComps) {
    for (auto it = connComps.begin(); it != connComps.end(); ++it) {
        ind compSize = *it;
        ind onto = *(++it);
        for (ind c = 0; c < compSize - 1; ++c) {
            ind from = *(++it);
            const std::vector<GOG>& oldRepsFrom = GOGs.getRepresentatives(from);
            const std::vector<VertexID> newRepsFrom = GOGs.mergeRepresentatives(from, onto);
            assert(oldRepsFrom.size() == newRepsFrom.size() &&
                   "Size mismatch between old and new representatives");
            auto itOld = oldRepsFrom.begin();
            for (auto itNew = newRepsFrom.begin();
                 itNew != newRepsFrom.end() && itOld != oldRepsFrom.end(); itOld++, itNew++) {
                UnionFindSubBlock<GreenProcessor>* parent =
                    reinterpret_cast<UnionFindSubBlock<GreenProcessor>*>(itOld->ParentBlock);
                // Rep is a new rep for the merged cluster -> point to it
                if (*itNew == itOld->ID) {
                    parent->PointerBlock.setPointer(
                        vec3i::fromIndexOfTotal(itOld->ID.baseID(), parent->totalSize()), onto);
                }
                // There is already a rep in the same block -> point to the rep from onto
                else {
                    parent->PointerBlock.setPointer(
                        vec3i::fromIndexOfTotal(itOld->ID.baseID(), parent->totalSize()), *itNew);
                }
            }
        }
    }
}

void GlobalBlock::receiveData() {
    // Init / Reset datastructures
    NumNewClusters = 0;
    NumClustersLocal = 0;
    MaxVolumeLocal = 0;
    TotalVolumeLocal = 0;
    ind plogsAddedSoFar = 0;
    Merges.clear();
    int numGOGs = GOGs.volumes().size();
    std::vector<double> commVolumes(numGOGs);

    // Request for each process and boolean if it has been processed yet
    MPI_Request* requests = new MPI_Request[NumNodes.prod()];
    std::vector<int> toProcessRanks(NumNodes.prod());
    std::iota(toProcessRanks.begin(), toProcessRanks.end(), 0);
    // Memory for receiving the two vector sizes for each process
    ind* vectorSizes = new ind[2 * NumNodes.prod()];

    // Non-Blocking Receive of the sizes
    for (ind p = 0; p < NumNodes.prod(); p++) {
        ind processDataIndex = p;
        ind processIndex = p + 1;
        int rank = 0;
#ifdef SINGLENODE
        // There is just a single node we will be receiving from
        // Encode the theoretical process index within the tag
        rank = processIndex;
        processIndex = 0;
#endif  // SINGLENODE
        int rankTag = rank << MPICommunication::RANK_SHIFT;
        MPI_Irecv(&vectorSizes[2 * processDataIndex], 2, MPI_IND, processIndex,
                  MPICommunication::VECTORSIZES | rankTag, MPI_COMM_WORLD,
                  &requests[processDataIndex]);
    }

    // Until we have received everything from everyone
    while (!toProcessRanks.empty()) {
        // Has not been processed yet
        ind numElements = toProcessRanks.size();
        for (ind i = numElements - 1; i >= 0; --i) {
            ind p = toProcessRanks[i];
            ind processDataIndex = p;
            // Do we already have sized
            int received;
            MPI_Test(&requests[processDataIndex], &received, MPI_STATUS_IGNORE);

            if (received) {
                // Replace by last position of the vector
                toProcessRanks[i] = toProcessRanks[numElements - 1];
                toProcessRanks.pop_back();
                numElements--;

                // Receive everything from that node now
                ind processIndex = p + 1;
                int rank = 0;
#ifdef SINGLENODE
                // There is just a single node we will be receiving from
                // Encode the theoretical process index within the tag
                rank = processIndex;
                processIndex = 0;
#endif  // SINGLENODE
                int rankTag = rank << MPICommunication::RANK_SHIFT;

                MPI_Status status;

                // Receive PLOGs and add a new cluster for each and initialize with their
                // respective volume
                std::vector<ClusterData> commPLOGs;
                commPLOGs.resize(vectorSizes[2 * processDataIndex]);

                MPICommunication::handleError(MPI_Recv(
                    commPLOGs.data(), commPLOGs.size() * sizeof(ClusterData), MPI_BYTE,
                    processIndex, MPICommunication::PLOGS | rankTag, MPI_COMM_WORLD, &status));

                // Receive merges, nothing more to be done here with them
                std::vector<ClusterMerge>& merges = *ReceivedMerges[processDataIndex];
                merges.resize(vectorSizes[2 * processDataIndex + 1]);

                MPICommunication::handleError(MPI_Recv(
                    merges.data(), merges.size() * sizeof(ClusterMerge), MPI_BYTE, processIndex,
                    MPICommunication::MERGES | rankTag, MPI_COMM_WORLD, &status));

                // Receive updated red blocks
                MPICommunication::handleError(
                    MPI_Recv(PerProcessData[processDataIndex].MemoryLOG,
                             PerProcessData[processDataIndex].MemoryLOGSize * sizeof(ID), MPI_BYTE,
                             processIndex, MPICommunication::REDPOINTERS | rankTag, MPI_COMM_WORLD,
                             &status));

                // Record where PLOGS for this process start (as in: How many other PLOGs for
                // other processes where added before)
                PerProcessData[p].StartOfLocalPlog = plogsAddedSoFar;
                NumNewClusters += commPLOGs.size();
                plogsAddedSoFar += commPLOGs.size();
                for (auto cluster : commPLOGs) {
                    UnionFindSubBlock<GrayProcessor>* parentBlock;
                    vec3i cPos = vec3i::fromIndexOfTotal(cluster.Index.RawID, TotalSize);
                    for (auto& logBlock : LOGSubBlocks)
                        if (logBlock.contains(cPos)) {
                            parentBlock = &logBlock;
                            break;
                        }
                    assert(parentBlock && "No Block for PLOG found");
                    ClusterID newCluster = GOGs.addCluster(cluster.Volume);
                    parentBlock->PointerBlock.setPointer(cPos, newCluster);
                }

#ifndef COLLECTIVES
                // Receive number of local Clusters and maxVolume and totalVolume
                ind numClusters;
                MPICommunication::handleError(MPI_Recv(&numClusters, 1, MPI_IND, processIndex,
                                                       MPICommunication::NUMCLUSTERS | rankTag,
                                                       MPI_COMM_WORLD, &status));

                NumClustersLocal += numClusters;
                double maxVolume;
                MPICommunication::handleError(MPI_Recv(&maxVolume, 1, MPI_DOUBLE, processIndex,
                                                       MPICommunication::MAXVOLUME | rankTag,
                                                       MPI_COMM_WORLD, &status));
                MaxVolumeLocal = std::max(maxVolume, MaxVolumeLocal);
                double volume;
                MPICommunication::handleError(MPI_Recv(&volume, 1, MPI_DOUBLE, processIndex,
                                                       MPICommunication::TOTALVOLUME | rankTag,
                                                       MPI_COMM_WORLD, &status));
                TotalVolumeLocal += volume;

                // Receive volumes for LOGs (Additional volume) and update volumes (Should be
                // some collected for this too)
                MPICommunication::handleError(
                    MPI_Recv(commVolumes.data(), commVolumes.size(), MPI_DOUBLE, processIndex,
                             MPICommunication::VOLUMES | rankTag, MPI_COMM_WORLD, &status));

                ind counter = 0;
                for (double vol : commVolumes) {
                    if (vol != 0) {
                        GOGs.extendCluster(ClusterID(counter, false), vol);
                    }
                    counter++;
                }
#endif  // !COLLECTIVES
            }
        }
    }

#ifndef SINGLENODE
    // Cleanup our requests
    MPI_Waitall(NumNodes.prod(), requests, MPI_STATUS_IGNORE);
#endif  // SINGLENODE

    delete[] requests;
    delete[] vectorSizes;

#ifdef COLLECTIVES
    ind numClustersDummy = 0;
    MPICommunication::handleError(
        MPI_Reduce(&numClustersDummy, &NumClustersLocal, 1, MPI_IND, MPI_SUM, 0, MPI_COMM_WORLD));
    double maxVolumeDummy = 0.0;
    MPICommunication::handleError(
        MPI_Reduce(&maxVolumeDummy, &MaxVolumeLocal, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD));
    double totalVolumeDummy = 0.0;
    MPICommunication::handleError(MPI_Reduce(&totalVolumeDummy, &TotalVolumeLocal, 1, MPI_DOUBLE,
                                             MPI_SUM, 0, MPI_COMM_WORLD));

    std::vector<double> commVolumesDummy(numGOGs);

    std::fill(commVolumesDummy.begin(), commVolumesDummy.end(), 0.0);
    MPICommunication::handleError(MPI_Reduce(commVolumesDummy.data(), commVolumes.data(), numGOGs,
                                             MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD));

    ind counter = 0;
    for (double vol : commVolumes) {
        if (vol != 0) {
            GOGs.extendCluster(ClusterID(counter, false), vol);
        }
        counter++;
    }
#endif
}

void GlobalBlock::sendData() {

#ifdef COLLECTIVES
    // Broadcast number of new Clusters and updated Merges
    MPI_Bcast(&NumNewClusters, 1, MPI_IND, 0, MPI_COMM_WORLD);
    ind mergesSize = Merges.size();
    MPI_Bcast(&mergesSize, 1, MPI_IND, 0, MPI_COMM_WORLD);
    MPI_Bcast(Merges.data(), Merges.size() * sizeof(ind), MPI_BYTE, 0, MPI_COMM_WORLD);
#endif

    for (ind p = 0; p < NumNodes.prod(); p++) {
        ind processDataIndex = p;
        ind processIndex = p + 1;
        int rank = 0;
#ifdef SINGLENODE
        // There is just a single node we will be receiving from
        // but the tag encodes with process it would be
        rank = processIndex;
        processIndex = 0;

#endif  // SINGLENODE
        int rankTag = rank << MPICommunication::RANK_SHIFT;

        // NumNewClusters, Merges Vector, PLOG range,
        ind numMessages = 1 + PerProcessData[processDataIndex].GreenAdjacent.size();
#ifndef COLLECTIVES
        numMessages += 2;
#endif

        MPI_Request* requests = new MPI_Request[numMessages];
        int messageID = 0;
#ifndef COLLECTIVES
        MPICommunication::handleError(MPI_Isend(&NumNewClusters, 1, MPI_IND, processIndex,
                                                MPICommunication::NUMNEWCLUSTERS | rankTag,
                                                MPI_COMM_WORLD, &requests[messageID++]));
        MPICommunication::handleError(
            MPICommunication::IsendVector(Merges, processIndex, MPICommunication::MERGES | rankTag,
                                          MPI_COMM_WORLD, &requests[messageID++]));
#endif  // COLLECTIVES

        // Plog Range (This might want to be a scatter operation) (TODO)
        MPICommunication::handleError(MPI_Isend(
            &PerProcessData[processDataIndex].StartOfLocalPlog, 1, MPI_IND, processIndex,
            MPICommunication::STARTOFPLOG | rankTag, MPI_COMM_WORLD, &requests[messageID++]));

        // Updated green blocks
        const std::vector<ind>& greenIndices = PerProcessData[processDataIndex].GreenAdjacent;
        int counter = 0;
        for (ind id : greenIndices) {
            auto& gogBlock = GOGSubBlocks[id];
            MPICommunication::handleError(MPI_Isend(
                gogBlock.PointerBlock.PointerBlock, gogBlock.blockSize().prod() * sizeof(ID),
                MPI_BYTE, processIndex, MPICommunication::GREENPOINTERS | rankTag | counter,
                MPI_COMM_WORLD, &requests[messageID++]));
            counter++;
        }

#ifndef SINGLENODE
        MPI_Waitall(numMessages, requests, MPI_STATUS_IGNORE);
#endif  // SINGLENODE
    }
}  // namespace perc

void GlobalBlock::checkConsistency() const {
#ifndef NDEBUG
    for (auto& gog : GOGSubBlocks) gog.checkConsistency();

    for (auto& merge : GOGs.Merges) {
        assert(GOGs.getClusterVolume(merge.From) >= 0 &&
               "Cluster recorded to merge from does not exist.");
        assert(GOGs.getClusterVolume(merge.Onto) > 0 &&
               "Cluster recorded to merge onto does not exist.");
    }
#endif
}

ind GlobalBlock::memEstimate() const {
    ind memSize = GOGs.memEstimate();

    for (auto& log : LOGSubBlocks) memSize += log.memEstimate();
    for (auto& gog : GOGSubBlocks) memSize += gog.memEstimate();

    return memSize;
}

std::vector<std::pair<vec3i, double>> GlobalBlock::getVoluminaForAddedVertices(double maxVal) {
    std::vector<std::pair<vec3i, double>> result;

    for (auto& gog : GOGSubBlocks) {
        gog.getVoluminaForAddedVertices(maxVal, result);
    }

    std::sort(result.begin(), result.end(), [this](auto& a, auto& b) {
        return a.first.toIndexOfTotal(TotalSize) > b.first.toIndexOfTotal(TotalSize);
    });

    return result;
}

}  // namespace perc
