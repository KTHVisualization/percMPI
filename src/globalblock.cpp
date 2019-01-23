#include "globalblock.h"
#include "performancetimer.h"
#include "mpicommuncation.h"
#include "interfaceblockbuilder.h"

namespace perc {

GlobalBlock::GlobalBlock(const vec3i& blockSize, const vec3i& blockOffset, const vec3i& totalSize,
                         const vec3i& numNodes)
    : UnionFindBlock(totalSize)
    , GOGs(GLOBAL_LIST)
    , NumNodes(numNodes)
    , NumNewClusters(0)
    , NumClustersLocal(0)
    , MaxVolumeLocal(0)
    , TotalVolumeLocal(0) {

    // TODO: Actually load green blocks in correct locations, and set up number of red blocks
    int numBlocks = 0;  // Currently set to 0 for communications test
    float fraction = 1.0 / numBlocks;
    vec3i subBlocksize = vec3i(blockSize.x, blockSize.y, ind(blockSize.z * fraction));
    vec3i subOffset = vec3i(0, 0, ind(blockSize.z * fraction));

    GOGSubBlocks.reserve(numBlocks);

    for (int i = 0; i < numBlocks; i++) {
        int offsetScale = numBlocks - i;
        GOGSubBlocks.emplace_back(subBlocksize, subOffset * i, totalSize, *this,
                                  GreenProcessor(GOGs));
    }

    for (auto& gogBlock : GOGSubBlocks) {
        gogBlock.loadData();
    }

    // Set up data for sending and receiving as well as red blocks
    ReceivedMerges.resize(NumNodes.prod());

    for (ind p = 0; p < numNodes.prod(); ++p) {
        InfoPerProcess dataPerProcess;

        dataPerProcess.Merges = &ReceivedMerges[p];

        // TODO: Add correct greenIndices
        std::vector<ind> greenIndices;
        dataPerProcess.GreenIndices = greenIndices;

        auto constructor = [this]() -> GrayProcessor { return GrayProcessor(); };
        vec3i whiteBlockSize, whiteBlockOffset;

        interfaceblockbuilder::buildRedBlocks<GrayProcessor>(
            blockSize, blockOffset, totalSize, LOGSubBlocks, dataPerProcess.MemoryLOG,
            dataPerProcess.MemoryLOGSize, whiteBlockSize, whiteBlockOffset, *this, constructor);

        PerProcessData.push_back(dataPerProcess);
    }
}

GlobalBlock::GlobalBlock(const vec3i& totalSize)
    : UnionFindBlock(totalSize)
    , GOGs(GLOBAL_LIST)
    , NumNodes(vec3i(1, 1, 1))
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
    block->ReceivedMerges.resize(block->NumNodes.prod());

    InfoPerProcess dataPerProcess;
    dataPerProcess.Merges = &block->ReceivedMerges[0];

    // No green indices
    std::vector<ind> greenIndices;
    dataPerProcess.GreenIndices = greenIndices;

    dataPerProcess.MemoryLOGSize = sliceSize.prod() * 3;
    dataPerProcess.MemoryLOG = new ID[dataPerProcess.MemoryLOGSize];

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
    block->ReceivedMerges.resize(1);

    InfoPerProcess dataPerProcess;
    dataPerProcess.Merges = &block->ReceivedMerges[0];

    std::vector<ind> greenIndices = {1, 2};
    dataPerProcess.GreenIndices = greenIndices;

    dataPerProcess.MemoryLOGSize = sliceSize.prod() * 2;
    dataPerProcess.MemoryLOG = new ID[dataPerProcess.MemoryLOGSize];
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
    ReceivedMerges.push_back(GOGs.Merges);
    // PerformanceTimer timer;
    // timer.Reset();
    // std::cout << "Merge ";
    Merges = ClusterMerge::mergeClusterAsList(ClusterMerge::mergeClustersFromLists(ReceivedMerges));
    // std::cout << "took " << timer.ElapsedTimeAndReset() << " seconds." << std::endl;

    // Merges are processed, clear for next step (remove from received again)
    GOGs.Merges.clear();
    ReceivedMerges.pop_back();

    // Do the representatives merge and repointering first
    repointerMultipleMerges(Merges);
    // Merge Clusters and Representatives in the list (this process will be recomputed by all other
    // blocks that receive the merges)
    GOGs.mergeClusterFromList(Merges);

    checkConsistency();
}

ClusterID* GlobalBlock::findClusterID(const vec3i& idx, vec3i& lastClusterID) {
    // TODO: Do this more cleverly, and for all types of subblocks:
    // One might want to do this more cleverly, especialy in the sheet tree.
    for (auto& gogBlock : GOGSubBlocks)
        if (gogBlock.contains(idx)) return gogBlock.findClusterID(idx, lastClusterID);
    for (auto& logBlock : LOGSubBlocks) {
        if (logBlock.contains(idx)) return logBlock.findClusterID(idx, lastClusterID);
    }
    assert(false && "Can not find block containing this idx.");
    return nullptr;
}

ID* GlobalBlock::setID(const vec3i& idx, const ID& id) {
    // One might want to do this more cleverly, especialy in the sheet tree.
    ID* ptr = nullptr;
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
            int from = *(++it);
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

    for (ind p = 0; p < NumNodes.prod(); p++) {

        ind processDataIndex = p;
        ind processIndex = p + 1;
#ifdef SINGLENODE
        // There is just a single node we will be receiving from
        processIndex = 0;
#endif  // SINGLENODE

        ind numMessages = 7;
        MPI_Status status;
        int err;

        // Receive number of local Clusters and maxVolume and totalVolume (TODO: Might want to be
        // gathers)
        int numClusters;
        err = MPI_Recv(&numClusters, 1, MPI_INT, processIndex, MPICommunication::NUMCLUSTERS,
                       MPI_COMM_WORLD, &status);
        NumClustersLocal += numClusters;
        double maxVolume;
        err = MPI_Recv(&maxVolume, 1, MPI_DOUBLE, processIndex, MPICommunication::MAXVOLUME,
                       MPI_COMM_WORLD, &status);
        MaxVolumeLocal = std::max(maxVolume, MaxVolumeLocal);
        double volume;
        err = MPI_Recv(&volume, 1, MPI_DOUBLE, processIndex, MPICommunication::TOTALVOLUME,
                       MPI_COMM_WORLD, &status);
        TotalVolumeLocal += volume;

        // Receive volumes for LOGs (Additional volume) and update volumes (Should be some collected
        // for this too)
        std::vector<double> commVolumes(GOGs.volumes().size());
        err = MPI_Recv(commVolumes.data(), commVolumes.size(), MPI_DOUBLE, processIndex,
                       MPICommunication::VOLUMES, MPI_COMM_WORLD, &status);

        ind counter = 0;
        for (double vol : commVolumes) {
            if (vol != 0) {
                GOGs.extendCluster(ClusterID(counter, false), vol);
            }
            counter++;
        }

        // Receive PLOGs and add a new cluster for each and initialize with their respective
        // volume
        std::vector<ClusterData> commPLOGs;
        err = MPICommunication::RecvVectorUknownSize(commPLOGs, p, MPICommunication::PLOGS,
                                                     MPI_COMM_WORLD, &status);

        // Receive merges, nothing more to be done here with them
        std::vector<ClusterMerge>& merges = ReceivedMerges[processDataIndex];
        err = MPICommunication::RecvVectorUknownSize(merges, processIndex, MPICommunication::MERGES,
                                                     MPI_COMM_WORLD, &status);

        // Receive updated red blocks
        err = MPI_Recv(PerProcessData[processDataIndex].MemoryLOG,
                       PerProcessData[processDataIndex].MemoryLOGSize * sizeof(ID), MPI_BYTE,
                       processIndex, MPICommunication::REDPOINTERS, MPI_COMM_WORLD, &status);

        // Record where PLOGS for this process start (as in: How many other PLOGs for other
        // processes where added before)
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
            ClusterID newCluster = GOGs.addCluster(cluster.Index, cluster.Volume, parentBlock);
            parentBlock->PointerBlock.setPointer(cPos, newCluster);
        }
    }
}

void GlobalBlock::sendData() {

#ifdef COLLECTIVES
    // Broadcast number of new Clusters and updated Merges
    MPI_Bcast(&NumNewClusters, 1, MPI_INT, 0, MPI_COMM_WORLD);
    ind mergesSize = Merges.size();
    MPI_Bcast(&mergesSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&Merges.data(), Merges.size() * sizeof(ind), MPI_BYTE, 0, MPI_COMM_WORLD);
#endif

    for (ind p = 0; p < NumNodes.prod(); p++) {
        ind processDataIndex = p;
        ind processIndex = p + 1;
#ifdef SINGLENODE
        // There is just a single node we will be receiving from
        processIndex = 0;
#endif  // SINGLENODE

        // NumNewClusters, Merges Vector, PLOG range,
        ind numMessages = 1 + PerProcessData[processDataIndex].GreenIndices.size();
#ifndef COLLECTIVES
        numMessages += 2;
#endif

        MPI_Request* requests = new MPI_Request[numMessages];

        // Updated green blocks
        const std::vector<ind>& greenIndices = PerProcessData[processDataIndex].GreenIndices;
        ind messageID = 0;
        for (ind id : greenIndices) {
            auto& gogBlock = GOGSubBlocks[id];
            MPI_Isend(gogBlock.PointerBlock.PointerBlock, gogBlock.blockSize().prod() * sizeof(ID),
                      MPI_BYTE, processIndex, MPICommunication::GREENPOINTERS & messageID,
                      MPI_COMM_WORLD, &requests[messageID++]);
        }

        // Plog Range (This might want to be a scatter operation) (TODO)
        MPI_Isend(&PerProcessData[processDataIndex].StartOfLocalPlog, 1, MPI_INT, processIndex,
                  MPICommunication::STARTOFPLOG, MPI_COMM_WORLD, &requests[messageID++]);

#ifndef COLLECTIVES
        MPI_Isend(&NumNewClusters, 1, MPI_INT, processIndex, MPICommunication::NUMNEWCLUSTERS,
                  MPI_COMM_WORLD, &requests[messageID++]);
        MPICommunication::IsendVector(Merges, processIndex, MPICommunication::MERGES,
                                      MPI_COMM_WORLD, &requests[messageID++]);
#endif  // COLLECTIVES

#ifndef SINGLENODE
        // This should free our requests as well??? -> Would mean blocking after each process, might
        // not be what we want
        MPI_Waitall(numMessages, requests, MPI_STATUS_IGNORE);
#endif  // SINGLENODE
    }
}

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