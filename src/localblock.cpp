#include "localblock.h"
#include "interfaceblockbuilder.h"
#include "mpicommuncation.h"

namespace perc {

LocalBlock::LocalBlock(const vec3i& blockSize, const vec3i& blockOffset, const vec3i& totalSize)
    : UnionFindBlock(totalSize)
    , RefPLOGs(10000, &ClusterID::hash)
    , LOLs(LOCAL_LIST)
    , LOGs(GLOBAL_LIST)
    , LOGSubBlocks() {

    auto constructor = [this]() -> RedProcessor { return RedProcessor(LOLs, LOGs, RefPLOGs); };
    vec3i whiteBlockSize, whiteBlockOffset;
    interfaceblockbuilder::buildRedBlocks<RedProcessor>(
        blockSize, blockOffset, totalSize, LOGSubBlocks, MemoryLOG, MemoryLOGSize, whiteBlockSize,
        whiteBlockOffset, *this, constructor);

    LOLSubBlock = new UnionFindSubBlock<WhiteProcessor>(
        whiteBlockSize, whiteBlockOffset, totalSize, *this, WhiteProcessor(LOLs, LOGs, RefPLOGs));

    LOLSubBlock->loadData();
    for (auto& red : LOGSubBlocks) red.loadData();
}

LocalBlock::LocalBlock(const vec3i& totalSize)
    : UnionFindBlock(totalSize)
    , RefPLOGs(10000, &ClusterID::hash)
    , LOLs(LOCAL_LIST)
    , LOGs(GLOBAL_LIST) {}

LocalBlock* LocalBlock::makeGroundtruth(const vec3i& blockSize, const vec3i& blockOffset,
                                        const vec3i& totalSize) {
    LocalBlock* block = new LocalBlock(totalSize);
    block->LOLSubBlock = new UnionFindSubBlock<WhiteProcessor>(
        blockSize, blockOffset, totalSize, *block,
        WhiteProcessor(block->LOLs, block->LOGs, block->RefPLOGs));
    block->LOLSubBlock->loadData();

    // Id Block for Red is empty
    block->MemoryLOG = nullptr;
    block->MemoryLOGSize = 0;

    return block;
}

LocalBlock* LocalBlock::makeWhiteRedTest(const vec3i& blockSize, const vec3i& blockOffset,
                                         const vec3i& totalSize) {
    LocalBlock* block = new LocalBlock(totalSize);

    vec3i min = blockOffset;
    vec3i max = blockOffset + blockSize;

    // Create 3 red slices: one left, two right.
    vec3i sliceSize = vec3i(blockSize.x, blockSize.y, 10);
    min[2] += sliceSize.z;
    max[2] -= sliceSize.z * 2;

    block->LOLSubBlock = new UnionFindSubBlock<WhiteProcessor>(
        max - min, min, totalSize, *block,
        WhiteProcessor(block->LOLs, block->LOGs, block->RefPLOGs));
    block->LOLSubBlock->loadData();

    block->MemoryLOGSize = sliceSize.prod() * 3;
    block->MemoryLOG = new ID[block->MemoryLOGSize];

    block->LOGSubBlocks.reserve(3);
    // Left slice.
    block->LOGSubBlocks.emplace_back(sliceSize, blockOffset, totalSize, *block,
                                     RedProcessor(block->LOLs, block->LOGs, block->RefPLOGs),
                                     block->MemoryLOG);

    // Two slices right.
    block->LOGSubBlocks.emplace_back(sliceSize, vec3i(blockOffset.x, blockOffset.y, max.z),
                                     totalSize, *block,
                                     RedProcessor(block->LOLs, block->LOGs, block->RefPLOGs),
                                     block->MemoryLOG + sliceSize.prod());
    block->LOGSubBlocks.emplace_back(
        sliceSize, vec3i(blockOffset.x, blockOffset.y, max.z + sliceSize.z), totalSize, *block,
        RedProcessor(block->LOLs, block->LOGs, block->RefPLOGs),
        block->MemoryLOG + 2 * sliceSize.prod());
    for (auto& red : block->LOGSubBlocks) red.loadData();

    return block;
}

LocalBlock* LocalBlock::makeWhiteRedGreenTest(const vec3i& blockSize, const vec3i& blockOffset,
                                              const vec3i& totalSize) {

    LocalBlock* block = new LocalBlock(totalSize);

    vec3i min = blockOffset;
    vec3i max = blockOffset + blockSize;

    // Create 3 green slices: two left, one right.
    // Create 2 red slices: one left, one right,
    vec3i sliceSize = vec3i(blockSize.x, blockSize.y, 10);
    min[2] += sliceSize.z * 3;
    max[2] -= sliceSize.z * 2;

    block->LOLSubBlock = new UnionFindSubBlock<WhiteProcessor>(
        max - min, min, totalSize, *block,
        WhiteProcessor(block->LOLs, block->LOGs, block->RefPLOGs));
    block->LOLSubBlock->loadData();

    block->MemoryLOGSize = sliceSize.prod() * 2;
    block->MemoryLOG = new ID[block->MemoryLOGSize];

    // Red data.
    block->LOGSubBlocks.reserve(2);
    // Left slice.
    block->LOGSubBlocks.emplace_back(
        sliceSize, vec3i(blockOffset.x, blockOffset.y, min.z - sliceSize.z), totalSize, *block,
        RedProcessor(block->LOLs, block->LOGs, block->RefPLOGs), block->MemoryLOG);

    // Right right.
    block->LOGSubBlocks.emplace_back(sliceSize, vec3i(blockOffset.x, blockOffset.y, max.z),
                                     totalSize, *block,
                                     RedProcessor(block->LOLs, block->LOGs, block->RefPLOGs),
                                     block->MemoryLOG + sliceSize.prod());
    for (auto& red : block->LOGSubBlocks) red.loadData();

    // Green data.
    block->GOGSubBlocks.reserve(2);

    // Left slice (only the right most one of the two left slices that the global block holds is
    // contained in this local block
    block->GOGSubBlocks.emplace_back(
        sliceSize, vec3i(blockOffset.x, blockOffset.y, blockOffset.z + sliceSize.z), totalSize,
        *block, GrayProcessor());

    // Right slice.
    block->GOGSubBlocks.emplace_back(sliceSize,
                                     vec3i(blockOffset.x, blockOffset.y, max.z + sliceSize.z),
                                     totalSize, *block, GrayProcessor());

    return block;
}

void LocalBlock::doWatershed(const double minVal) {
    LOLSubBlock->doWatershed(minVal);
    for (auto& log : LOGSubBlocks) log.doWatershed(minVal);
}

ClusterID* LocalBlock::findClusterID(const vec3i& idx, vec3i& lastClusterID) {
    // One might want to do this more cleverly, especialy in the sheet tree.
    if (LOLSubBlock->contains(idx)) return LOLSubBlock->findClusterID(idx, lastClusterID);
    for (auto& log : LOGSubBlocks)
        if (log.contains(idx)) return log.findClusterID(idx, lastClusterID);

    for (auto& gog : GOGSubBlocks)
        if (gog.contains(idx)) return gog.findClusterID(idx, lastClusterID);

    assert(false && "Can not find block containing this idx.");
    return nullptr;
}

ID* LocalBlock::setID(const vec3i& idx, const ID& id) {
    // One might want to do this more cleverly, especialy in the sheet tree.
    ID* ptr = nullptr;
    if (LOLSubBlock->contains(idx))
        ptr = LOLSubBlock->PointerBlock.getPointer(idx);
    else {
        for (auto& log : LOGSubBlocks)
            if (log.contains(idx)) {
                ptr = log.PointerBlock.getPointer(idx);
                break;
            }
    }

    // No need to look into GOGBlocks, we should not be writing there anyways

    assert(ptr && "Can not find block containing this idx.");
    *ptr = id;

    return ptr;
}

double LocalBlock::getClusterVolume(ClusterID cluster) {
    return cluster.isGlobal() ? LOGs.getClusterVolume(cluster) : LOLs.getClusterVolume(cluster);
}

void LocalBlock::receiveData() {
    int numNewLOGs;
    ind startOfLocalPlog;
    std::vector<int> merges;
#ifdef COMMUNICATION
    MPI_Status status;
    int err;
#ifdef COLLECTIVES
    // TODO: Broadcasts
#else
    err = MPI_Recv(&numNewLOGs, 1, MPI_INT, 0, MPICommunication::NUMNEWCLUSTERS, MPI_COMM_WORLD,
                   &status);
    err = MPI_Recv(&startOfLocalPlog, 1, MPI_INT, 0, MPICommunication::STARTOFPLOG, MPI_COMM_WORLD,
                   &status);
    err = MPICommunication::RecvVectorUknownSize(merges, 0, MPICommunication::MERGES,
                                                 MPI_COMM_WORLD, &status);
#endif  // COLLECTIVES
    // Receive updated green blocks
    int counter = 0;
    for (int counter = 0; counter < GOGSubBlocks.size(); ++counter) {
        auto& gogBlock = GOGSubBlocks[counter];
        // Should this potentially be Non-Blocking?
        err = MPI_Recv(gogBlock.PointerBlock.PointerBlock, gogBlock.blockSize().prod() * sizeof(ID),
                       MPI_BYTE, 0, MPICommunication::GREENPOINTERS & counter, MPI_COMM_WORLD,
                       &status);
    }
#else   // !COMMUNICATION
    numNewLOGs = CommPLOGs.size();
    startOfLocalPlog = 0;
    merges = ClusterMerge::mergeClusterAsList(ClusterMerge::mergeClustersFromLists({LOGs.Merges}));
#endif  // COMMUNICATION

    // Add some incognito clusters of other compute nodes.
    LOGs.addClusters(startOfLocalPlog);

    // For each PLOG, add a new cluster and reference it.
    for (ClusterData& c : CommPLOGs) {
        vec3i cPos = vec3i::fromIndexOfTotal(c.Index.RawID, TotalSize);
        ClusterID newID = LOGs.addCluster(c.Volume);
        // Add representative and repointer PLOG (now LOG)
        setID(cPos, newID);
        LOGs.setRepresentative(newID, c.Index);
    }
    // Add Remaining new Clusters
    LOGs.addClusters(numNewLOGs - startOfLocalPlog - CommPLOGs.size());
    LOGs.clearVolumesAndMerges();
    CommPLOGs.clear();

    // First change pointers, then merge (cluster representative information is lost on merge).
    repointerMultipleMerges(merges);
    LOGs.mergeClusterFromList(merges);

    checkConsistency();
}

void LocalBlock::sendData() {
    checkConsistency();

    CommPLOGs.clear();
    CommPLOGs.reserve(RefPLOGs.size());

    for (ClusterID plog : RefPLOGs) {
        Cluster c = LOLs.getCluster(plog);
        CommPLOGs.emplace_back(c.Index, c.Volume);
        LOLs.removeCluster(plog);
    }
    RefPLOGs.clear();

#ifdef COMMUNICATION
    ind numMessages = 7;
    MPI_Request* requests = new MPI_Request[numMessages];
    int err;
    ind messageId = 0;

    // Number of local Clusters, maxVolume and totalVolume
    int numClustersLocal = numClusters();
    err = MPI_Isend(&numClustersLocal, 1, MPI_INT, 0, MPICommunication::NUMCLUSTERS, MPI_COMM_WORLD,
                    &requests[messageId++]);
    double maxVolumeLocal = maxVolume();
    err = MPI_Isend(&maxVolumeLocal, 1, MPI_DOUBLE, 0, MPICommunication::MAXVOLUME, MPI_COMM_WORLD,
                    &requests[messageId++]);
    double totalVolumeLocal = totalVolume();
    err = MPI_Isend(&totalVolumeLocal, 1, MPI_DOUBLE, 0, MPICommunication::TOTALVOLUME,
                    MPI_COMM_WORLD, &requests[messageId++]);

    // Volumes for LOGs (Additional volume)
    const std::vector<double>& commVolumes = LOGs.volumes();
    err = MPI_Isend(commVolumes.data(), commVolumes.size(), MPI_DOUBLE, 0,
                    MPICommunication::VOLUMES, MPI_COMM_WORLD, &requests[messageId++]);

    // New global clusters in red
    err = MPICommunication::IsendVector(CommPLOGs, 0, MPICommunication::PLOGS, MPI_COMM_WORLD,
                                        &requests[messageId++]);

    // Merges, nothing more to be done here with them
    std::vector<ClusterMerge>& merges = LOGs.Merges;
    err = MPICommunication::IsendVector(merges, 0, MPICommunication::MERGES, MPI_COMM_WORLD,
                                        &requests[messageId++]);

    // Send updated red blocks
    err = MPI_Isend(MemoryLOG, MemoryLOGSize * sizeof(ID), MPI_BYTE, 0,
                    MPICommunication::REDPOINTERS, MPI_COMM_WORLD, &requests[messageId++]);

#ifndef SINGLENODE
    MPI_Waitall(7, requests, MPI_STATUS_IGNORE);
#endif  // SINGLENODE
#endif  // COMMUNCATION
}

void LocalBlock::repointerMultipleMerges(const std::vector<ind>& connComps) {
    for (auto it = connComps.begin(); it != connComps.end(); ++it) {
        ind compSize = *it;
        ClusterID ontoCluster(*(++it), false);
        for (ind c = 0; c < compSize - 1; ++c) {
            ClusterID fromCluster(*(++it), false);
            VertexID fromRep = LOGs.getCluster(fromCluster).Index;
            VertexID ontoRep = LOGs.getCluster(ontoCluster).Index;
            // The from rep may not be part of this node (-> Empty Rep)
            if (fromRep.RawID != -1) {
                // Onto Rep not part of this node, make fromRep the Representative
                if (ontoRep.RawID == -1) {
                    setID(fromRep, ontoCluster);
                    LOGs.setRepresentative(ontoCluster, fromRep);
                } else  // Both reps are part of this node
                {
                    assert(fromRep != ontoRep && "Self pointing rep");
                    setID(fromRep, ontoRep);
                }
            }
        }
    }
}

void LocalBlock::checkConsistency() const {
#ifndef NDEBUG
    LOLSubBlock->checkConsistency();
    for (auto& log : LOGSubBlocks) log.checkConsistency();

    for (auto& merge : LOGs.Merges) {
        assert(LOGs.getClusterVolume(merge.From) >= 0 &&
               "Cluster recorded to merge from does not exist.");
        assert(LOGs.getClusterVolume(merge.Onto) > 0 &&
               "Cluster recorded to merge onto does not exist.");
    }
#endif
}

std::vector<std::pair<vec3i, double>> LocalBlock::getVoluminaForAddedVertices(double maxVal) {
    std::vector<std::pair<vec3i, double>> result;

    LOLSubBlock->getVoluminaForAddedVertices(maxVal, result);
    for (auto& log : LOGSubBlocks) {
        log.getVoluminaForAddedVertices(maxVal, result);
    }

    std::sort(result.begin(), result.end(), [this](auto& a, auto& b) {
        return a.first.toIndexOfTotal(TotalSize) > b.first.toIndexOfTotal(TotalSize);
    });

    return result;
}

}  // namespace perc