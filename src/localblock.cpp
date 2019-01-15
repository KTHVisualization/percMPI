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
    // TODO: Load subblocks into data
    vec3i min = blockOffset;
    vec3i max = blockOffset + blockSize;
    // Data actually covered by this block:
    // If min is global min there is no red or green,
    // if not, offset by one (for red)
    // if max is global max, there is no red or green,
    // if not, offset by two (for red and green)
    for (ind d = 0; d < 3; ++d) {
        min[d] = min[d] > 0 ? min[d]++ : 0;
        max[d] = max[d] < totalSize[d] ? max[d] - 2 : totalSize[d];
    }
    block->LOLSubBlock = new UnionFindSubBlock<WhiteProcessor>(
        max - min, min, totalSize, *block,
        WhiteProcessor(block->LOLs, block->LOGs, block->RefPLOGs));
    block->LOLSubBlock->loadData();

    // TODO: create IDBlock for all SubBlocks so that Pointers are together (need to be send).
    block->MemoryLOG = nullptr;  // For now: prevent dtor fail.
    block->MemoryLOGSize = 0;

    return block;
}  // namespace perc

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
    MPI_Status status;
    int err;

    // TODO: Change this into a broadcast
    int numNewLOGs;
    ind startOfLocalPlog;
    std::vector<int> merges;
    err = MPI_Recv(&numNewLOGs, 1, MPI_INT, 0, MPICommunication::NUMNEWCLUSTERS, MPI_COMM_WORLD,
                   &status);
    err = MPI_Recv(&startOfLocalPlog, 1, MPI_INT, 0, MPICommunication::STARTOFPLOG, MPI_COMM_WORLD,
                   &status);
    err = MPICommunication::RecvVectorUknownSize(merges, 0, MPICommunication::MERGES,
                                                 MPI_COMM_WORLD, &status);

    // Add some incognito clusters of other compute nodes.
    LOGs.addClusters(startOfLocalPlog);

    // For each PLOG, add a new cluster and reference it.
    for (ClusterData& c : CommPLOGs) {
        vec3i cPos = vec3i::fromIndexOfTotal(c.Index.RawID, TotalSize);
        ClusterID newID = LOGs.addCluster();
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

    // Receive updated green blocks
    int counter = 0;
    for (int counter = 0; counter < GOGSubBlocks.size(); ++counter) {
        auto gogBlock = GOGSubBlocks[counter];
        // Should this potentially be Non-Blocking?
        err = MPI_Recv(gogBlock.PointerBlock.PointerBlock, gogBlock.blockSize().prod() * sizeof(ID),
                       MPI_BYTE, 0, MPICommunication::GREENPOINTERS & counter, MPI_COMM_WORLD,
                       &status);
    }

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

    MPI_Request requests[7];
    int err;

    // Send number of local Clusters, maxVolume and totalVolume
    int numClustersLocal = numClusters();
    err == MPI_Isend(&numClustersLocal, 1, MPI_INT, 0, MPICommunication::NUMCLUSTERS,
                     MPI_COMM_WORLD, &requests[0]);
    double maxVolumeLocal = maxVolume();
    err = MPI_Isend(&maxVolumeLocal, 1, MPI_DOUBLE, 0, MPICommunication::MAXVOLUME, MPI_COMM_WORLD,
                    &requests[1]);
    double totalVolumeLocal = totalVolume();
    err = MPI_Isend(&totalVolumeLocal, 1, MPI_DOUBLE, 0, MPICommunication::TOTALVOLUME,
                    MPI_COMM_WORLD, &requests[2]);

    // Send volumes for LOGs (Additional volume)
    const std::vector<double>& commVolumes = LOGs.volumes();
    if (commVolumes.size() > 0) {
        err = MPI_Isend(commVolumes.data(), commVolumes.size(), MPI_DOUBLE, 0,
                        MPICommunication::VOLUMES, MPI_COMM_WORLD, &requests[3]);
    } else {
        err =
            MPI_Isend(0, 0, MPI_DOUBLE, 0, MPICommunication::VOLUMES, MPI_COMM_WORLD, &requests[3]);
    }
    err = MPICommunication::IsendVector(CommPLOGs, 0, MPICommunication::PLOGS, MPI_COMM_WORLD,
                                        &requests[4]);

    // Send merges, nothing more to be done here with them
    std::vector<ClusterMerge>& merges = LOGs.Merges;
    err = MPICommunication::IsendVector(merges, 0, MPICommunication::MERGES, MPI_COMM_WORLD,
                                        &requests[5]);

    // Send updated red blocks
    err = MPI_Isend(MemoryLOG, MemoryLOGSize * sizeof(ID), MPI_BYTE, 0,
                    MPICommunication::REDPOINTERS, MPI_COMM_WORLD, &requests[6]);

    // TODO: Confirm that everything arrived
    // MPI_Waitall(7, requests, MPI_STATUS_IGNORE);*/
}

void LocalBlock::repointerMultipleMerges(const std::vector<ind>& connComps) {
    for (auto it = connComps.begin(); it != connComps.end(); ++it) {
        ind compSize = *it;
        ind ontoCluster = *(++it);
        for (ind c = 0; c < compSize - 1; ++c) {
            ClusterID fromCluster = *(++it);
            VertexID fromRep = LOGs.getCluster(fromCluster).Index;
            VertexID ontoRep = LOGs.getCluster(ontoCluster).Index;
            setID(fromRep, ontoRep);
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