#include "whiteblock.h"
#include "localprocessor.h"

namespace perc {

#ifndef TEST_LOCAL

WhiteBlock::WhiteBlock(const vec3i& blockSize, const vec3i& blockOffset, const vec3i& totalSize)
    : UnionFindBlock(totalSize)
    , RefPLOGs(10000, &ClusterID::hash)
    , LOLs(LOCAL_LIST)
    , LOGs(GLOBAL_LIST) {
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
    LOLSubBlock = new UnionFindSubBlock<LocalLocalProcessor>(
        max - min, min, totalSize, *this, LocalLocalProcessor(LOLs, LOGs, RefPLOGs));
    LOLSubBlock->loadData();

    // TODO: create IDBlock for all SubBlocks so that Pointers are together (need to be send).
    MemoryLOG = nullptr;  // For now: prevent dtor fail.
}

#else

WhiteBlock::WhiteBlock(const vec3i& blockSize, const vec3i& blockOffset, const vec3i& totalSize)
    : UnionFindBlock(totalSize)
    , RefPLOGs(10000, &ClusterID::hash)
    , LOLs(LOCAL_LIST)
    , LOGs(GLOBAL_LIST) {
    vec3i min = blockOffset;
    vec3i max = blockOffset + blockSize;

    // Create 3 red slices: one left, two right.
    vec3i sliceSize = vec3i(blockSize.x, blockSize.y, 10);
    min[2] += sliceSize.z;
    max[2] -= sliceSize.z * 2;

    LOLSubBlock = new UnionFindSubBlock<LocalLocalProcessor>(
        max - min, min, totalSize, *this, LocalLocalProcessor(LOLs, LOGs, RefPLOGs));
    LOLSubBlock->loadData();

    MemoryLOG = new ID[blockSize.prod() * 3];

    LOGSubBlocks.reserve(4);
    // Left slice.
    LOGSubBlocks.emplace_back(sliceSize, blockOffset, totalSize, *this,
                              LocalGlobalProcessor(LOLs, LOGs, RefPLOGs), MemoryLOG);

    // Two slices right.
    LOGSubBlocks.emplace_back(sliceSize, vec3i(0, 0, max.z), totalSize, *this,
                              LocalGlobalProcessor(LOLs, LOGs, RefPLOGs),
                              MemoryLOG + sliceSize.prod());
    LOGSubBlocks.emplace_back(sliceSize, vec3i(0, 0, max.z + sliceSize.z), totalSize, *this,
                              LocalGlobalProcessor(LOLs, LOGs, RefPLOGs),
                              MemoryLOG + 2 * sliceSize.prod());
    for (auto& red : LOGSubBlocks) red.loadData();
}

WhiteBlock::WhiteBlock(const vec3i& totalSize)
    : UnionFindBlock(totalSize)
    , RefPLOGs(10000, &ClusterID::hash)
    , LOLs(LOCAL_LIST)
    , LOGs(GLOBAL_LIST) {}

WhiteBlock* WhiteBlock::makeGroundtruth(const vec3i& blockSize, const vec3i& blockOffset,
                                        const vec3i& totalSize) {
    WhiteBlock* block = new WhiteBlock(totalSize);
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
    block->LOLSubBlock = new UnionFindSubBlock<LocalLocalProcessor>(
        max - min, min, totalSize, *block,
        LocalLocalProcessor(block->LOLs, block->LOGs, block->RefPLOGs));
    block->LOLSubBlock->loadData();

    // TODO: create IDBlock for all SubBlocks so that Pointers are together (need to be send).
    block->MemoryLOG = nullptr;  // For now: prevent dtor fail.

    return block;
}  // namespace perc

void WhiteBlock::doWatershed(const double minVal) {
    LOLSubBlock->doWatershed(minVal);
    for (auto& log : LOGSubBlocks) log.doWatershed(minVal);
}

ClusterID* WhiteBlock::findClusterID(const vec3i& idx, vec3i& lastClusterID) {
    // One might want to do this more cleverly, especialy in the sheet tree.
    if (LOLSubBlock->contains(idx)) return LOLSubBlock->findClusterID(idx, lastClusterID);
    for (auto& log : LOGSubBlocks)
        if (log.contains(idx)) return log.findClusterID(idx, lastClusterID);

    // TODO: Do the same for GOGs.

    assert(false && "Can not find block containing this idx.");
    return nullptr;
}

ID* WhiteBlock::setID(const vec3i& idx, const ID& id) {
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

    assert(ptr && "Can not find block containing this idx.");
    *ptr = id;

    return ptr;
}

double WhiteBlock::getClusterVolume(ClusterID cluster) {
    return cluster.isGlobal() ? LOGs.getClusterVolume(cluster) : LOLs.getClusterVolume(cluster);
}

// Sketch.
void WhiteBlock::receiveData() {
    /*  TODO: [ ] Receive merge sets
     *          [x]    Merge locally
     */
    std::vector<std::vector<ClusterMerge>> merges = {LOGs.Merges};
    auto graphs = ClusterMerge::mergeClustersFromLists(merges);
    auto listAsItWouldArriveFromGreenBlock = ClusterMerge::mergeClusterAsList(graphs);
    // First change pointers, then merge (cluster representative information is lost on merge).
    repointerMultipleMerges(listAsItWouldArriveFromGreenBlock);
    LOGs.mergeClusterFromList(listAsItWouldArriveFromGreenBlock);

    /*  TODO: [ ] Receive numNewLOGs and startOfLocalPLOGs
     *          [x]    Do LOGs.addCluster() that often
     *        [x] Go through CommPLOGS:
     *          [x] Add the representative to the respective LOG (Reveive where my PLOGS are)
     */
    ind numNewLOGs = CommPLOGs.size();  // TODO: Change to MPI recv.
    ind startOfLocalPlog = 0;           // TODO: Change to MPI recv.

    // Add some incognito clusters of other compute nodes.
    LOGs.addClusters(startOfLocalPlog);

    // For each PLOG, add a new cluster and reference it.
    for (ClusterData& c : CommPLOGs) {
        vec3i cPos = vec3i::fromIndexOfTotal(c.Index.RawID, TotalSize);
        ClusterID newID = LOGs.addCluster();

        setID(cPos, newID);
        LOGs.setRepresentative(newID, c.Index);
    }
    LOGs.addClusters(numNewLOGs - startOfLocalPlog - CommPLOGs.size());
    LOGs.clearVolumesAndMerges();
    CommPLOGs.clear();

    checkConsistency();
}

void WhiteBlock::sendData() {
    checkConsistency();

    /*  TODO: [ ] Fill CommPLOGs:
     *          [x] Move clusters into CommPLOGs
     *          [x] Remove from LOL list, clear RefPLOGs
     *          [ ] Send over LOL TotalVolume and number of clusters
     *          [ ] Send over CommPLOGs
     */
    CommPLOGs.clear();
    CommPLOGs.reserve(RefPLOGs.size());

    for (ClusterID plog : RefPLOGs) {
        Cluster c = LOLs.getCluster(plog);
        CommPLOGs.emplace_back(c.Index, c.Volume);
        LOLs.removeCluster(plog);
    }
    RefPLOGs.clear();
}

void WhiteBlock::repointerMultipleMerges(const std::vector<ind>& connComps) {
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

void WhiteBlock::checkConsistency() const {
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

}  // namespace perc