#include "whiteblock.h"
#include "localprocessor.h"

namespace perc {

#ifndef TEST_LOCAL

WhiteBlock::WhiteBlock(const vec3i& blockSize, const vec3i& blockOffset, const vec3i& totalSize)
    : TotalSize(totalSize), RefPLOGs(10000, &ClusterID::hash) {
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
    : TotalSize(totalSize), RefPLOGs(10000, &ClusterID::hash) {
    vec3i min = blockOffset;
    vec3i max = blockOffset + blockSize;

    // Create 3 red slices: one left, two right.
    vec3i sliceSize = vec3i(50, blockSize.y, blockSize.z);
    min[0] += sliceSize.x;
    max[0] -= sliceSize.x * 2;

    LOLSubBlock = new UnionFindSubBlock<LocalLocalProcessor>(
        max - min, min, totalSize, *this, LocalLocalProcessor(LOLs, LOGs, RefPLOGs));
    LOLSubBlock->loadData();

    MemoryLOG = new ID[blockSize.prod() * 3];

    LOGSubBlocks.reserve(4);
    // Left slice.
    LOGSubBlocks.emplace_back(sliceSize, blockOffset, totalSize, *this,
                              LocalGlobalProcessor(LOLs, LOGs, RefPLOGs), MemoryLOG);

    // Two slices right.
    LOGSubBlocks.emplace_back(sliceSize, vec3i(max.x, 0, 0), totalSize, *this,
                              LocalGlobalProcessor(LOLs, LOGs, RefPLOGs),
                              MemoryLOG + sliceSize.prod());
    LOGSubBlocks.emplace_back(sliceSize, vec3i(max.x + sliceSize.x, 0, 0), totalSize, *this,
                              LocalGlobalProcessor(LOLs, LOGs, RefPLOGs),
                              MemoryLOG + 2 * sliceSize.prod());
    for (auto& red : LOGSubBlocks) red.loadData();
}

#endif

void WhiteBlock::doWatershed(const double minVal) {
    // TODO: for each subblock....
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

// Sketch.
void WhiteBlock::receiveData() {
    /*  TODO: [ ] Receive merge sets
     *          [x]    Merge locally :tick:
     */
    std::vector<std::vector<ClusterMerge>> merges = {LOGs.Merges};
    auto graphs = ClusterMerge::mergeClustersFromLists(merges);
    auto listAsItWouldArriveFromGreenBlock = ClusterMerge::mergeClusterAsList(graphs);
    LOGs.mergeClusterFromList(listAsItWouldArriveFromGreenBlock);

    /*  TODO: [ ] Receive numNewLOGs and startOfLocalPLOGs
     *          [x]    Do LOGs.addCluster() that often
     *        [x] Go through CommPLOGS:
     *          [x] Add the representative to the respective LOG (Reveive where my PLOGS are)
     */
    ind numNewLOGs = RefPLOGs.size();  // TODO: Change to MPI recv.
    ind startOfLocalPlog = 0;          // TODO: Change to MPI recv.

    // Add some incognito clusters of other compute nodes.
    LOGs.addClusters(startOfLocalPlog);

    // For each PLOG, add a new cluster and reference it.
    for (ClusterData& c : CommPLOGs) {
        // Search for green subblock containing this PLOG.
        for (auto& log : LOGSubBlocks) {
            vec3i cPos = vec3i::fromIndexOfTotal(c.Index.RawID, TotalSize);
            ClusterID newID = LOGs.addCluster();
            if (log.contains(cPos)) {
                *log.PointerBlock.getPointer(cPos) = newID;
                break;
            }
        }
    }
    LOGs.addClusters(numNewLOGs - startOfLocalPlog - CommPLOGs.size());

    CommPLOGs.clear();
}
void WhiteBlock::sendData() {
    /*  TODO: [ ] Fill CommPLOGs:
     *          [x] Move clusters into CommPLOGs
     *          [x] Remove from LOL list, clear RefPLOGs
     *          [ ] Send over LOL TotalVolume and number of clusters
     *          [ ] Send over CommPLOGs
     */
    CommPLOGs.clear();

    for (ClusterID plog : RefPLOGs) {
        Cluster c = LOGs.getCluster(plog);
        CommPLOGs.emplace_back(c.Index, c.Volume);
        LOLs.removeCluster(plog);
    }
    RefPLOGs.clear();
}

}  // namespace perc