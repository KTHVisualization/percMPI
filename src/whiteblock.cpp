#include "whiteblock.h"
#include "localprocessor.h"

namespace perc {

WhiteBlock::WhiteBlock(const vec3i& blockSize, const vec3i& blockOffset, const vec3i& totalSize) {
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
    // TODO: create IDBlock for all SubBlocks so that Pointers are together (need to be send)
    LOLSubBlock = new UnionFindSubBlock<LocalLocalProcessor>(
        max - min, min, totalSize, *this, LocalLocalProcessor(LOLs, LOGs, RefPLOGs));
    LOLSubBlock->loadData();
}

void WhiteBlock::doWatershed(const double minVal) {
    // TODO: for each subblock....
    LOLSubBlock->doWatershed(minVal);
}

ClusterID* WhiteBlock::findClusterID(const vec3i& idx, vec3i& lastClusterID) {
    // TODO: Do this more cleverly, and for all types of subblocks:
    /*
    ClusterID* UnionFindBlock::findClusterID(const vec3i& idx, vec3i& lastClusterID) {
        // One might want to do this more cleverly, especialy in the sheet tree.
        for (UnionFindSubBlock* sub : SubBlocks)
            if (sub->contains(idx)) return sub->findClusterID(idx, lastClusterID);

        assert(false && "Can not find block containing this idx.");
        return nullptr;
    }
    */
    return nullptr;
}

// Sketch.
void WhiteBlock::receiveData() {
    /*  TODO: Receive merge sets
     *      Merge locally
     *
     *  TODO: Receive numNewLOGs
     *      Do LOGs.addCluster() that often
     *      Go through CommPLOGS:
     *          Add the representative to the respective LOG (Reveive where my PLOGS are)
     */
}
void WhiteBlock::sendData() {
    /*  TODO: Fill CommPLOGs:
     *      Move clusters into CommPLOGs
     *      Remove from LOL list, clear RefPLOGs
     *      Send over LOL TotalVolume and number of clusters
     *      Send over CommPLOGs
     */
}

}  // namespace perc