#include "greenblock.h"

namespace perc {

GreenBlock::GreenBlock(const vec3i& blockSize, const vec3i& blockOffset, const vec3i& totalSize) {

    // TODO: Actually load green blocks in correct locations, and set up number of red blocks
    int numBlocks = 1;
    float fraction = 1.0 / numBlocks;
    vec3i subBlocksize = vec3i(ind(blockSize.x * fraction), blockSize.y, blockSize.z);

    for (int i = 0; i < numBlocks; i++) {
        auto newGOGBlock = new UnionFindSubBlock<GlobalProcessor>(
            subBlocksize, {subBlocksize.x * i, subBlocksize.y * i, subBlocksize.z * i}, totalSize,
            *this, GlobalProcessor(GOGs));
        GOGSubBlocks.push_back(newGOGBlock);
    }

    for (UnionFindSubBlock<GlobalProcessor>* gogBlock : GOGSubBlocks) {
        gogBlock->loadData();
    }
}

void GreenBlock::doWatershed(const double minVal) {
    for (UnionFindSubBlock<GlobalProcessor>* gogBlock : GOGSubBlocks) {
        gogBlock->doWatershed(minVal);
    }
    // TODO: Apply merge list
    std::vector<ind> clusterMergeList =
        ClusterMerge::mergeClusterAsList(ClusterMerge::mergeClustersFromLists({GOGs.Merges}));
    // For each cluster merge
    ind from, onto;
    GOGs.mergeClustersForReal(from, onto);
    // Do all of the messy repointering and representative merging

    // All merges are processed, clear for next step
    GOGs.Merges.clear();
}

ClusterID* GreenBlock::findClusterID(const vec3i& idx, vec3i& lastClusterID) {
    // TODO: Do this more cleverly, and for all types of subblocks:
    // One might want to do this more cleverly, especialy in the sheet tree.
    for (UnionFindSubBlock<GlobalProcessor>* gogBlock : GOGSubBlocks)
        if (gogBlock->contains(idx)) return gogBlock->findClusterID(idx, lastClusterID);
    for (UnionFindSubBlock<LocalGlobalProcessor>* logBlock : LOGSubBlocks) {
        if (logBlock->contains(idx)) return logBlock->findClusterID(idx, lastClusterID);
    }
    assert(false && "Can not find block containing this idx.");
    return nullptr;
}

// Sketch.
void GreenBlock::receiveData() {
    /* TODO:
     * Receive Merges, number of clusters, PLOGs, volumes for LOGs and PLOGs and updated red
     * blocks from all nodes Turn PLOGs into LOGs (record range per node), initialize with their
     * respective volume Update Volume of existing global clusters
     */
}
void GreenBlock::sendData() {
    /* TODO:
     * Send number of new Clusters, PLOG range, updated Merges and updated green blocks to each
     * node
     */
}

}  // namespace perc