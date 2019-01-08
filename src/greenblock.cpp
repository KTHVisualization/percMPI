#include "greenblock.h"
#include "performancetimer.h"

namespace perc {

GreenBlock::GreenBlock(const vec3i& blockSize, const vec3i& blockOffset, const vec3i& totalSize)
    : UnionFindBlock(totalSize), GOGs(LOCAL_LIST) {

    // TODO: Actually load green blocks in correct locations, and set up number of red blocks
    int numBlocks = 1;
    float fraction = 1.0 / numBlocks;
    vec3i subBlocksize = vec3i(ind(blockSize.x * fraction), blockSize.y, blockSize.z);
    vec3i subOffset = vec3i(ind(blockSize.x * fraction), 0, 0);

    GOGSubBlocks.reserve(numBlocks);

    for (int i = 0; i < numBlocks; i++) {
        int offsetScale = numBlocks - i;
        auto newGOGBlock = new UnionFindSubBlock<GlobalProcessor>(
            subBlocksize, {subOffset.x * i, subOffset.y * i, subOffset.z * i}, totalSize, *this,
            GlobalProcessor(GOGs));
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

    // Merges have only been recorded
    // TODO: Include received merges here
    // PerformanceTimer timer;
    // timer.Reset();
    // std::cout << "Merge ";

    std::vector<ind> connComps =
        ClusterMerge::mergeClusterAsList(ClusterMerge::mergeClustersFromLists({GOGs.Merges}));

    // std::cout << "took " << timer.ElapsedTimeAndReset() << " seconds." << std::endl;

    // Do the representatives merge and repointering first
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
                UnionFindSubBlock<GlobalProcessor>* parent =
                    reinterpret_cast<UnionFindSubBlock<GlobalProcessor>*>(itOld->ParentBlock);
                parent->PointerBlock.setPointer(
                    vec3i::fromIndexOfTotal(itOld->ID.baseID(), parent->totalSize()), *itNew);
            }
        }
    }

    // Merge Clusters and Representatives in the list
    GOGs.mergeClusterFromList(connComps);
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

ID* GreenBlock::setID(const vec3i& idx, const ID& id) {
    // One might want to do this more cleverly, especialy in the sheet tree.
    ID* ptr = nullptr;
    for (UnionFindSubBlock<GlobalProcessor>* gogBlock : GOGSubBlocks)
        if (gogBlock->contains(idx)) {
            ptr = gogBlock->PointerBlock.getPointer(idx);
            break;
        }

    if (!ptr) {
        for (UnionFindSubBlock<LocalGlobalProcessor>* logBlock : LOGSubBlocks)
            if (logBlock->contains(idx)) {
                ptr = logBlock->PointerBlock.getPointer(idx);
                break;
            }
    }

    assert(ptr && "Can not find block containing this idx.");
    *ptr = id;

    return ptr;
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

    // Merges are processed, clear for next step
    GOGs.Merges.clear();

    checkConsistency();
}

void GreenBlock::checkConsistency() const {
#ifndef NDEBUG
    for (auto gog : GOGSubBlocks) gog->checkConsistency();

    for (auto& merge : GOGs.Merges) {
        assert(GOGs.getClusterVolume(merge.From) >= 0 &&
               "Cluster recorded to merge from does not exist.");
        assert(GOGs.getClusterVolume(merge.Onto) > 0 &&
               "Cluster recorded to merge onto does not exist.");
    }
#endif
}

}  // namespace perc