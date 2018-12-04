#include "unionfindsubblock.h"
namespace perc {

template <typename ClusterProcessor>
UnionFindSubBlock<ClusterProcessor>::UnionFindSubBlock(const vec3i& size, const vec3i& offset,
                                                       const vec3i& total, UnionFindBlock& parent,
                                                       ClusterProcessor&& neighProcessor)
    : PointerBlock(size, offset, total)
    , Parent(parent)
    , CurrentWatershedIndex(-1)
    , NeighborProcessor(neighProcessor)
    , Data(nullptr) {

    NeighborCache.reserve(6);
    NeighborProcessor.setParent(this);
}

template <typename ClusterProcessor>
void UnionFindSubBlock<ClusterProcessor>::loadData() {
    assert(!Data && "Data was already set.");
    Data = new DataBlock(PointerBlock.BlockSize, PointerBlock.BlockOffset, PointerBlock.TotalSize);

    CurrentWatershedIndex = 0;
    Data->sort();
}

template <typename ClusterProcessor>
void UnionFindSubBlock<ClusterProcessor>::doWatershed(double maxVal) {
    if (!Data) return;

    NeighborCache.clear();

    // Watershed until a threshold is reached.
    ind dataIdx = Data->Indices[CurrentWatershedIndex];
    while (Data->Scalars[dataIdx] <= maxVal) {
        // Get cluster ID and representative vertex for each neighbor.
        vec3i globIdx = Data->BlockOffset + vec3i::fromIndexOfTotal(dataIdx, Data->BlockSize);
        for (int dim = 0; dim < 3; ++dim)
            for (int sign = -1; sign <= 1; sign += 2) {
                vec3i neighIdx = globIdx;
                neighIdx[dim] += sign;
                if (neighIdx[dim] < 0 || neighIdx[dim] > Data->BlockSize[dim]) continue;

                vec3i neighRep;
                ClusterID* neighCluster = findClusterID(neighIdx, neighRep);

                //
                if (neighCluster) {
                    bool addCluster = true;
                    //
                    for (Neighbor& neigh : NeighborCache)
                        if (neigh.Cluster == *neighCluster) {
                            addCluster = false;
                            break;
                        }
                    if (addCluster) NeighborCache.emplace_back(*neighCluster, neighRep);
                }
            }

        *(PointerBlock.getPointer(globIdx)) =
            NeighborProcessor.doWatershed(globIdx, Data->Volumes[dataIdx], NeighborCache);

        // Watershed on.
        ind dataIdx = Data->Indices[++CurrentWatershedIndex];
    }  // namespace perc
}

/* Return final cluster ID (or nullptr if none yet).
 * Last vertex pointing directly to cluster ID is returned as second argument.
 */
template <typename ClusterProcessor>
ClusterID* UnionFindSubBlock<ClusterProcessor>::findClusterID(const vec3i& idx,
                                                              vec3i& lastVertexID) {
    vec3i curIdx = idx;
    assert(contains(idx) && "First index is expected to be in here.");
    ID* firstPointer = PointerBlock.getPointer(idx);
    ID* curPointer = firstPointer;
    vec3i lastPointer = idx;

    ClusterID* finalID;
    while (curPointer && curPointer->isVertex()) {
        curIdx = vec3i::fromIndexOfTotal(curPointer->RawID, PointerBlock.TotalSize);

        if (!PointerBlock.contains(curIdx)) break;
        lastPointer = curIdx;

        ID* curPointer = PointerBlock.getPointer(curIdx);
        if (curPointer->isCluster()) {
            finalID = curPointer->asCluster();
            break;
        }
    }
    // No cluster yet.
    if (!curPointer) return nullptr;
    // Not there yet, give it back to parent block.
    if (curPointer->isVertex()) finalID = Parent.findClusterID(idx, lastPointer);

    assert(finalID && finalID->baseID() >= 0 && "Invalid result of union find.");
    *firstPointer = lastPointer.toIndexOfTotal(PointerBlock.TotalSize);

    return finalID;
}

}  // namespace perc