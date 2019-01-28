#pragma once
#include "unionfindsubblock.h"
#include "unionfindblock.h"
#include <iostream>

namespace perc {

template <typename ClusterProcessor>
UnionFindSubBlock<ClusterProcessor>::UnionFindSubBlock(const vec3i& size, const vec3i& offset,
                                                       const vec3i& total, UnionFindBlock& parent,
                                                       ClusterProcessor&& neighProcessor,
                                                       ID* memory)
    : PointerBlock(size, offset, total, memory)
    , Parent(&parent)
    , CurrentWatershedIndex(-1)
    , NeighborProcessor(neighProcessor)
    , Data(nullptr) {

    assert(vec3i::min(offset + size, total) == offset + size &&
           "Given block would lay over the maximum.");
    assert(vec3i::max(offset + size, vec3i(0)) == offset + size &&
           "Given block would lay below the minimum.");
    NeighborCache.reserve(6);
    NeighborProcessor.setParent(this);
    std::cout << " block\t" << offset << "\t - " << offset + size << std::endl;
}

template <typename ClusterProcessor>
UnionFindSubBlock<ClusterProcessor>::UnionFindSubBlock(UnionFindSubBlock<ClusterProcessor>&& other)
    : PointerBlock(std::move(other.PointerBlock))
    , Parent(nullptr)
    , CurrentWatershedIndex(other.CurrentWatershedIndex)
    , NeighborProcessor(std::move(other.NeighborProcessor)) {
    Data = other.Data;
    other.Data = nullptr;
    // We copied the parent pointer over from other, thus not calling 'setParent' here.
    NeighborProcessor.setParent(this);
}

// template <typename ClusterProcessor>
// void UnionFindSubBlock<ClusterProcessor>::setParent(UnionFindBlock& parent) {
//     assert(!Parent && "Parent already set!");
//     Parent = parent;
// }

template <typename ClusterProcessor>
void UnionFindSubBlock<ClusterProcessor>::loadData() {
    assert(!Data && "Data was already set.");
    Data = new DataBlock(PointerBlock.BlockSize, PointerBlock.BlockOffset, PointerBlock.TotalSize);
    // TODO: Hardcoded time slice and data directory for now
    Data->loadData(1, "../../Data/Percolation/P3", "uv_000");
    CurrentWatershedIndex = 0;
    Data->sort();
}

template <typename ClusterProcessor>
void UnionFindSubBlock<ClusterProcessor>::doWatershed(const double minVal) {
    if (!Data) return;

    // Watershed until a threshold is reached.
    ind dataIdx = Data->Indices[CurrentWatershedIndex];
    ind finalWaterShedIndex = Data->BlockSize.prod();
    while (CurrentWatershedIndex < finalWaterShedIndex && Data->Scalars[dataIdx] >= minVal) {
        // Get cluster ID and representative vertex for each neighbor.
        vec3i globIdx = Data->BlockOffset + vec3i::fromIndexOfTotal(dataIdx, Data->BlockSize);
        NeighborCache.clear();
        for (int dim = 0; dim < 3; ++dim)
            for (int sign = -1; sign <= 1; sign += 2) {
                vec3i neighIdx = globIdx;
                neighIdx[dim] += sign;
                if (neighIdx[dim] < 0 || neighIdx[dim] >= Data->TotalSize[dim]) continue;

                vec3i neighRep;
                ClusterID* neighCluster;
                if (contains(neighIdx))
                    neighCluster = findClusterID(neighIdx, neighRep);
                else
                    neighCluster = Parent->findClusterID(neighIdx, neighRep);

                if (neighCluster) {
                    bool addCluster = true;
                    for (Neighbor& neigh : NeighborCache)
                        if (neigh.Cluster == *neighCluster) {
                            addCluster = false;
                            break;
                        }
                    if (addCluster) NeighborCache.emplace_back(*neighCluster, neighRep);
                }
            }

        ID newIdx = NeighborProcessor.doWatershed(globIdx.toIndexOfTotal(Data->TotalSize),
                                                  Data->Volumes[dataIdx], NeighborCache);

        PointerBlock.setPointer(globIdx, newIdx);

        // Watershed on.
        dataIdx = Data->Indices[++CurrentWatershedIndex];
    }
}

/* Return final cluster ID (or nullptr if none yet).
 * Last vertex pointing directly to cluster ID is returned as second argument.
 */
template <typename ClusterProcessor>
ClusterID* UnionFindSubBlock<ClusterProcessor>::findClusterID(const vec3i& idx,
                                                              vec3i& lastVertexID) {
    vec3i curIdx = idx;
    assert(contains(curIdx) && "First index is expected to be in here.");
    ID* firstPointer = PointerBlock.getPointer(idx);
    ID* curPointer = firstPointer;
    lastVertexID = idx;

    while (curPointer && curPointer->isVertex()) {
        curIdx = vec3i::fromIndexOfTotal(curPointer->RawID, PointerBlock.TotalSize);

        if (!PointerBlock.contains(curIdx)) break;
        lastVertexID = curIdx;

        curPointer = PointerBlock.getPointer(curIdx);
        assert(curPointer && "Should not point nowhere.");
        if (curPointer->isCluster()) {
            break;
        }
    }
    // No cluster yet.
    if (!curPointer) return nullptr;

    ClusterID* finalID;
    // Not there yet, give it back to parent block.
    if (curPointer->isVertex()) {
        finalID = Parent->findClusterID(curIdx, lastVertexID);
        assert(finalID && "Vertex is pointed to but not pointing to anything.");
    }
    // It is a cluster, set final ID
    else {
        finalID = curPointer->asCluster();
    }

    assert(finalID && finalID->baseID() >= 0 && "Invalid result of union find.");

    // Path compression, if path compression to be done
    if (firstPointer != curPointer) {
        *firstPointer = lastVertexID.toIndexOfTotal(PointerBlock.TotalSize);
    }

    return finalID;
}

template <typename ClusterProcessor>
void UnionFindSubBlock<ClusterProcessor>::checkConsistency() const {
#ifndef NDEBUG
    NeighborProcessor.checkConsistency();
#endif
}

template <typename ClusterProcessor>
void UnionFindSubBlock<ClusterProcessor>::getVoluminaForAddedVertices(
    double maxVal, std::vector<VolumeStat>& stats) {
    vec3i dummy(-1, -1, -1);
    ind currID = CurrentWatershedIndex - 1;
    if (currID < 0) return;
    ind dataIdx = Data->Indices[currID];

    while (currID >= 0 && Data->Scalars[dataIdx] < maxVal) {
        vec3i globIdx = Data->BlockOffset + vec3i::fromIndexOfTotal(dataIdx, Data->BlockSize);

        // Find cluster we belong to and get volume.
        ClusterID& cluster = *findClusterID(globIdx, dummy);
        double volume = Parent->getClusterVolume(cluster);
        stats.push_back({globIdx, volume});

        currID--;
        dataIdx = Data->Indices[currID];
    }
}

}  // namespace perc