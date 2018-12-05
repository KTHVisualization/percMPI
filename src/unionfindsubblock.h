#pragma once
#include <unordered_set>
#include <vector>
#include "vec.h"
#include "datablock.h"
#include "unionfind.h"

namespace perc {

class UnionFindBlock;

template <typename ClusterProcessor>
class UnionFindSubBlock {
public:
    friend ClusterProcessor;

    UnionFindSubBlock<ClusterProcessor>();

    UnionFindSubBlock<ClusterProcessor>(const vec3i& size, const vec3i& offset, const vec3i& total,
                                        UnionFindBlock& parent, ClusterProcessor&& neighProcessor);

    ~UnionFindSubBlock() { delete Data; }
    void loadData();
    void doWatershed(double minVal);

    ClusterID* findClusterID(const vec3i& idx, vec3i& lastClusterID);

    bool contains(const vec3i& idx) { return PointerBlock.contains(idx); }

    const vec3i& blockSize() { return Data->BlockSize; }
    const vec3i& blockOffset() { return Data->BlockOffset; }
    const vec3i& totalSize() { return Data->TotalSize; }

public:
    DataBlock* Data;
    UnionFind PointerBlock;
    UnionFindBlock& Parent;
    // Until which index we have watershedded in the sorted list.
    ind CurrentWatershedIndex;

private:
    mutable std::vector<Neighbor> NeighborCache;
    ClusterProcessor NeighborProcessor;
};

}  // namespace perc

#include "unionfindsubblock.inl"