#pragma once
#include <unordered_set>
#include <vector>
#include "vec.h"
#include "datablock.h"
#include "idblock.h"

namespace perc {

class UnionFindBlock;

template <typename ClusterProcessor>
class UnionFindSubBlock {
public:
    friend ClusterProcessor;

    UnionFindSubBlock<ClusterProcessor>(const vec3i& size, const vec3i& offset, const vec3i& total,
                                        UnionFindBlock& parent, ClusterProcessor&& neighProcessor,
                                        ID* memory = nullptr);

    ~UnionFindSubBlock() { delete Data; }

    void loadData();
    void doWatershed(const double minVal);

    ClusterID* findClusterID(const vec3i& idx, vec3i& lastClusterID);

    bool contains(const vec3i& idx) { return PointerBlock.contains(idx); }

    const vec3i& blockSize() { return PointerBlock.BlockSize; }
    const vec3i& blockOffset() { return PointerBlock.BlockOffset; }
    const vec3i& totalSize() { return PointerBlock.TotalSize; }

    void checkConsistency() const;
    using VolumeStat = std::pair<vec3i, double>;
    void getVoluminaForAddedVertices(double maxVal, std::vector<VolumeStat>& stats);

public:
    DataBlock* Data;
    IDBlock PointerBlock;
    UnionFindBlock& Parent;
    // Until which index we have watershedded in the sorted list.
    ind CurrentWatershedIndex;

private:
    mutable std::vector<Neighbor> NeighborCache;
    ClusterProcessor NeighborProcessor;
};  // namespace perc

}  // namespace perc

#include "unionfindsubblock.inl"