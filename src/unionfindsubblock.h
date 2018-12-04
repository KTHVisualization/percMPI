#pragma once
#include <unordered_set>
#include "vec.h"
#include "datablock.h"
#include "unionfind.h"
#include "unionfindblock.h"

namespace perc {

class UnionFindSubBlock {
public:
    UnionFindSubBlock(const vec3i& size, const vec3i& offset, const vec3i& total,
                      UnionFindBlock& parent);

    ~UnionFindSubBlock() { delete Data; }
    void loadData();
    void doWatershed(double maxVal);

    ClusterID* findClusterID(const vec3i& idx, vec3i& lastClusterID);

    bool contains(const vec3i& idx) { return PointerBlock.contains(idx); }

public:
    DataBlock* Data;
    UnionFind PointerBlock;
    UnionFindBlock& Parent;
    ind CurrentWatershedIndex;

    struct Neighbor {
        Neighbor(ClusterID cluster, VertexID representative)
            : Cluster(cluster), Representative(representative) {}
        ClusterID Cluster;
        VertexID Representative;
    };

private:
    mutable std::vector<Neighbor> NeighborCache;
    //    mutable std::unordered_set<ClusterID> NeighborClusterCache;
};

}  // namespace perc