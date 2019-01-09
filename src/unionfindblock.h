#pragma once
#include <vector>
#include "vec.h"
#include "handles.h"

namespace perc {

class UnionFindBlock {
public:
    UnionFindBlock(const vec3i& totalSize) : TotalSize(totalSize) {}
    virtual void doWatershed(const double minVal) = 0;
    virtual ClusterID* findClusterID(const vec3i& idx, vec3i& lastClusterID) = 0;
    virtual ID* setID(const vec3i& idx, const ID& id) = 0;
    ID* setID(VertexID idx, const VertexID& id) {
        return setID(vec3i::fromIndexOfTotal(idx.RawID, TotalSize), id);
    }
    virtual double getClusterVolume(ClusterID cluster) = 0;
    virtual void receiveData() = 0;
    virtual void sendData() = 0;
    virtual ind numClusters() = 0;
    virtual double totalVolume() = 0;
    virtual double maxVolume() = 0;

    virtual void checkConsistency() const = 0;

    const vec3i TotalSize;
};

}  // namespace perc