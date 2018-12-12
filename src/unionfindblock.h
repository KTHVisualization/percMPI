#pragma once
#include <vector>
#include "vec.h"
#include "handles.h"

namespace perc {

class UnionFindBlock {
public:
    virtual void doWatershed(const double minVal) = 0;
    virtual ClusterID* findClusterID(const vec3i& idx, vec3i& lastClusterID) = 0;
    virtual void receiveData() = 0;
    virtual void sendData() = 0;
    virtual ind numClusters() = 0;
    virtual double totalVolume() = 0;
    virtual double maxVolume() = 0;
};

}  // namespace perc