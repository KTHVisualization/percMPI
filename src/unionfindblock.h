#pragma once
#include <vector>
#include "vec.h"
#include "handles.h"
#include "clusterlist.h"
#include "clusterlistrecording.h"

namespace perc {

class UnionFindBlock {
public:
    virtual void doWatershed(double maxVal) = 0;
    virtual ClusterID* findClusterID(const vec3i& idx, vec3i& lastClusterID) = 0;
    virtual void receiveData() = 0;
    virtual void sendData() = 0;
};

class LocalBlock : public UnionFindBlock {
public:
    LocalBlock(const vec3i& blockSize, const vec3i& blockOffset, const vec3i& totalSize);

    virtual void doWatershed(double maxVal) override;
    virtual ClusterID* findClusterID(const vec3i& idx, vec3i& lastClusterID) override;

    // Sketch.
    virtual void receiveData() override;
    virtual void sendData() override;

private:
    // Local representations of local clusters.
    ClusterList LOLs;
    // Local representations of global clusters.
    ClusterListRecording LOGs;
    // Potential LOGs: LOLs that touch the boundary.
    std::vector<ClusterID> RefPLOGs;
    // List of PLOGS that will be used for sending and receiving data.
    std::vector<Cluster> CommPLOGs;
};  // namespace perc

}  // namespace perc