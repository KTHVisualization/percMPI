#pragma once
#include <vector>
#include "vec.h"
#include "handles.h"
#include "clusterlist.h"
#include "clusterlistrecording.h"
#include "globalprocessor.h"
#include "localprocessor.h"
#include "unionfindsubblock.h"

namespace perc {

class GreenBlock : public UnionFindBlock {
public:
    GreenBlock(const vec3i& blockSize, const vec3i& blockOffset, const vec3i& totalSize);

    virtual void doWatershed(const double minVal) override;
    virtual ClusterID* findClusterID(const vec3i& idx, vec3i& lastClusterID) override;

    virtual ID* setID(const vec3i& idx, const ID& id) override;
    using UnionFindBlock::setID;

    virtual double getClusterVolume(ClusterID cluster) override;

    // Sketch.
    virtual void receiveData() override;
    virtual void sendData() override;
    virtual ind numClusters() override { return GOGs.numClusters(); };
    virtual double totalVolume() override { return GOGs.totalVolume(); };
    virtual double maxVolume() override { return GOGs.maxVolume(); };

    virtual void checkConsistency() const override;

private:
    // The global (actual) part for this block
    std::vector<UnionFindSubBlock<GlobalProcessor>*> GOGSubBlocks;
    // The local global part for this block, just for lookup
    std::vector<UnionFindSubBlock<LocalGlobalProcessor>*> LOGSubBlocks;
    // Local representations of local clusters.
    ClusterListRecordingMultiple GOGs;
};  // namespace perc

}  // namespace perc