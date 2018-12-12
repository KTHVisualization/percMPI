#pragma once
#include <vector>
#include "vec.h"
#include "handles.h"
#include "clusterlist.h"
#include "clusterlistrecording.h"
#include "globalprocessor.h"
#include "unionfindsubblock.h"

namespace perc {

class GreenBlock : public UnionFindBlock {
public:
    GreenBlock(const vec3i& blockSize, const vec3i& blockOffset, const vec3i& totalSize);

    virtual void doWatershed(const double minVal) override;
    virtual ClusterID* findClusterID(const vec3i& idx, vec3i& lastClusterID) override;

    // Sketch.
    virtual void receiveData() override;
    virtual void sendData() override;
    virtual ind numClusters() override;
    virtual double totalVolume() override;
    virtual double maxVolume() override;

private:
    // The local local part for this block
    // std::vector<UnionFindSubBlock<GlobalGlobalProcessor>*> GOGSubBlocks;
    // TODO: The local global part for this block
    // std::vector<UnionFindSubBlock<LocalGlobalProcessor>*> LOGSubBlocks;
    // TODO: The global part for this block, this is only for lookup for the local node
    // std::vector<UnionFindSubBlock<GlobalProcessor>*> GOGSubBlocks;
    // Local representations of local clusters.
    // ClusterListMultiple LOLs;

    std::vector<ClusterMerge> Merges;
};

}  // namespace perc