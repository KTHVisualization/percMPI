#pragma once
#include <vector>
#include "vec.h"
#include "handles.h"
#include "clusterlist.h"
#include "clusterlistrecording.h"
#include "localprocessor.h"
#include "unionfindsubblock.h"

namespace perc {

class WhiteBlock : public UnionFindBlock {
public:
    WhiteBlock(const vec3i& blockSize, const vec3i& blockOffset, const vec3i& totalSize);
    ~WhiteBlock() { delete[] MemoryLOG; }

    virtual void doWatershed(const double minVal) override;
    virtual ClusterID* findClusterID(const vec3i& idx, vec3i& lastClusterID) override;
    virtual ID* setID(const vec3i& idx, const ID& id) override;
    using UnionFindBlock::setID;

    // Sketch.
    virtual void receiveData() override;
    virtual void sendData() override;
    virtual ind numClusters() override { return LOLs.numClusters(); }
    ind totalNumClusters() {
        if (LOLs.numClusters() + LOGs.numClusters() < 0)
            std::cout << "WTF? " << LOLs.numClusters() << " + " << LOGs.numClusters() << std::endl;
        return LOLs.numClusters() + LOGs.numClusters();
    }
    virtual double totalVolume() override { return LOLs.totalVolume(); }
    double totalTotalVolume() { return LOLs.totalVolume() + LOGs.totalVolume(); }
    virtual double maxVolume() override { return LOLs.maxVolume(); }
    double totalMaxVolume() { return std::max(LOLs.maxVolume(), LOGs.maxVolume()); }

    void checkConsistency() const;


private:
    // The local local part for this block
    UnionFindSubBlock<LocalLocalProcessor>* LOLSubBlock;

    // Memory for LOG ID blocks.
    ID* MemoryLOG;
    // The local global part for this block
    std::vector<UnionFindSubBlock<LocalGlobalProcessor>> LOGSubBlocks;

    // TODO: The global part for this block, this is only for lookup for the local node
    // std::vector<UnionFindSubBlock<GlobalProcessor>> GOGSubBlocks;

    // Local representations of local clusters.
    ClusterList LOLs;
    // Local representations of global clusters.
    ClusterListRecordingSingle LOGs;
    // Potential LOGs: LOLs that touch the boundary.
    std::unordered_set<ClusterID, ClusterID::hash_type> RefPLOGs;
    // List of PLOGS that will be used for sending and receiving data.
    std::vector<ClusterData> CommPLOGs;
};

}  // namespace perc