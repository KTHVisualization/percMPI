#pragma once
#include <vector>
#include "vec.h"
#include "handles.h"
#include "clusterlist.h"
#include "clusterlistrecording.h"
#include "unionfindsubblock.h"
#include "greenprocessor.h"
#include "grayprocessor.h"

namespace perc {

class GlobalBlock : public UnionFindBlock {
public:
    struct InfoPerProcess {
        InfoPerProcess(std::vector<ind> greenIndices, ID* memoryLOG, ind memoryLOGSize,
                       std::vector<ClusterMerge>* merges)
            : GreenIndices(greenIndices)
            , MemoryLOG(memoryLOG)
            , MemoryLOGSize(memoryLOGSize)
            , Merges(merges)
            , StartOfLocalPlog(0) {}
        InfoPerProcess()
            : GreenIndices()
            , MemoryLOG(nullptr)
            , MemoryLOGSize(0)
            , Merges(nullptr)
            , StartOfLocalPlog(0) {}
        std::vector<ind> GreenIndices;
        ID* MemoryLOG;
        ind MemoryLOGSize;
        const std::vector<ClusterMerge>* Merges;
        ind StartOfLocalPlog;
    };

    GlobalBlock(const vec3i& blockSize, const vec3i& blockOffset, const vec3i& totalSize,
                const vec3i& numNodes);

    virtual void doWatershed(const double minVal) override;
    virtual ClusterID* findClusterID(const vec3i& idx, vec3i& lastClusterID) override;

    virtual ID* setID(const vec3i& idx, const ID& id) override;
    using UnionFindBlock::setID;

    virtual double getClusterVolume(ClusterID cluster) override;

    virtual void receiveData() override;
    virtual void sendData() override;
    virtual ind numClusters() override { return GOGs.numClusters(); };
    ind totalNumClusters() { return GOGs.numClusters() + NumClustersLocal; }
    virtual double totalVolume() override { return GOGs.totalVolume(); };
    double totalTotalVolume() { return GOGs.totalVolume() + TotalVolumeLocal; }
    virtual double maxVolume() override { return GOGs.maxVolume(); };
    double totalMaxVolume() { return std::max(GOGs.maxVolume(), MaxVolumeLocal); };

    virtual void checkConsistency() const override;
    virtual std::vector<std::pair<vec3i, double>> getVoluminaForAddedVertices(
        double maxVal) override;

protected:
    void repointerMultipleMerges(const std::vector<ind>& connComps);

private:
    // The global (actual) part for this block
    std::vector<UnionFindSubBlock<GreenProcessor>> GOGSubBlocks;
    // The local global part for this block, just for lookup
    std::vector<UnionFindSubBlock<GrayProcessor>> LOGSubBlocks;

    // Global representations of global clusters.
    ClusterListRecordingMultiple GOGs;

    // ****** Structures for sending and receiving ******
    vec3i NumNodes;
    // *** Sending ***
    // All global merges that happen in one step
    std::vector<ind> Merges;
    // Number of (global) clusters created in red and green
    ind NumNewClusters;
    // *** Receiving ****
    ind NumClustersLocal;
    double MaxVolumeLocal;
    double TotalVolumeLocal;
    // Merges received
    std::vector<std::vector<ClusterMerge>> ReceivedMerges;
    // *** Mixed ***
    std::vector<InfoPerProcess> PerProcessData;
};

}  // namespace perc