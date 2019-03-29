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
    struct InfoPerProcess;

    GlobalBlock(const vec3i& blockSize, const vec3i& totalSize, const vec3i& numNodes);
    GlobalBlock(const GlobalBlock&) = delete;
    GlobalBlock operator=(const GlobalBlock&) = delete;

    ~GlobalBlock() {
        for (auto& processData : PerProcessData) {
            delete processData.MemoryLOG;
        }
    }

private:
    GlobalBlock(const vec3i& totalSize);

public:
    static GlobalBlock* makeGreenTest(const vec3i& blockSize, const vec3i& blockOffset,
                                      const vec3i& totalSize);
    static GlobalBlock* makeWhiteRedTest(const vec3i& blockSize, const vec3i& blockOffset,
                                         const vec3i& totalSize);
    static GlobalBlock* makeWhiteRedGreenTest(const vec3i& blockSize, const vec3i& blockOffset,
                                              const vec3i& totalSize);

    virtual void doWatershed(const double minVal) override;
    virtual ClusterID* findClusterID(const vec3i& idx, vec3i& lastClusterID) override;

    virtual ID* setID(const vec3i& idx, const ID& id) override;
    using UnionFindBlock::setID;

    virtual double getClusterVolume(ClusterID cluster) override;

    virtual void reset() override;
    virtual void loadData() override;
    virtual void sortData(bool useBuckets) override;

    virtual void receiveData() override;
    virtual void sendData() override;
    virtual ind numClusters() override { return GOGs.numClusters(); };
    virtual ind numClustersCombined() override { return GOGs.numClusters() + NumClustersLocal; }
    virtual double totalVolume() override { return GOGs.totalVolume(); };
    virtual double totalVolumeCombined() override { return GOGs.totalVolume() + TotalVolumeLocal; }
    virtual double maxVolume() override { return GOGs.maxVolume(); };
    virtual double maxVolumeCombined() override {
        return std::max(GOGs.maxVolume(), MaxVolumeLocal);
    };
    ind greenSize() const { return MemoryGreenSize; };
    ind redSize() const { return MemoryRedSize; };

    virtual ind memEstimate() const override;

    virtual void checkConsistency() const override;
    virtual std::vector<std::pair<vec3i, double>> getVoluminaForAddedVertices(
        double maxVal) override;

protected:
    void repointerMultipleMerges(const std::vector<ind>& connComps);

public:
    struct InfoPerProcess {
        InfoPerProcess() : MemoryLOG(nullptr), MemoryLOGSize(0), StartOfLocalPlog(-1) {}

        InfoPerProcess(const std::vector<ind>& greenAdjacent, ID* memoryLOG, ind memoryLOGSize,
                       std::vector<ClusterMerge> merges = std::vector<ClusterMerge>())
            : GreenAdjacent(greenAdjacent)
            , MemoryLOG(memoryLOG)
            , MemoryLOGSize(memoryLOGSize)
            , Merges(merges)
            , StartOfLocalPlog(ind(0)) {}

        std::vector<ind> GreenAdjacent;
        ind LOGstartIndex;
        ID* MemoryLOG;
        ind MemoryLOGSize;
        std::vector<ClusterMerge> Merges;
        ind StartOfLocalPlog;
    };

private:
    // The global (actual) part for this block
    std::vector<UnionFindSubBlock<GreenProcessor>> GOGSubBlocks;
    // The local global part for this block, just for lookup
    std::vector<UnionFindSubBlock<GrayProcessor>> LOGSubBlocks;

    // Global representations of global clusters.
    ClusterListRecordingMultiple GOGs;

    // Needed for efficient lookup
    vec3i BlockSize;

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
    std::vector<std::vector<ClusterMerge>*> ReceivedMerges;
    // *** Mixed ***
    std::vector<InfoPerProcess> PerProcessData;
    // Green pointer blocks data and its size
    ID* MemoryGreen;
    ind MemoryGreenSize;
    ind MemoryRedSize;
};

}  // namespace perc