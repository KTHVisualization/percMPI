#pragma once
#include <vector>
#include "vec.h"
#include "handles.h"
#include "clusterlist.h"
#include "clusterlistrecording.h"
#include "whiteprocessor.h"
#include "redprocessor.h"
#include "grayprocessor.h"
#include "unionfindsubblock.h"

namespace perc {

class LocalBlock : public UnionFindBlock {
public:
    LocalBlock(const vec3i& blockSize, const vec3i& blockOffset, const vec3i& totalSize,
               const int rank = 0);
    LocalBlock(const LocalBlock&) = delete;
    LocalBlock operator=(const LocalBlock&) = delete;
    LocalBlock(LocalBlock&& other);

private:
    LocalBlock(const vec3i& totalSize);

public:
    static LocalBlock* makeGroundtruth(const vec3i& blockSize, const vec3i& blockOffset,
                                       const vec3i& totalSize);
    static LocalBlock* makeWhiteRedTest(const vec3i& blockSize, const vec3i& blockOffset,
                                        const vec3i& totalSize);
    static LocalBlock* makeWhiteRedGreenTest(const vec3i& blockSize, const vec3i& blockOffset,
                                             const vec3i& totalSize);
    // 0: All, 1: Only White, 2: Only Red, 3: Only Green
    void outputFrontBlocks(std::vector<char>& field, ind slice, ind selection);

    ~LocalBlock();

    virtual void doWatershed(const double minVal) override;
    virtual ClusterID* findClusterID(const vec3i& idx, vec3i& lastClusterID) override;
    virtual ID* setID(const vec3i& idx, const ID& id) override;
    using UnionFindBlock::setID;
    virtual double getClusterVolume(ClusterID cluster) override;

    virtual void receiveData() override;
    virtual void sendData() override;
    virtual ind numClusters() override { return LOLs->numClusters(); }
    virtual ind numClustersCombined() override {
        assert(LOLs->numClusters() + LOGs->numClusters() >= 0 &&
               "Combined number of clusters smaller then 0. Not initialized?");
        return LOLs->numClusters() + LOGs->numClusters();
    }
    virtual double totalVolume() override { return LOLs->totalVolume(); }
    virtual double totalVolumeCombined() override {
        return LOLs->totalVolume() + LOGs->totalVolume();
    }
    virtual double maxVolume() override { return LOLs->maxVolume(); }
    virtual double maxVolumeCombined() override {
        return std::max(LOLs->maxVolume(), LOGs->maxVolume());
    }
    virtual ind memEstimate() const override;

    void checkConsistency() const override;
    std::vector<std::pair<vec3i, double>> getVoluminaForAddedVertices(double maxVal) override;

protected:
    void repointerMultipleMerges(const std::vector<ind>& connComps);

private:
    // Data send to the global node
    struct DataToSend {
        ind NumClusters;
        double MaxVolume;
        double TotalVolume;
        // List of PLOGS
        std::vector<ClusterData> PLOGs;
        ind VectorSizes[2];
    } CommData;

    // This blocks rank (if on a single node)
    int Rank;

    // The local local part for this block
    UnionFindSubBlock<WhiteProcessor>* LOLSubBlock;

    // Memory for LOG ID blocks.
    ID* MemoryLOG;
    // Size of that memory
    ind MemoryLOGSize;

    // The local global part for this block
    std::vector<UnionFindSubBlock<RedProcessor>> LOGSubBlocks;

    // The global part for this block, this is only for lookup for the local node
    std::vector<UnionFindSubBlock<GrayProcessor>> GOGSubBlocks;

    // Local representations of local clusters.
    ClusterList* LOLs;
    // Local representations of global clusters.
    ClusterListRecordingSingle* LOGs;

    using RefPLOGtype = std::unordered_set<ClusterID, ClusterID::hash_type>;
    // Potential LOGs: LOLs that touch the boundary.
    RefPLOGtype* RefPLOGs;
};

}  // namespace perc