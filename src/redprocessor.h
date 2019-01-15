#pragma once
#include <unordered_set>
#include "handles.h"
#include "clusterlist.h"
#include "clusterlistrecording.h"

namespace perc {

template <typename ClusterProcessor>
class UnionFindSubBlock;

struct RedProcessor {

    RedProcessor(ClusterList& lols, ClusterListRecording<ClusterList>& logs,
                 std::unordered_set<ClusterID, ClusterID::hash_type>& plogs)
        : LOLs(lols), LOGs(logs), PLOGs(plogs), Parent(nullptr) {}

    ID doWatershed(VertexID pos, double volume, std::vector<Neighbor>& neighClusters);

    void setParent(UnionFindSubBlock<RedProcessor>* parent);

    void checkConsistency() const;

private:
    UnionFindSubBlock<RedProcessor>* Parent;
    ClusterList& LOLs;
    ClusterListRecording<ClusterList>& LOGs;
    std::unordered_set<ClusterID, ClusterID::hash_type>& PLOGs;
};

}  // namespace perc