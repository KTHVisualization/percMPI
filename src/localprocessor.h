#pragma once
#include "handles.h"
#include "clusterlist.h"
#include "clusterlistrecording.h"

namespace perc {

template <typename ClusterProcessor>
class UnionFindSubBlock;

struct LocalLocalProcessor {

    LocalLocalProcessor(ClusterList& lols, ClusterListRecording& logs)
        : LOLs(lols), LOGs(logs), Parent(nullptr) {}

    ID doWatershed(VertexID pos, double volume, std::vector<Neighbor>& neighClusters);

    void setParent(UnionFindSubBlock<LocalLocalProcessor>* parent);

private:
    UnionFindSubBlock<LocalLocalProcessor>* Parent;
    ClusterList& LOLs;
    ClusterListRecording& LOGs;
};

}  // namespace perc