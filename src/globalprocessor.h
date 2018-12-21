#pragma once
#include "handles.h"
#include "clusterlistmultiple.h"
#include "clusterlistrecording.h"

namespace perc {

template <typename ClusterProcessor>
class UnionFindSubBlock;

struct GlobalProcessor {

    GlobalProcessor(ClusterListRecordingMultiple& gogs) : GOGs(gogs), Parent(nullptr) {}

    ID doWatershed(VertexID pos, double volume, std::vector<Neighbor>& neighClusters);

    void setParent(UnionFindSubBlock<GlobalProcessor>* parent);

private:
    UnionFindSubBlock<GlobalProcessor>* Parent;
    ClusterListRecordingMultiple& GOGs;
};

}  // namespace perc