#pragma once
#include "handles.h"
#include "clusterlist.h"
#include "clusterlistrecording.h"

namespace perc {

template <typename ClusterProcessor>
class UnionFindSubBlock;

struct LocalLocalProcessor {

    LocalLocalProcessor(ClusterList& lols, ClusterListRecording& logs,
                        std::vector<ClusterID>& plogs)
        : LOLs(lols), LOGs(logs), PLOGs(plogs), Parent(nullptr) {}

    ID doWatershed(VertexID pos, double volume, std::vector<Neighbor>& neighClusters);

    void setParent(UnionFindSubBlock<LocalLocalProcessor>* parent);

private:
    UnionFindSubBlock<LocalLocalProcessor>* Parent;
    ClusterList& LOLs;
    ClusterListRecording& LOGs;
    std::vector<ClusterID>& PLOGs;
};

}  // namespace perc