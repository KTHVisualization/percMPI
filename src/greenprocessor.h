#pragma once
#include "handles.h"
#include "clusterlistmultiple.h"
#include "clusterlistrecording.h"

namespace perc {

template <typename ClusterProcessor>
class UnionFindSubBlock;

struct GreenProcessor {

    GreenProcessor(ClusterListRecordingMultiple& gogs) : GOGs(gogs), Parent(nullptr) {
        std::cout << "Green";
    }

    ID doWatershed(VertexID pos, double volume, std::vector<Neighbor>& neighClusters);

    void setParent(UnionFindSubBlock<GreenProcessor>* parent);

    void checkConsistency() const;

private:
    UnionFindSubBlock<GreenProcessor>* Parent;
    ClusterListRecordingMultiple& GOGs;
};

}  // namespace perc