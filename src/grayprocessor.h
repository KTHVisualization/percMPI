#pragma once
#include "handles.h"
#include <vector>

namespace perc {

template <typename ClusterProcessor>
class UnionFindSubBlock;

// The red block seen from the global node,
// the green block seen from the local node.
struct GrayProcessor {

    GrayProcessor() { std::cout << "Gray"; }

    ID doWatershed(VertexID pos, double volume, std::vector<Neighbor>& neighClusters) {
        assert(false && "Read only.");
    }

    void setParent(UnionFindSubBlock<GrayProcessor>* parent) {}

    void checkConsistency() const {}
};
}  // namespace perc