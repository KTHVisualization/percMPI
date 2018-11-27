#pragma once
#include <vector>
#include <stack>
#include "vec.h"
#include "handles.h"

namespace perc {

struct Cluster {
    Cluster(VertexID idx) : Index(idx), Volume(0) {}
    VertexID Index;
    double Volume;
};

class ClusterList {
    Cluster& getCluster(ClusterID cluster) { return Clusters[cluster.localID()]; }
    inline ClusterID addCluster(VertexID id);
    inline void removeCluster(ClusterID cluster);

private:
    std::vector<Cluster> Clusters;
    std::stack<size_t, std::vector<size_t>> Holes;
};

// ========= Inline Definitions ========= //

ClusterID ClusterList::addCluster(VertexID id) {
    if (Holes.empty()) {
        Clusters.emplace_back(id);
        return ClusterID(Clusters.size() - 1);
    } else {
        size_t holeIdx = Holes.top();
        Holes.pop();
        Clusters[holeIdx] = Cluster(id);
        return ClusterID(holeIdx);
    }
}
void ClusterList::removeCluster(ClusterID cluster) {
    ind locID = cluster.localID();
    if (locID == Clusters.size() - 1)
        Clusters.pop_back();
    else {
        Clusters[locID] = Cluster(-1);
        Holes.push(locID);
    }
}

}  // namespace perc