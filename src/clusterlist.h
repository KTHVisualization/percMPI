#pragma once
#include <vector>
#include <stack>
#include "vec.h"
#include "handles.h"

namespace perc {

struct Cluster {
    Cluster(VertexID& idx, double& vol) : Index(idx), Volume(vol) {}
    VertexID& Index;
    double& Volume;
};

class ClusterList {
public:
    ClusterList(ind size = 100) : Indices(size), Volumes(size) {}

    Cluster getCluster(ClusterID cluster) {
        return Cluster(Indices[cluster.localID()], Volumes[cluster.localID()]);
    }
    inline ClusterID addCluster(VertexID id);
    inline void removeCluster(ClusterID cluster);
    inline void mergeClusters(ClusterID from, ClusterID onto);

    void clearVolumes() { std::fill(Volumes.data(), Volumes.data() + Volumes.size(), 0); }

private:
    std::vector<VertexID> Indices;
    std::vector<double> Volumes;
    std::vector<size_t> Holes;
};

// ========= Inline Definitions ========= //

ClusterID ClusterList::addCluster(VertexID id) {
    if (Holes.empty()) {
        Indices.push_back(id);
        Volumes.push_back(0.0);
        return ClusterID(Indices.size() - 1);
    } else {
        size_t holeIdx = Holes.back();
        Holes.pop_back();
        Indices[holeIdx] = id;
        Volumes[holeIdx] = 0.0;
        return ClusterID(holeIdx);
    }
}
void ClusterList::removeCluster(ClusterID cluster) {
    ind locID = cluster.localID();
    if (locID == Indices.size() - 1) {
        Indices.pop_back();
        Volumes.pop_back();
    } else {
        Holes.push_back(locID);
    }
}

inline void ClusterList::mergeClusters(ClusterID from, ClusterID onto) {
    ind locFrom = from.localID();
    ind locOnto = onto.localID();

    Volumes[locOnto] += Volumes[locFrom];
    removeCluster(locFrom);
}

}  // namespace perc