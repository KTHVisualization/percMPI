#pragma once
#include <vector>
#include <stack>
#include "vec.h"
#include "handles.h"
#include "clusterlist.h"

namespace perc {

struct ClusterMerge {
    ClusterMerge(ClusterID from, ClusterID onto) : From(from), Onto(onto) {}
    ClusterID From, Onto;
};

class ClusterListRecording {
public:
    ClusterListRecording(ind size = 100) : Clusters(size), Merges(20) {}

    Cluster getCluster(ClusterID cluster) { return Clusters.getCluster(cluster); }
    inline ClusterID addCluster(VertexID id) { return Clusters.addCluster(id); }
    inline void removeCluster(ClusterID cluster) { Clusters.removeCluster(cluster); }
    inline void mergeClusters(ClusterID from, ClusterID onto);

    inline void clearVolumesAndMerges();

private:
    ClusterList Clusters;

    std::vector<ClusterMerge> Merges;
};

// ========= Inline Definitions ========= //
inline void ClusterListRecording::mergeClusters(ClusterID from, ClusterID onto) {
    Merges.push_back(ClusterMerge(from, onto));
}

inline void ClusterListRecording::clearVolumesAndMerges() {
    Clusters.clearVolumes();
    Merges.clear();
}

}  // namespace perc