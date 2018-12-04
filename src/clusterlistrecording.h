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
    ClusterListRecording(ind size = 100) : Clusters(size), Merges() { Merges.reserve(20); }

    Cluster getCluster(ClusterID cluster) { return Clusters.getCluster(cluster); }
    void addClusters(ind numNewCLusters);
    void removeCluster(ClusterID cluster) { Clusters.removeCluster(cluster); }
    void mergeClusters(ClusterID from, ClusterID onto);
    void mergeClustersForReal(ClusterID from, ClusterID onto);
    void extendCluster(ClusterID id, double volume) { Clusters.extendCluster(id, volume); }

    void clearVolumesAndMerges();
    ind numClusters() { return Clusters.numClusters(); }

private:
    ClusterList Clusters;

    std::vector<ClusterMerge> Merges;
};

// ========= Inline Definitions ========= //
inline void ClusterListRecording::mergeClusters(ClusterID from, ClusterID onto) {
    Merges.push_back(ClusterMerge(from, onto));
}

inline void ClusterListRecording::mergeClustersForReal(ClusterID from, ClusterID onto) {
    Clusters.mergeClusters(from, onto);
}

inline void ClusterListRecording::clearVolumesAndMerges() {
    Clusters.clearVolumes();
    Merges.clear();
}

inline void ClusterListRecording::addClusters(ind numNewClusters) {
    for (ind n = 0; n < numNewClusters; ++n) Clusters.addCluster(VertexID(-1), 0.0);
}

}  // namespace perc