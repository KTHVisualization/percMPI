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

template <typename CL>
class ClusterListRecording {
public:
    ClusterListRecording(ind size = 100) : Clusters(size), Merges() { Merges.reserve(20); }

    double getClusterVolume(ClusterID cluster) { return Clusters.getClusterVolume(cluster); }
    void addClusters(ind numNewCLusters);
    void removeCluster(ClusterID cluster) { Clusters.removeCluster(cluster); }
    void mergeClusters(ClusterID from, ClusterID onto);
    void mergeClustersForReal(ClusterID from, ClusterID onto);
    void extendCluster(ClusterID id, double volume) { Clusters.extendCluster(id, volume); }

    void clearVolumesAndMerges();
    ind numClusters() { return Clusters.numClusters(); }

private:
    CL Clusters;

    std::vector<ClusterMerge> Merges;
};

// ========= Inline Definitions ========= //
template <typename CL>
inline void ClusterListRecording<CL>::mergeClusters(ClusterID from, ClusterID onto) {
    Merges.push_back(ClusterMerge(from, onto));
}

template <typename CL>
inline void ClusterListRecording<CL>::mergeClustersForReal(ClusterID from, ClusterID onto) {
    Clusters.mergeClusters(from, onto);
}

template <typename CL>
inline void ClusterListRecording<CL>::clearVolumesAndMerges() {
    Clusters.clearVolumes();
    Merges.clear();
}

template <typename CL>
inline void ClusterListRecording<CL>::addClusters(ind numNewClusters) {
    for (ind n = 0; n < numNewClusters; ++n) Clusters.addCluster(VertexID(-1), 0.0);
}

}  // namespace perc