#pragma once
#include <vector>
#include <unordered_set>
#include "vec.h"
#include "handles.h"
#include "clusterlist.h"

namespace perc {

struct ClusterMerge {
    ClusterMerge(ClusterID from, ClusterID onto) : From(from), Onto(onto) {}

    static inline std::vector<std::vector<ind>> mergeClustersFromLists(
        const std::vector<std::vector<ClusterMerge>>& merges);
    static inline std::vector<ind> mergeClusterAsList(
        const std::vector<std::vector<ind>>& mergeClusters);

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

std::vector<std::vector<ind>> ClusterMerge::mergeClustersFromLists(
    const std::vector<std::vector<ClusterMerge>>& merges) {

    // Resulting list of components.
    std::vector<std::vector<ind>> mergeClusters;
    // Clusters that are not yet merged.
    std::unordered_set<ind> todoClusters;

    for (const std::vector<ClusterMerge>& mergeList : merges)
        for (const ClusterMerge& merge : mergeList) {
            todoClusters.insert(merge.From.RawID);
            todoClusters.insert(merge.Onto.RawID);
        }

    while (!todoClusters.empty()) {
        ind cluster = *todoClusters.begin();
        todoClusters.erase(todoClusters.begin());

        // Add new list of clusters to be merged.
        mergeClusters.emplace_back();
        std::vector<ind>& graph = mergeClusters.back();

        graph.push_back(cluster);

        bool addedSomething = false;
        do {
            addedSomething = false;
            for (const std::vector<ClusterMerge>& mergeList : merges)
                for (const ClusterMerge& merge : mergeList) {
                    bool fromInGraph =
                        std::find(graph.begin(), graph.end(), merge.From.RawID) != graph.end();
                    bool ontoInGraph =
                        std::find(graph.begin(), graph.end(), merge.Onto.RawID) != graph.end();

                    // Either both in graph already, or neither.
                    if (fromInGraph == ontoInGraph) continue;

                    if (fromInGraph) {
                        graph.push_back(merge.Onto.RawID);

                        assert(std::find(todoClusters.begin(), todoClusters.end(),
                                         merge.Onto.RawID) != todoClusters.end() &&
                               "Invalid constellation in cluster merge.");
                        todoClusters.erase(
                            std::find(todoClusters.begin(), todoClusters.end(), merge.Onto.RawID));
                    } else {
                        graph.push_back(merge.From.RawID);

                        assert(std::find(todoClusters.begin(), todoClusters.end(),
                                         merge.From.RawID) != todoClusters.end() &&
                               "Invalid constellation in cluster merge.");
                        todoClusters.erase(
                            std::find(todoClusters.begin(), todoClusters.end(), merge.From.RawID));
                    }
                    addedSomething = true;
                }
        } while (addedSomething);
    }
    return mergeClusters;
}

std::vector<ind> ClusterMerge::mergeClusterAsList(
    const std::vector<std::vector<ind>>& mergeClusters) {
    std::vector<ind> clusterList;
    for (auto& graph : mergeClusters) {
        clusterList.push_back(graph.size());
        for (ind val : graph) clusterList.push_back(val);
    }
    return clusterList;
}

}  // namespace perc