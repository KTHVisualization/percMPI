#pragma once
#include <vector>
#include <unordered_set>
#include <map>
#include <list>
#include "vec.h"
#include "handles.h"
#include "clusterlist.h"
#include "clusterlistmultiple.h"

namespace perc {

struct ClusterMerge {
    ClusterMerge(ClusterID from, ClusterID onto) : From(from), Onto(onto) {}
    ClusterMerge() : From(-1), Onto(-1) {}

    static inline std::vector<std::vector<ind>> mergeClustersFromLists(
        const std::vector<std::vector<ClusterMerge>*>& merges);
    static inline std::vector<std::vector<ind>> mergeClustersFromLists(
        std::vector<std::vector<ClusterMerge>>& merges);
    static inline std::vector<ind> mergeClusterAsList(
        const std::vector<std::vector<ind>>& mergeClusters);

    ClusterID From, Onto;
};

template <typename CL>
class ClusterListRecording {
public:
    ClusterListRecording(bool isLocal = false, ind size = 100) : Clusters(isLocal, size), Merges() {
        Merges.reserve(20);
    }

    double getClusterVolume(ClusterID cluster) const { return Clusters.getClusterVolume(cluster); }
    void addClusters(ind numNewClusters);
    ClusterID addCluster(double volume);
    ClusterID addCluster(VertexID id, double volume, void* parentBlock = nullptr);
    VertexID getRepresentative(ClusterID cluster) { return Clusters.getCluster(cluster).Index; }
    VertexID setRepresentative(ClusterID cluster, VertexID newID, bool replace = true,
                               void* parentBlock = nullptr);
    void removeCluster(ClusterID cluster) { Clusters.removeCluster(cluster); }
    void mergeClusters(ClusterID from, ClusterID onto);
    void mergeClustersForReal(ClusterID from, ClusterID onto);
    void mergeClusterFromList(const std::vector<ind>& connectedComps);
    void extendCluster(ClusterID id, double volume) { Clusters.extendCluster(id, volume); }

    void clearVolumesAndMerges();
    const std::vector<double>& volumes();
    ind numClusters() { return Clusters.numClusters(); }
    double totalVolume() { return Clusters.totalVolume(); };
    double maxVolume() { return Clusters.maxVolume(); };

protected:
    CL Clusters;

    std::vector<ClusterMerge> Merges;
};

class ClusterListRecordingSingle : public ClusterListRecording<ClusterList> {
    friend class LocalBlock;

public:
    ClusterListRecordingSingle(bool isLocal = true, ind size = 100)
        : ClusterListRecording<ClusterList>(isLocal, size) {}
    Cluster getCluster(ClusterID id) { return Clusters.getCluster(id); }
};

class ClusterListRecordingMultiple : public ClusterListRecording<ClusterListMultiple> {
    friend class GlobalBlock;

public:
    ClusterListRecordingMultiple(bool isLocal = true, ind size = 100)
        : ClusterListRecording<ClusterListMultiple>(isLocal, size) {}

    const std::vector<GOG> getRepresentatives(ClusterID cluster) {
        return Clusters.getRepresentatives(cluster);
    }

    const std::vector<VertexID> mergeRepresentatives(ClusterID from, ClusterID onto) {
        return Clusters.mergeRepresentatives(from, onto);
    }
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
inline void ClusterListRecording<CL>::mergeClusterFromList(const std::vector<ind>& connComps) {

    for (auto it = connComps.begin(); it != connComps.end(); ++it) {
        ind compSize = *it;
        ind onto = *(++it);
        for (ind c = 0; c < compSize - 1; ++c) mergeClustersForReal(*(++it), onto);
    }
}

template <typename CL>
inline void ClusterListRecording<CL>::clearVolumesAndMerges() {
#ifdef COMMUNICATION
    Clusters.clearVolumes();
#endif
    Merges.clear();
}

template <typename CL>
const std::vector<double>& ClusterListRecording<CL>::volumes() {
    return Clusters.volumes();
}

template <typename CL>
inline void ClusterListRecording<CL>::addClusters(ind numNewClusters) {
    for (ind n = 0; n < numNewClusters; ++n) Clusters.addCluster(VertexID(-1), 0.0);
}

template <typename CL>
inline ClusterID ClusterListRecording<CL>::addCluster(double volume) {
    return Clusters.addCluster(VertexID(-1), volume);
}

template <typename CL>
inline ClusterID ClusterListRecording<CL>::addCluster(VertexID id, double volume,
                                                      void* parentBlock) {
    return Clusters.addCluster(id, volume, parentBlock);
}

template <typename CL>
inline VertexID ClusterListRecording<CL>::setRepresentative(ClusterID cluster, VertexID newID,
                                                            bool replace, void* parentBlock) {
    return Clusters.setRepresentative(cluster, newID, replace, parentBlock);
}

std::vector<std::vector<ind>> ClusterMerge::mergeClustersFromLists(
    std::vector<std::vector<ClusterMerge>>& merges) {
    std::vector<std::vector<ClusterMerge>*> mergesPointers;
    mergesPointers.reserve(merges.size());
    for (std::vector<ClusterMerge>& mergeList : merges) mergesPointers.push_back(&mergeList);
    return mergeClustersFromLists(mergesPointers);
}

std::vector<std::vector<ind>> ClusterMerge::mergeClustersFromLists(
    const std::vector<std::vector<ClusterMerge>*>& merges) {

    // Resulting list of components.
    std::vector<std::vector<ind>> mergeClusters;
    // Clusters that are not yet merged.
    std::unordered_set<ind> todoClusters;

    // Convert edge list (potentially with duplicates) into a graph
    std::map<ind, std::unordered_set<ind>> mergeGraph;
    for (const std::vector<ClusterMerge>* pMergeList : merges) {
        const std::vector<ClusterMerge>& mergeList = *pMergeList;
        for (const ClusterMerge& merge : mergeList) {
            ind from = merge.From.RawID;
            ind onto = merge.Onto.RawID;
            auto itVertex = mergeGraph.lower_bound(from);
            if (itVertex != mergeGraph.end() && itVertex->first == from) {
                itVertex->second.insert(onto);
            } else {
                mergeGraph.emplace_hint(itVertex, from, std::initializer_list<ind>({onto}));
            }
            todoClusters.insert(from);
            itVertex = mergeGraph.lower_bound(onto);
            if (itVertex != mergeGraph.end() && itVertex->first == onto) {
                itVertex->second.insert(from);
            } else {
                mergeGraph.emplace_hint(itVertex, onto, std::initializer_list<ind>({from}));
            }
            todoClusters.insert(onto);
        }
    }

    while (!todoClusters.empty()) {
        ind cluster = *todoClusters.begin();
        todoClusters.erase(todoClusters.begin());

        // Add new list of clusters to be merged.
        mergeClusters.emplace_back();
        std::vector<ind>& component = mergeClusters.back();

        component.push_back(cluster);
        std::list<ind> queue;
        ind current;
        queue.push_back(cluster);

        while (!queue.empty()) {
            current = queue.front();
            queue.pop_front();

            auto neighbors = mergeGraph.at(current);
            for (auto& neigh : neighbors) {
                auto neighIt = std::find(todoClusters.begin(), todoClusters.end(), neigh);
                if (neighIt != todoClusters.end()) {
                    todoClusters.erase(neighIt);
                    queue.push_back(neigh);
                    component.push_back(neigh);
                }
            }
        }
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