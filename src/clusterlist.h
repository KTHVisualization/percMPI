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
    ClusterList(ind size = 100) : TotalVolume(0), MaxVolume(0) {
        Indices.reserve(100);
        Volumes.reserve(100);
    }

    Cluster getCluster(ClusterID cluster) {
        return Cluster(Indices[cluster.localID()], Volumes[cluster.localID()]);
    }
    ClusterID addCluster(VertexID id, double volume);
    void removeCluster(ClusterID cluster);
    void mergeClusters(ClusterID from, ClusterID onto);
    void extendCluster(ClusterID id, double volume);

    void clearVolumes();
    ind numClusters() { return Indices.size() - Holes.size(); }
    double totalVolume() { return TotalVolume; }
    double maxVolume() { return MaxVolume; }

private:
    std::vector<VertexID> Indices;
    std::vector<double> Volumes;
    std::vector<size_t> Holes;
    double TotalVolume;
    double MaxVolume;
};

// ========= Inline Definitions ========= //

inline ClusterID ClusterList::addCluster(VertexID id, double volume) {
    TotalVolume += volume;
    // The newly added cluster has a larger volume that other volumes so far.
    if (volume > MaxVolume) {
        MaxVolume = volume;
    }
    // Place the new cluster into a new hole or to the back when no holes exist.
    if (Holes.empty()) {
        Indices.push_back(id);
        Volumes.push_back(volume);
        return ClusterID(Indices.size() - 1);
    } else {
        size_t holeIdx = Holes.back();
        Holes.pop_back();
        Indices[holeIdx] = id;
        Volumes[holeIdx] = volume;
        return ClusterID(holeIdx);
    }
}
inline void ClusterList::removeCluster(ClusterID cluster) {
    ind locID = cluster.localID();
    TotalVolume -= Volumes[locID];
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
    if (Volumes[locOnto] > MaxVolume) {
        MaxVolume = Volumes[locOnto];
    }
    removeCluster(locFrom);
}

inline void ClusterList::extendCluster(ClusterID id, double volume) {
    Volumes[id.localID()] += volume;
    TotalVolume += volume;
    if (Volumes[id.localID()] > MaxVolume) {
        MaxVolume = Volumes[id.localID()];
    }
}

inline void ClusterList::clearVolumes() {
    std::fill(Volumes.data(), Volumes.data() + Volumes.size(), 0);
    TotalVolume = 0;
    MaxVolume = 0;
}

}  // namespace perc