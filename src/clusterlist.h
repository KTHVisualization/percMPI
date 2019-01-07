#pragma once
#include <vector>
#include <stack>
#include <cassert>

#include <iostream>
#include "vec.h"
#include "handles.h"

namespace perc {

struct Cluster {
    Cluster(VertexID& idx, double& vol) : Index(idx), Volume(vol) {}
    VertexID& Index;
    double& Volume;
};

struct ClusterData {
    ClusterData(VertexID& idx, double& vol) : Index(idx), Volume(vol) {}
    VertexID Index;
    double Volume;
};

class ClusterList {
public:
    ClusterList(bool isLocal = true, ind size = 100)
        : TotalVolume(0), MaxVolume(0), IsLocal(isLocal) {
        Indices.reserve(100);
        Volumes.reserve(100);
    }

    double getClusterVolume(ClusterID cluster) const;
    VertexID setRepresentative(ClusterID cluster, VertexID newID, bool replace = true,
                               void* parentBlock = nullptr);
    ClusterID addCluster(VertexID id, double volume, void* parentBlock = nullptr);
    void removeCluster(ClusterID cluster);
    void mergeClusters(ClusterID from, ClusterID onto);
    void extendCluster(ClusterID id, double volume, void* parentBlock = nullptr);
    Cluster getCluster(ClusterID id);

    void clearVolumes();
    ind numClusters() { return Indices.size() - Holes.size(); }
    double totalVolume() { return TotalVolume; }
    double maxVolume() { return MaxVolume; }

protected:
    void checkCluster(ClusterID cluster) const;

private:
    std::vector<VertexID> Indices;
    std::vector<double> Volumes;
    std::vector<size_t> Holes;
    double TotalVolume;
    double MaxVolume;
    const bool IsLocal;
};

// ========= Inline Definitions ========= //

inline double ClusterList::getClusterVolume(ClusterID cluster) const {
    checkCluster(cluster);

    return Volumes[cluster.localID()];
}

inline VertexID ClusterList::setRepresentative(ClusterID cluster, VertexID newID, bool replace,
                                               void*) {
    checkCluster(cluster);

    if (replace) {
        Indices[cluster.localID()] = newID;
        return newID;
    } else {
        return Indices[cluster.localID()];
    }
}

inline ClusterID ClusterList::addCluster(VertexID id, double volume, void*) {
    TotalVolume += volume;
    // The newly added cluster has a larger volume that other volumes so far.
    if (volume > MaxVolume) {
        MaxVolume = volume;
    }
    // Place the new cluster into a new hole or to the back when no holes exist.
    if (Holes.empty()) {
        Indices.push_back(id);
        Volumes.push_back(volume);
        return ClusterID(Indices.size() - 1, IsLocal);
    } else {
        size_t holeIdx = Holes.back();
        Holes.pop_back();
        Indices[holeIdx] = id;
        Volumes[holeIdx] = volume;
        return ClusterID(holeIdx, IsLocal);
    }
}
inline void ClusterList::removeCluster(ClusterID cluster) {
    checkCluster(cluster);

    ind locID = cluster.localID();
    if (locID == Indices.size() - 1) {
        Indices.pop_back();
        Volumes.pop_back();
    } else {
        Holes.push_back(locID);
    }
}

inline void ClusterList::mergeClusters(ClusterID from, ClusterID onto) {
    assert(from != onto && "Can't merge cluster onto itself.");
    checkCluster(from);
    checkCluster(onto);

    ind locFrom = from.localID();
    ind locOnto = onto.localID();

    Volumes[locOnto] += Volumes[locFrom];
    if (Volumes[locOnto] > MaxVolume) {
        MaxVolume = Volumes[locOnto];
    }
    removeCluster(from);
}

inline void ClusterList::extendCluster(ClusterID id, double volume, void*) {
    checkCluster(id);

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

inline Cluster ClusterList::getCluster(ClusterID id) {
    checkCluster(id);

    return Cluster(Indices[id.localID()], Volumes[id.localID()]);
}

inline void ClusterList::checkCluster(ClusterID cluster) const {
    assert(std::find(Holes.begin(), Holes.end(), cluster.localID()) == Holes.end() &&
           "Trying to access non-existent cluster.");
    assert(cluster.isGlobal() != IsLocal && "Local/global ID disagreement.");
}

}  // namespace perc