#pragma once
#include <vector>
#include <stack>
#include "vec.h"
#include "handles.h"
#include "clusterlist.h"
#include "unionfindsubblock.h"

namespace perc {

class GlobalProcessor;

struct GOG {
    VertexID ID;
    // TODO: Is a pointer to the parent better here?
    vec3i* ParentOffset;

    GOG(VertexID idx, vec3i* offset) : ID(idx), ParentOffset(offset){};
};

class ClusterListMultiple {
public:
    ClusterListMultiple(ind size = 100) : TotalVolume(0), MaxVolume(0) {
        IndicesPerCluster.reserve(100);
        Volumes.reserve(100);
    }

    double getClusterVolume(ClusterID cluster) {
        assert(std::find(Holes.begin(), Holes.end(), cluster) != Holes.end() &&
               "Trying to access non-existent cluster.");
        return Volumes[cluster.localID()];
    }

    void setRepresentative(ClusterID cluster, VertexID newID, bool replace = true,
                           vec3i* parentOffset = nullptr);

    ClusterID addCluster(VertexID id, double volume, vec3i* parentOffset = nullptr);
    void removeCluster(ClusterID cluster);
    void mergeClusters(ClusterID from, ClusterID onto);

    void extendCluster(ClusterID id, double volume, vec3i* parentOffset = nullptr);

    void clearVolumes();
    ind numClusters() { return IndicesPerCluster.size() - Holes.size(); }
    double totalVolume() { return TotalVolume; }
    double maxVolume() { return MaxVolume; }

private:
    std::vector<std::vector<GOG>> IndicesPerCluster;
    std::vector<double> Volumes;
    std::vector<size_t> Holes;
    double TotalVolume;
    double MaxVolume;
};

// ========= Inline Definitions ========= //

void ClusterListMultiple::setRepresentative(ClusterID cluster, VertexID newID, bool replace,
                                            vec3i* parentOffset) {
    assert(std::find(Holes.begin(), Holes.end(), cluster) != Holes.end() &&
           "Trying to access non-existent cluster.");
    assert(parentOffset && "No parent given.");
    // Check if cluster already has a representative for the given parent
    // TODO: Maybe this can be done outside?
    auto& GOGs = IndicesPerCluster[cluster.localID()];
    std::vector<GOG>::iterator it;
    for (it = GOGs.begin(); it != GOGs.end(); it++) {
        if (*it->ParentOffset == *parentOffset) {
            break;
        }
    }

    if (it != GOGs.end() && replace) {
        *it = {newID, parentOffset};
    } else {
        GOGs.emplace_back(newID, parentOffset);
    }
}

inline ClusterID ClusterListMultiple::addCluster(VertexID id, double volume, vec3i* parentOffset) {
    TotalVolume += volume;
    // The newly added cluster has a larger volume that other volumes so far.
    if (volume > MaxVolume) {
        MaxVolume = volume;
    }
    // Place the new cluster into a new hole or to the back when no holes exist.
    if (Holes.empty()) {
        IndicesPerCluster.emplace_back(id, parentOffset);
        Volumes.push_back(volume);
        return ClusterID(IndicesPerCluster.size() - 1);
    } else {
        size_t holeIdx = Holes.back();
        Holes.pop_back();
        std::vector<GOG> newGOG = {{id, parentOffset}};
        std::move(newGOG);
        IndicesPerCluster[holeIdx] = newGOG;
        Volumes[holeIdx] = volume;
        return ClusterID(holeIdx);
    }
}

inline void ClusterListMultiple::removeCluster(ClusterID cluster) {
    ind locID = cluster.localID();
    if (locID == IndicesPerCluster.size() - 1) {
        IndicesPerCluster.pop_back();
        Volumes.pop_back();
    } else {
        Holes.push_back(locID);
    }
}

inline void ClusterListMultiple::mergeClusters(ClusterID from, ClusterID onto) {
    ind locFrom = from.localID();
    ind locOnto = onto.localID();

    Volumes[locOnto] += Volumes[locFrom];
    if (Volumes[locOnto] > MaxVolume) {
        MaxVolume = Volumes[locOnto];
    }

    // Merge Representatives
    // TODO

    // All from clusters that are not part of the representatives anymore need to point to the new
    // cluster representative in their block-> Just make all of them point directly

    removeCluster(locFrom);
}

inline void ClusterListMultiple::extendCluster(ClusterID id, double volume, vec3i* parentOffset) {
    Volumes[id.localID()] += volume;
    TotalVolume += volume;
    if (Volumes[id.localID()] > MaxVolume) {
        MaxVolume = Volumes[id.localID()];
    }
    // Part of a block that does not have a GOG yet?
}

inline void ClusterListMultiple::clearVolumes() {
    std::fill(Volumes.data(), Volumes.data() + Volumes.size(), 0);
    TotalVolume = 0;
    MaxVolume = 0;
}

}  // namespace perc