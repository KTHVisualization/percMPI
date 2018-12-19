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
    void* ParentBlock;

    GOG(VertexID idx, void* parent) : ID(idx), ParentBlock(parent){};

    bool operator<(const GOG& rhs) const { return ParentBlock < rhs.ParentBlock; }
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

    const std::vector<GOG> getRepresentatives(ClusterID cluster);
    VertexID setRepresentative(ClusterID cluster, VertexID newID, bool replace = true,
                               void* parentBlock = nullptr);

    const std::vector<VertexID> mergeRepresentatives(ClusterID from, ClusterID onto);

    ClusterID addCluster(VertexID id, double volume, void* parentBlock = nullptr);
    void removeCluster(ClusterID cluster);
    void mergeClusters(ClusterID from, ClusterID onto);

    void extendCluster(ClusterID id, double volume, void* parentBlock = nullptr);

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

inline const std::vector<GOG> ClusterListMultiple::getRepresentatives(ClusterID cluster) {
    return IndicesPerCluster[cluster.localID()];
}

inline VertexID ClusterListMultiple::setRepresentative(ClusterID cluster, VertexID newID,
                                                       bool replace, void* parentBlock) {
    assert(std::find(Holes.begin(), Holes.end(), cluster) != Holes.end() &&
           "Trying to access non-existent cluster.");
    assert(parentBlock && "No parent given.");
    GOG newRep = {newID, parentBlock};
    auto& GOGs = IndicesPerCluster[cluster.localID()];
    auto it = std::lower_bound(GOGs.begin(), GOGs.end(), newRep);

    if (it != GOGs.end() && it->ParentBlock == newRep.ParentBlock) {
        if (replace) {
            *it = newRep;
        } else {
            return it->ID;
        }
    } else {
        GOGs.insert(it, newRep);
    }
    return newID;
}

inline const std::vector<VertexID> ClusterListMultiple::mergeRepresentatives(ClusterID from,
                                                                             ClusterID onto) {
    ind locFrom = from.localID();
    ind locOnto = onto.localID();
    auto& ontoIDs = IndicesPerCluster[locOnto];
    auto& fromIDs = IndicesPerCluster[locFrom];
    std::vector<GOG> newRepsOnto;
    std::vector<VertexID> newRepsFrom;
    newRepsFrom.reserve(ontoIDs.size());

    auto itOnto = ontoIDs.begin();
    auto itFrom = fromIDs.begin();

    // Both lists of representatives are sorted, we merge the lists
    // If both lists have a representative for a block, we keep the rep in ontoIDs
    // and mark it as what the id in fromIDs should point to
    while (itOnto != ontoIDs.end() && itFrom != fromIDs.end()) {
        if (*itOnto < *itFrom) {
            newRepsOnto.push_back(*itOnto);
            itOnto++;
        } else if (*itFrom < *itOnto) {
            newRepsFrom.push_back(itFrom->ID);
            newRepsOnto.push_back(*itFrom);
            itFrom++;
        } else {
            newRepsFrom.push_back(itOnto->ID);
            newRepsOnto.push_back(*itOnto);
            itOnto++;
            itFrom++;
        }
    }

    // Copy remaining ones over
    while (itOnto != ontoIDs.end()) {
        newRepsOnto.push_back(*itOnto);
        itOnto++;
    }

    while (itFrom != fromIDs.end()) {
        newRepsFrom.push_back(itFrom->ID);
        newRepsOnto.push_back(*itFrom);
        itFrom++;
    }

    IndicesPerCluster[locOnto] = std::move(newRepsOnto);

    return newRepsFrom;
}

inline ClusterID ClusterListMultiple::addCluster(VertexID id, double volume, void* parentBlock) {
    TotalVolume += volume;
    // The newly added cluster has a larger volume that other volumes so far.
    if (volume > MaxVolume) {
        MaxVolume = volume;
    }
    std::vector<GOG> newGOG = {GOG(id, parentBlock)};
    // Place the new cluster into a new hole or to the back when no holes exist.
    if (Holes.empty()) {
        IndicesPerCluster.push_back(std::move(newGOG));
        Volumes.push_back(volume);
        return ClusterID(IndicesPerCluster.size() - 1);
    } else {
        size_t holeIdx = Holes.back();
        Holes.pop_back();
        IndicesPerCluster[holeIdx] = std::move(newGOG);
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

    removeCluster(locFrom);
}

inline void ClusterListMultiple::extendCluster(ClusterID id, double volume, void* parentBlock) {
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