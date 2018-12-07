#pragma once
#include "vec.h"

namespace perc {

struct ClusterID;
struct VertexID;

struct ID {
    static const ind CLUSTER_FLAG = 1 << (sizeof(ind) * 8 - 2);

    ID(ind id) : RawID(id) {}
    ID() : RawID(-1) {}

    bool isCluster() const { return RawID & CLUSTER_FLAG; }
    bool isVertex() const { return !isCluster(); }
    bool isValid() const { return RawID >= 0; }
    ClusterID* asCluster() { return isCluster() ? reinterpret_cast<ClusterID*>(this) : nullptr; }
    const ClusterID* asCluster() const {
        return isCluster() ? reinterpret_cast<const ClusterID*>(this) : nullptr;
    }

    ind baseID() const { return RawID & ~CLUSTER_FLAG; }
    ind RawID;
};

struct ClusterID : public ID {
    static const ind GLOBAL_FLAG = 1 << (sizeof(ind) * 8 - 3);

    ClusterID(ind id) : ID(id) { RawID = RawID | CLUSTER_FLAG; }
    ClusterID(ind id, bool isLocal) : ID(id) {
        RawID = RawID | CLUSTER_FLAG;
        if (!isLocal) {
            RawID = RawID | GLOBAL_FLAG;
        }
    }
    ClusterID() : ID() {}

    bool isGlobal() const { return RawID & GLOBAL_FLAG; }
    ind localID() const { return baseID() & ~GLOBAL_FLAG; }
};

struct VertexID : public ID {
    VertexID(ind id) : ID(id) {}
    VertexID() {}
};

struct Neighbor {
    Neighbor(ClusterID cluster, const vec3i& representative)
        : Cluster(cluster), Representative(representative) {}
    ClusterID Cluster;
    vec3i Representative;
};

inline bool operator==(const ID& a, const ID& b) { return a.RawID == b.RawID; }
inline bool operator==(const ClusterID& a, const ID& b) { return a.RawID == b.RawID; }
inline bool operator==(const ClusterID& a, const ClusterID& b) { return a.RawID == b.RawID; }
inline bool operator==(const VertexID& a, const ID& b) { return a.RawID == b.RawID; }
inline bool operator==(const VertexID& a, const VertexID& b) { return a.RawID == b.RawID; }

inline bool operator!=(const ID& a, const ID& b) { return a.RawID != b.RawID; }
inline bool operator!=(const ClusterID& a, const ID& b) { return a.RawID != b.RawID; }
inline bool operator!=(const ClusterID& a, const ClusterID& b) { return a.RawID != b.RawID; }
inline bool operator!=(const VertexID& a, const ID& b) { return a.RawID != b.RawID; }
inline bool operator!=(const VertexID& a, const VertexID& b) { return a.RawID != b.RawID; }

}  // namespace perc