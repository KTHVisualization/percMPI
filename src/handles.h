#pragma once
#include "vec.h"

namespace perc {

constexpr bool LOCAL_LIST = true;
constexpr bool GLOBAL_LIST = false;

struct ClusterID;
struct VertexID;

struct ID {
    static const ind CLUSTER_FLAG = ind(1) << (sizeof(ind) * 8 - 2);

    ID(ind id = -1) : RawID(id) {}

    bool isCluster() const { return RawID & CLUSTER_FLAG; }
    bool isVertex() const { return !isCluster(); }
    bool isValid() const { return RawID >= 0; }
    ClusterID* asCluster() { return isCluster() ? reinterpret_cast<ClusterID*>(this) : nullptr; }
    const ClusterID* asCluster() const {
        return isCluster() ? reinterpret_cast<const ClusterID*>(this) : nullptr;
    }

    ind baseID() const { return RawID & ~CLUSTER_FLAG; }
    ind RawID;

    static size_t hash(const ID& c) { return std::hash<ind>()(c.RawID); }
    typedef decltype(&hash) hash_type;
};

struct ClusterID : public ID {
    static const ind GLOBAL_FLAG = ind(1) << (sizeof(ind) * 8 - 3);

    ClusterID(ind id = -1, bool isLocal = true) : ID(id) {
        RawID = RawID | CLUSTER_FLAG;
        if (!isLocal) {
            RawID = RawID | GLOBAL_FLAG;
        }
    }

    bool isGlobal() const { return RawID & GLOBAL_FLAG; }
    ind localID() const { return baseID() & ~GLOBAL_FLAG; }
};

struct VertexID : public ID {
    VertexID(ind id = -1) : ID(id) {}
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