#pragma once
#include "vec.h"

namespace perc {

struct ClusterID;

struct ID {
    static const ind CLUSTER_FLAG = 1 << (sizeof(ind) * 8 - 2);

    ID(ind id) : RawID(id) {}
    ID() : RawID(-1) {}

    bool isCluster() { return RawID & CLUSTER_FLAG; }
    bool isVertex() { return !isCluster(); }
    bool isValid() { return RawID >= 0; }
    ClusterID* asCluster() { return isCluster() ? reinterpret_cast<ClusterID*>(this) : nullptr; }

    ind baseID() { return RawID & ~CLUSTER_FLAG; }
    ind RawID;
};

struct ClusterID : public ID {
    static const ind GLOBAL_FLAG = 1 << (sizeof(ind) * 8 - 3);
    ClusterID(ind id) : ID(id) {}

    bool isGlobal() { return RawID & GLOBAL_FLAG; }
    ind localID() { return baseID() & ~GLOBAL_FLAG; }
};

struct VertexID : public ID {
    VertexID(ind id) : ID(id) {}
};

}  // namespace perc