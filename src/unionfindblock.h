#pragma once
#include <vector>
#include "vec.h"
#include "handles.h"

namespace perc {

class UnionFindSubBlock;

class UnionFindBlock {
public:
    ClusterID* findClusterID(const vec3i& idx, vec3i& lastClusterID);

protected:
    std::vector<UnionFindSubBlock*> SubBlocks;
};

}  // namespace perc