#include "unionfindblock.h"
#include "unionfindsubblock.h"

namespace perc {

ClusterID* UnionFindBlock::findClusterID(const vec3i& idx, vec3i& lastClusterID) {
    // One might want to do this more cleverly, especialy in the sheet tree.
    for (UnionFindSubBlock* sub : SubBlocks)
        if (sub->contains(idx)) return sub->findClusterID(idx, lastClusterID);

    assert(false && "Can not find block containing this idx.");
    return nullptr;
}

}  // namespace perc