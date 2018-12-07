#include "idblock.h"
#include "clusterlist.h"

namespace perc {

IDBlock::IDBlock(const vec3i& size, const vec3i& offset, const vec3i& total)
    : BlockSize(size), BlockOffset(offset), TotalSize(total) {
    PointerBlock = new ID[BlockSize.prod()];
}

}  // namespace perc