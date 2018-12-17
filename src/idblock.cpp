#include "idblock.h"
#include "clusterlist.h"

namespace perc {

IDBlock::IDBlock(const vec3i& size, const vec3i& offset, const vec3i& total, ID* memory)
    : BlockSize(size), BlockOffset(offset), TotalSize(total) {
    if (memory) {
        PointerBlock = memory;
        ownsMemory = false;
    } else {
        PointerBlock = new ID[BlockSize.prod()];
        ownsMemory = true;
    }
}

IDBlock::~IDBlock() {
    if (ownsMemory) delete[] PointerBlock;
}

}  // namespace perc