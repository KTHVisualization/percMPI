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
    std::fill_n(PointerBlock, BlockSize.prod(), ind(-1));
}

IDBlock::IDBlock(IDBlock&& other) noexcept
    : BlockSize(other.BlockSize)
    , BlockOffset(other.BlockOffset)
    , TotalSize(other.TotalSize)
    , ownsMemory(other.ownsMemory) {
    PointerBlock = other.PointerBlock;
    other.PointerBlock = nullptr;
    other.ownsMemory = false;
}

IDBlock::~IDBlock() {
    if (ownsMemory) delete[] PointerBlock;
}

void IDBlock::reset() { std::fill_n(PointerBlock, BlockSize.prod(), ind(-1)); }

}  // namespace perc