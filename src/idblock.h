#pragma once
#include "vec.h"
#include "datablock.h"
#include "handles.h"
#include <cassert>

namespace perc {

class IDBlock {
public:
    IDBlock(const vec3i& size, const vec3i& offset, const vec3i& total);
    ~IDBlock() { delete[] PointerBlock; }

    bool contains(const vec3i& idx) { return idx.liesWithin(BlockOffset, BlockOffset + BlockSize); }

    ID* getPointer(const vec3i& idx) {
        assert(contains(idx) && "Index not within this block. Please check before!\n");
        return getPointerLocal(idx - BlockOffset);
    }

    ID* getPointerLocal(const vec3i& locIdx) {
        ID* ptr = PointerBlock + locIdx.toIndexOfTotal(TotalSize);
        return ptr->isValid() ? ptr : nullptr;
    }

    void setPointer(const vec3i& idx, const ID& id) {
        assert(contains(idx) && "Index not within this block. Please check before!\n");
        return setPointerLocal(idx - BlockOffset, id);
    }

    void setPointerLocal(const vec3i& locIdx, const ID& id) {
        ID* ptr = PointerBlock + locIdx.toIndexOfTotal(TotalSize);
        *ptr = id;
    }

public:
    const vec3i BlockSize, BlockOffset, TotalSize;
    ID* PointerBlock;
};

}  // namespace perc