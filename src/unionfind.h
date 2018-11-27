#pragma once
#include "vec.h"
#include "datablock.h"
#include "handles.h"
#include <cassert>

namespace perc {

class UnionFind {
public:
    UnionFind(const vec3i& size, const vec3i& offset, const vec3i& total);
    ~UnionFind() { delete[] PointerBlock; }

    bool contains(const vec3i& idx) { return idx.liesWithin(BlockOffset, BlockOffset + BlockSize); }
    ID* getPointer(const vec3i& idx) {
        assert(contains(idx) && "Index not within this block. Please check before!\n");
        ID* ptr = PointerBlock + (idx - BlockOffset).toIndexOfTotal(TotalSize);
        return ptr->isValid() ? ptr : nullptr;
    }

public:
    vec3i BlockSize, BlockOffset, TotalSize;
    ID* PointerBlock;
};

}  // namespace perc