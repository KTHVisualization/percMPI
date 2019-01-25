#pragma once
#include "vec.h"
#include "datablock.h"
#include "handles.h"
#include <cassert>

namespace perc {

class IDBlock {
public:
    IDBlock(const vec3i& size, const vec3i& offset, const vec3i& total, ID* memory = nullptr);
    ~IDBlock();

    bool contains(const vec3i& idx) {
        vec3i max = BlockOffset + BlockSize;
        return idx.liesWithin(BlockOffset, max);
    }

    ID* getPointer(const vec3i& idx) {
        assert(contains(idx) && "Index not within this block. Please check before!\n");
        return getPointerLocal(idx - BlockOffset);
    }

    ID* getPointerLocal(const vec3i& locIdx) {
        ID* ptr = PointerBlock + locIdx.toIndexOfTotal(BlockSize);
        return ptr->isValid() ? ptr : nullptr;
    }

    void setPointer(const vec3i& idx, const ID& id) {
        assert(id.isCluster() || VertexID(idx.toIndexOfTotal(TotalSize)) != id &&
                                     "Vertex attempting to point to itself.\n");
        assert(contains(idx) && "Index not within this block. Please check before!\n");
        return setPointerLocal(idx - BlockOffset, id);
    }

    void setPointerLocal(const vec3i& locIdx, const ID& id) {
        ID* ptr = PointerBlock + locIdx.toIndexOfTotal(BlockSize);
        *ptr = id;
    }

public:
    const vec3i BlockSize, BlockOffset, TotalSize;
    ID* PointerBlock;
    bool ownsMemory;
};

}  // namespace perc