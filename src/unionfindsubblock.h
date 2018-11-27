#pragma once
#include "vec.h"
#include "datablock.h"
#include "unionfind.h"
#include "unionfindblock.h"

namespace perc {

class UnionFindSubBlock {
public:
    UnionFindSubBlock(const vec3i& size, const vec3i& offset, const vec3i& total,
                      UnionFindBlock& parent)
        : PointerBlock(size, offset, total), Parent(parent), CurrentWatershedIndex(-1) {}

    ~UnionFindSubBlock() { delete Data; }

    void loadData() {
        assert(!Data && "Data was already set.");
        Data =
            new DataBlock(PointerBlock.BlockSize, PointerBlock.BlockOffset, PointerBlock.TotalSize);

        CurrentWatershedIndex = 0;
        Data->sort();
    }

    void doWatershed(double maxVal) {
        if (!Data) return;
        // TODO
    }

    ClusterID* findClusterID(const vec3i& idx, vec3i& lastClusterID) {
        vec3i curIdx = idx;
        assert(contains(idx) && "First index is expected to be in here.");
        ID* firstPointer = PointerBlock.getPointer(idx);
        ID* curPointer = firstPointer;
        vec3i lastPointer = idx;

        ClusterID* finalID;
        while (curPointer && curPointer->isVertex()) {
            curIdx = vec3i::fromIndexOfTotal(curPointer->RawID, PointerBlock.TotalSize);

            if (!PointerBlock.contains(curIdx)) break;
            lastPointer = curIdx;

            ID* curPointer = PointerBlock.getPointer(curIdx);
            if (curPointer->isCluster()) {
                finalID = curPointer->asCluster();
                break;
            }
        }
        if (curPointer->isVertex()) finalID = Parent.findClusterID(idx, lastPointer);

        assert(finalID && finalID->baseID() >= 0 && "Invalid result of union find.");
        *firstPointer = lastPointer.toIndexOfTotal(PointerBlock.TotalSize);

        return finalID;
    }

    bool contains(const vec3i& idx) { return PointerBlock.contains(idx); }

public:
    DataBlock* Data;
    UnionFind PointerBlock;
    UnionFindBlock& Parent;
    ind CurrentWatershedIndex;
};

}  // namespace perc