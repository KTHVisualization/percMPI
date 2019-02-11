#pragma once
#include <string>
#include <numeric>
#include "vec.h"

namespace perc {

class DataBlock {
public:
    DataBlock(const vec3i& size, const vec3i& offset, const vec3i& total)
        : BlockSize(size), BlockOffset(offset), TotalSize(total), Indices(nullptr) {}
    ~DataBlock();

    bool loadData();
    void sort();
    vec3i toGlobalIndex(ind locIdx);

public:
    double *Scalars, *Volumes;
    ind* Indices;
    const vec3i BlockSize, BlockOffset, TotalSize;
};

}  // namespace perc