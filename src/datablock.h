#pragma once
#include <string>
#include "vec.h"

namespace perc {

class DataBlock {
public:
    DataBlock(const vec3i& size, const vec3i& offset, const vec3i& total)
        : BlockSize(size), BlockOffset(offset), TotalSize(total) {}
    ~DataBlock();

    bool loadData(ind timeSlice, const std::string& directory, const std::string& rmsFilename);

public:
    double *Scalar, *Volume;
    int* Indices;
    vec3i BlockSize, BlockOffset, TotalSize;
    // TODO: Neighboring node ids
};

}  // namespace perc