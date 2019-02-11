#include "datablock.h"
#include <algorithm>
#include <cmath>
#include "percolationloader.h"

namespace perc {

DataBlock::~DataBlock() {
    delete[] Scalars;
    delete[] Volumes;
    delete[] Indices;
}

bool DataBlock::loadData() {

    // Init loader with size.
    PercolationLoader loader(BlockSize, BlockOffset);
    Scalars = loader.loadScalarData();
    if (!Scalars) return false;

    // Volume just set to uniform for now.
    Volumes = new double[BlockSize.prod()];
    std::fill_n(Volumes, BlockSize.prod(), 1.0);

    return true;
}

void DataBlock::sort() {
    assert(!Indices && "Was already sorted?");

    // Create Indices, fill with [0, numElements) and sort by Scalar value.
    ind numElements = BlockSize.prod();
    Indices = new ind[numElements];
    std::iota(Indices, Indices + numElements, 0);
    // Sorts from largest to smallest value.
    std::sort(Indices, Indices + numElements,
              [this](ind a, ind b) { return Scalars[a] > Scalars[b]; });
}

vec3i DataBlock::toGlobalIndex(ind locIdx) {
    assert(locIdx >= 0 && locIdx < BlockSize.prod() && "Index outside the block size.");
    vec3i relIdx = vec3i::fromIndexOfTotal(locIdx, BlockSize);
    return BlockOffset + relIdx;
}

}  // namespace perc