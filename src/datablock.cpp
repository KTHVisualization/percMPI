#include "datablock.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include "percolationloader.h"

namespace perc {

DataBlock::~DataBlock() {
    delete[] Scalars;
    delete[] Volumes;
}

bool DataBlock::loadData(ind timeSlice, const std::string& directory,
                         const std::string& rmsFilename) {

    // Init loader with size.
    PercolationLoader loader(BlockSize, BlockOffset, TotalSize);
    Scalars = loader.loadScalarData(timeSlice, directory + "/VELOCITY/", directory + "/STAT/",
                                    directory + "/ZEXPORT_STAT_wall_correction/" + rmsFilename,
                                    directory + "/STAT/");
    if (!Scalars) return false;

    // Volume just set to uniform for now.
    Volumes = new double[BlockSize.prod()];
    std::fill_n(Volumes, BlockSize.prod(), 1.0);

    return true;
}

void DataBlock::sort() {
    delete[] Indices;

    // Create Indices, fill with [0, numElements) and sort by Scalar value.
    ind numElements = BlockSize.prod();
    Indices = new ind[numElements];
    std::iota(Indices, Indices + numElements, 0);
    std::sort(Indices, Indices + numElements,
              [this](ind a, ind b) { return Scalars[a] > Scalars[b]; });
}

}  // namespace perc