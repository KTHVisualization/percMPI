#include "datablock.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include "percolationloader.h"

namespace perc {

bool DataBlock::loadData(ind timeSlice, const std::string& directory,
                         const std::string& rmsFilename) {

    // Init loader with size.
    PercolationLoader loader(BlockSize, BlockOffset, TotalSize);
    Scalar = loader.loadScalarData(timeSlice, directory + "/VELOCITY/", directory + "/STAT/",
                                   directory + "/ZEXPORT_STAT_wall_correction/" + rmsFilename,
                                   directory + "/STAT/");
    if (!Scalar) return false;

    // Volume just set to uniform for now.
    Volume = new double[BlockSize.prod()];
    std::fill_n(Volume, BlockSize.prod(), 1.0);

    return true;
}

}  // namespace perc