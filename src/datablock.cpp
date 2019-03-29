#include "datablock.h"
#include <algorithm>
#include <cmath>
#include "percolationloader.h"
#include <iostream>

namespace perc {

ind DataBlock::NumThresholds = 10;
double DataBlock::ThresholdMin = 0.0;
double DataBlock::ThresholdMax = 1.0;

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

void DataBlock::sortData(bool useBuckets) {
    assert(!Indices && "Was already sorted?");

    // Create Indices, fill with [0, numElements) and sort by Scalar value.
    ind numElements = BlockSize.prod();
    Indices = new ind[numElements];

    if (!useBuckets || numElements < NumThresholds * 10) {
        std::iota(Indices, Indices + numElements, 0);
        // Sorts from largest to smallest value.
        std::sort(Indices, Indices + numElements,
                  [this](ind a, ind b) { return Scalars[a] > Scalars[b]; });
    } else {
        std::vector<std::pair<double, std::vector<ind>>> buckets(NumThresholds);

        // Setup vector of thresholds and buckets.
        // Yes, the original loo uses a mix of float and double.
        float hStep = (ThresholdMax - ThresholdMin) / (NumThresholds - 1);
        ind index = 0;
        for (float currentH = ThresholdMax; currentH >= ThresholdMin - 1e-5; currentH -= hStep) {
            buckets[index].first = currentH;
            ++index;
        }

        // Bucket all values.
        double hMin = buckets[NumThresholds - 1].first;
        ind indexOfAnyValueBelow = -1;
        for (ind e = 0; e < numElements; ++e) {
            double val = Scalars[e];
            if (val < hMin) {
                indexOfAnyValueBelow = e;
                continue;
            }
            if (val >= ThresholdMax) {
                buckets[0].second.push_back(e);
                continue;
            }
            ind probableIdx = std::ceil((ThresholdMax - val) / hStep);
            probableIdx = std::min(probableIdx, NumThresholds - 1);

            // Value too low.
            if (val < buckets[probableIdx].first) {
                probableIdx++;
            }
            // Value too large.
            else if (probableIdx != 0 && val >= buckets[probableIdx - 1].first) {
                probableIdx--;
            }

            assert(val >= buckets[probableIdx].first &&
                   (probableIdx == 0 || val < buckets[probableIdx - 1].first) &&
                   "Did not find correct bucket.");

            buckets[probableIdx].second.push_back(e);
        }

        ind indexIdx = 0;
        for (auto& bucket : buckets)
            for (ind idx : bucket.second) Indices[indexIdx++] = idx;

        if (indexIdx < numElements) Indices[indexIdx] = indexOfAnyValueBelow;
    }
}

vec3i DataBlock::toGlobalIndex(ind locIdx) {
    assert(locIdx >= 0 && locIdx < BlockSize.prod() && "Index outside the block size.");
    vec3i relIdx = vec3i::fromIndexOfTotal(locIdx, BlockSize);
    return BlockOffset + relIdx;
}

ind DataBlock::memEstimate() const {
    // Indices for Sorting
    ind memSize = BlockSize.prod() * sizeof(ind);
    // Volumes and Scalars
    memSize += BlockSize.prod() * sizeof(double) * 2;
    return memSize;
}

}  // namespace perc