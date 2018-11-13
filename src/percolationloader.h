#pragma once
#include "vec.h"
#include <string>
#include <unordered_map>

namespace perc {

class PercolationLoader {
public:
    PercolationLoader(const vec3i& blockSize, const vec3i& blockOffset, const vec3i& totalSize)
        : BlockSize(blockSize), BlockOffset(blockOffset), TotalSize(totalSize) {}

    double* loadScalarData(ind timeStep, const std::string& pathVelocity,
                           const std::string& pathAverage, const std::string& pathRms,
                           const std::string& pathVertex);

    double* loadBlock(const std::string& path, bool is2D = false) const;

    double* normalizedFromComponents(const std::array<double*, 3> velocity,
                                     const std::array<double*, 3> average, const double* rms);

    void getRmsTypeFromFilename(const std::string& rmsName);

    inline static bool fileExists(const std::string& name) {
        if (FILE* file = fopen(name.c_str(), "r")) {
            fclose(file);
            return true;
        } else {
            return false;
        }
    }

private:
    typedef double (*ScalarFunc)(const std::array<double, 3>&, const std::array<double, 3>&);
    static const std::unordered_map<std::string, ScalarFunc> ScalarVariants;

    vec3i BlockSize, BlockOffset, TotalSize;
    ScalarFunc RmsFunction;
};

}  // namespace perc