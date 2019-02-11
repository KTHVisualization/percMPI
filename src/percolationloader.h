#pragma once
#include "vec.h"
#include <string>
#include <unordered_map>

namespace perc {

class PercolationLoader {
public:
    PercolationLoader(const vec3i& blockSize, const vec3i& blockOffset)
        : BlockSize(blockSize), BlockOffset(blockOffset) {}

    double* loadScalarData();

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

    static void setSettings(vec3i totalSize, std::string directory, std::string rmsName,
                            ind timeStep) {
        TotalSizeFile = totalSize;
        RmsFilename = rmsName;
        Directory = directory;
        TimeStep = timeStep;
    }

    static void setTimeStep(ind timeStep) { TimeStep = timeStep; }

private:
    typedef double (*ScalarFunc)(const std::array<double, 3>&, const std::array<double, 3>&);
    static const std::unordered_map<std::string, ScalarFunc> ScalarVariants;

    vec3i BlockSize, BlockOffset;
    ScalarFunc RmsFunction;

    static vec3i TotalSizeFile;
    static std::string RmsFilename;
    static std::string Directory;
    static ind TimeStep;
};

}  // namespace perc