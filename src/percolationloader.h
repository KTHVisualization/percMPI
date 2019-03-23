#pragma once
#include "vec.h"
#include <string>
#include <unordered_map>

namespace perc {

enum class InputMode {
    INVALID = -1,
    // Velocity as  timestep.x, timestep.y, timestep.z,
    // Wall correction filename determines scalar
    // Averages for u/w/v
    // Rms file for product of two components uv/uw/..
    COMBINED_VELOCITY_AVG_RMS_FILE = 0,
    // Same as before except that we have two Rms files (e.g. components separately)
    COMBINED_VELOCITY_AVG_2RMS_FILE = 1,
    // Run on velocity product directly
    VELOCITY_FILE = 2,
    // Rms in each component and avgs as single value (Isotropic dataset)
    COMBINED_VELOCITY_AVG_RMS_VALUE = 3,

    // Directly load a scalar file (.vti)
    SCALAR = 10,
    // Uniformly random data between 0 and 1
    RANDOM_UNIFORM = 20,
    SHUFFLED_COMBINED_VELOCITY_AVG_RMS_FILE = 100,
    SHUFFLED_VELOCITY_FILE = 112,
    SHUFFLED_SCALAR = 110,
};

class PercolationLoader {
public:
    PercolationLoader(const vec3i& blockSize, const vec3i& blockOffset)
        : BlockSize(blockSize), BlockOffset(blockOffset) {}

    double* loadScalarData();

    double* loadBlock(const std::string& path, bool is2D = false) const;

    double* normalizedFromComponents(const std::array<double*, 3> velocity,
                                     const std::array<double*, 3> average, const double* rms) const;

    double* loadIsotrop() const;
    double* loadDuct();

    void getRmsTypeFromFilename(const std::string& rmsName);

    inline static bool fileExists(const std::string& name) {
        if (FILE* file = fopen(name.c_str(), "r")) {
            fclose(file);
            return true;
        } else {
            return false;
        }
    }

    static void setSettings(InputMode mode, vec3i totalSize, std::string directory,
                            std::string rmsName, ind timeStep) {
        Mode = mode;
        TotalSizeFile = totalSize;
        RmsFilename = rmsName;
        Directory = directory;
        TimeStep = timeStep;
        AvgValue = -1;
        RmsValue = -1;
    }

    static void setSettings(InputMode mode, vec3i totalSize, std::string filename, ind timeStep,
                            float avgValue, float rmsValue) {
        Mode = mode;
        TotalSizeFile = totalSize;
        Directory = filename;
        TimeStep = timeStep;
        AvgValue = avgValue;
        RmsValue = rmsValue;
    }

    static void setTimeStep(ind timeStep) { TimeStep = timeStep; }

private:
    typedef double (*ScalarFunc)(const std::array<double, 3>&, const std::array<double, 3>&);
    static const std::unordered_map<std::string, ScalarFunc> ScalarVariants;

    vec3i BlockSize, BlockOffset;
    ScalarFunc RmsFunction;

    static InputMode Mode;
    static vec3i TotalSizeFile;
    static std::string RmsFilename;
    static std::string Directory;
    static ind TimeStep;
    static float AvgValue;
    static float RmsValue;
};

}  // namespace perc