#include "percolationloader.h"
#include <iostream>
#include <fstream>
#include <cmath>

namespace perc {

// Static variables
ind PercolationLoader::TimeStep = 1;
std::string PercolationLoader::Directory = "";
std::string PercolationLoader::RmsFilename = "uv_000";
vec3i PercolationLoader::TotalSizeFile = vec3i(193, 194, 1000);
float PercolationLoader::AvgValue = 0.0f;
float PercolationLoader::RmsValue = 1.0f;
InputMode PercolationLoader::Mode = InputMode::INVALID;

std::mt19937 PercolationLoader::RandomGenerator;
std::uniform_real_distribution<> PercolationLoader::RandomDistribution;

const std::unordered_map<std::string, PercolationLoader::ScalarFunc>
    PercolationLoader::ScalarVariants(
        {{"uv", [](const std::array<double, 3>& data,
                   const std::array<double, 3>&) { return data[0] * data[1]; }},
         {"uw", [](const std::array<double, 3>& data,
                   const std::array<double, 3>&) { return data[0] * data[2]; }},
         {"vw", [](const std::array<double, 3>& data,
                   const std::array<double, 3>&) { return data[1] * data[2]; }},

         {"v2w2",
          [](const std::array<double, 3>& data, const std::array<double, 3>&) {
              return data[0] * data[0] - data[1] * data[1];
          }},
         {"K",
          [](const std::array<double, 3>&, const std::array<double, 3>& rawData) {
              return 0.5 * (rawData[0] * rawData[0] + rawData[1] * rawData[1]);
          }},  // No average!!! But divide by the "rms" file.
         {"k", [](const std::array<double, 3>& data, const std::array<double, 3>&) {
              return 0.5 * (data[0] * data[0] + data[1] * data[1] + data[2] * data[2]);
          }}});

double* PercolationLoader::loadBlock(const std::string& path, bool is2D) const {

    vec3i blockSize = BlockSize;
    vec3i totalSize = TotalSizeFile;
    vec3i blockOffset = BlockOffset;

    if (is2D) {
        blockSize.z = 1;
        totalSize.z = 1;
        blockOffset.z = 0;
    }

    ind totalN = totalSize.prod();
    ind blockN = blockSize.prod();

    // Load file.
    if (!fileExists(path)) {
        std::cerr << "File \"" << path << "\" does not exist.\n";
        return nullptr;
    }
    std::ifstream file(path, std::ios::binary);

    // Get length of file.
    ind fileLength;
    file.seekg(0, std::ios::end);
    fileLength = file.tellg();
    file.seekg(0, std::ios::beg);

    // Find out header size.
    ind headerSize = 0;
    if (fileLength != totalN * sizeof(double)) {
        ind expectedNumber = totalN * sizeof(double);
        ind header = 0;
        ind footer = 0;
        headerSize = 0;
        ind footerSize;
        for (footerSize = 4; footerSize <= 8; footerSize += 4) {

            file.seekg(fileLength - footerSize, std::ios::beg);
            file.read((char*)&footer, footerSize);
            if (footer == expectedNumber) {
                headerSize = fileLength - (totalN * sizeof(double) + footerSize);
                if (headerSize > 0) {
                    file.seekg(0, std::ios::beg);
                    file.read((char*)&header, headerSize);

                    if (header != totalN * sizeof(double)) continue;
                }
                break;
            }
            headerSize = 0;
        }

        if (!headerSize) {
            file.seekg(0, std::ios::beg);
            file.read((char*)&header, headerSize);

            if (header != totalN * sizeof(double)) {
                std::cout << "Could not determine header size.\n";
                std::cout << "\tHeader size: " << headerSize << " - " << header << std::endl;
                std::cout << "\tFooter size: " << footerSize << " - " << footer << std::endl;
                return nullptr;
            }
        }
    }

    // Copy data into buffer.
    char* buffer = new char[blockN * sizeof(double)];
    char* bufferIt = buffer;
    for (ind z = 0; z < blockSize.z; ++z) {
        for (ind y = 0; y < blockSize.y; ++y) {
            vec3i locIdx(0, y, z);
            vec3i globIdx = locIdx + blockOffset;

            file.seekg(globIdx.toIndexOfTotal(totalSize) * sizeof(double) + headerSize,
                       std::ios::beg);
            file.read(bufferIt, blockSize.x * sizeof(double));
            bufferIt += blockSize.x * sizeof(double);
        }
    }

    // Check whether reading was successful.
    if (!file) {
        std::cerr << "File \"" << path << "\" could not be read.\n";
        delete[] buffer;
        file.close();
        return nullptr;
    }
    file.close();

    return reinterpret_cast<double*>(buffer);
}

float* PercolationLoader::loadBlockFloat(const std::string& path, bool is2D) const {

    vec3i blockSize = BlockSize;
    vec3i totalSize = TotalSizeFile;
    vec3i blockOffset = BlockOffset;

    if (is2D) {
        blockSize.z = 1;
        totalSize.z = 1;
        blockOffset.z = 0;
    }

    ind totalN = totalSize.prod();
    ind blockN = blockSize.prod();

    // Load file.
    if (!fileExists(path)) {
        std::cerr << "File \"" << path << "\" does not exist.\n";
        return nullptr;
    }
    std::ifstream file(path, std::ios::binary);

    // Get length of file.
    ind fileLength;
    file.seekg(0, std::ios::end);
    fileLength = file.tellg();
    file.seekg(0, std::ios::beg);

    // Find out header size.
    ind headerSize = 0;
    if (fileLength != totalN * sizeof(float)) {
        // TODO: This is only due to a bug in filewriting, the header gives number of elements in
        // double!!
        ind expectedNumber = totalN * sizeof(double);
        // std::cout << "Filelength: " << fileLength << ", Expected: " << expectedNumber << ".\n";
        ind header = 0;
        ind footer = 0;
        headerSize = 0;
        ind footerSize;
        for (footerSize = 4; footerSize <= 8; footerSize += 4) {

            file.seekg(fileLength - footerSize, std::ios::beg);
            file.read((char*)&footer, footerSize);
            // std::cout << footer << ".\n";
            if (footer == expectedNumber) {
                headerSize = fileLength - (totalN * sizeof(float) + footerSize);
                if (headerSize > 0) {
                    file.seekg(0, std::ios::beg);
                    file.read((char*)&header, headerSize);

                    if (header != totalN * expectedNumber) continue;
                }
                break;
            }
            headerSize = 0;
        }

        if (!headerSize) {
            file.seekg(0, std::ios::beg);
            file.read((char*)&header, headerSize);

            if (header != expectedNumber) {
                std::cout << "Could not determine header size.\n";
                std::cout << "\tHeader size: " << headerSize << " - " << header << std::endl;
                std::cout << "\tFooter size: " << footerSize << " - " << footer << std::endl;
                return nullptr;
            }
        }
    }

    // Copy data into buffer.
    char* buffer = new char[blockN * sizeof(float)];
    char* bufferIt = buffer;
    for (ind z = 0; z < blockSize.z; ++z) {
        for (ind y = 0; y < blockSize.y; ++y) {
            vec3i locIdx(0, y, z);
            vec3i globIdx = locIdx + blockOffset;

            file.seekg(globIdx.toIndexOfTotal(totalSize) * sizeof(float) + headerSize,
                       std::ios::beg);
            file.read(bufferIt, blockSize.x * sizeof(float));
            bufferIt += blockSize.x * sizeof(float);
        }
    }

    // Check whether reading was successful.
    if (!file) {
        std::cerr << "File \"" << path << "\" could not be read.\n";
        delete[] buffer;
        file.close();
        return nullptr;
    }
    file.close();

    return reinterpret_cast<float*>(buffer);
}

double* PercolationLoader::normalizedFromComponents(const std::array<double*, 3> velocity,
                                                    const std::array<double*, 3> average,
                                                    const double* rms) const {
    if (!RmsFunction) {
        std::cerr << "No rms function defined.";
        return nullptr;
    }

    ind numElements = BlockSize.prod();
    ind numStatElements = BlockSize.x * BlockSize.y;

    double* scalar = new double[numElements];

#pragma omp parallel for
    for (ind xyz = 0; xyz < numElements; ++xyz) {
        ind xy = xyz % numStatElements;
        // Raw data pointer to write to.
        std::array<double, 3> components, normComponents;
        for (int n = 0; n < 3; ++n) {
            // Compute the percolation analysis scalar value.
            components[n] = velocity[n][xyz];
            normComponents[n] = fabs(velocity[n][xyz] - average[n][xy]);
            normComponents[n] = std::isfinite(normComponents[n]) ? normComponents[n] : 0;
        }
        scalar[xyz] = rms[xy] > 0 ? RmsFunction(normComponents, components) / rms[xy] : 0;
    }

    return scalar;
}

double* PercolationLoader::loadIsotrop() const {
    ind numElements = TotalSizeFile.prod();
    double* scalar;

    if (!fileExists(Directory)) {
        std::cerr << "File \"" << Directory << "\" does not exist.\n";
        return nullptr;
    }
    std::ifstream file(Directory, std::ios::binary);

    // Get length of file.
    ind fileLength;
    file.seekg(0, std::ios::end);
    fileLength = file.tellg();
    file.seekg(0, std::ios::beg);
    if (fileLength != numElements * sizeof(float)) {
        std::cerr << "File \"" << Directory << "\" does not have right size, expected "
                  << numElements << " floats.";
        return nullptr;
    }

    // Correct length. Just load all into file.
    ind blockNum = BlockSize.prod();

    char* buffer = new char[blockNum * sizeof(float)];
    char* bufferIt = buffer;
    for (ind z = 0; z < BlockSize.z; ++z) {
        for (ind y = 0; y < BlockSize.y; ++y) {
            vec3i locIdx(0, y, z);
            vec3i globIdx = locIdx + BlockOffset;

            file.seekg(globIdx.toIndexOfTotal(TotalSizeFile) * sizeof(float), std::ios::beg);
            file.read(bufferIt, BlockSize.x * sizeof(float));
            bufferIt += BlockSize.x * sizeof(float);
        }
    }

    // Reading done. Convert to normalized double.
    scalar = new double[blockNum];
    float* bufferData = reinterpret_cast<float*>(buffer);
    for (ind i = 0; i < blockNum; ++i) {
        scalar[i] = fabs(bufferData[i] - AvgValue) / RmsValue;
    }
    delete[] buffer;

    return scalar;
}

double* PercolationLoader::loadDuctWithoutNormalization() {
    ind numElements = BlockSize.prod();
    double* scalar;

    const std::string& pathVelocity = Directory + "/VELOCITY/";

    ind numZeros = 4 - (ind)std::log10(TimeStep);
    std::string zeros = std::string(numZeros, '0');
    zeros = zeros;

    std::array<double*, 3> dataBuffer = {
        loadBlock(pathVelocity + zeros + std::to_string(TimeStep) + ".vx"),
        loadBlock(pathVelocity + zeros + std::to_string(TimeStep) + ".vy"),
        loadBlock(pathVelocity + zeros + std::to_string(TimeStep) + ".vx")};

    const std::string& pathRms = Directory + "/ZEXPORT_STAT_wall_correction/" + RmsFilename;
    size_t nameBeginPos = pathRms.find_last_of('/');
    std::string filenameRms = pathRms.substr(nameBeginPos + 1, pathRms.length() - 2 - nameBeginPos);

    getRmsTypeFromFilename(filenameRms);

    if (!dataBuffer[0] || !dataBuffer[1] || !dataBuffer[2]) {
        for (int i = 0; i < 3; ++i) {
            delete[] dataBuffer[i];
        }
        std::cerr << "Could not load all needed files.\n";
        return nullptr;
    }

    scalar = new double[numElements];

    if (!RmsFunction) {
        std::cerr << "No rms function defined.";
        return nullptr;
    }

#pragma omp parallel for
    for (ind xyz = 0; xyz < numElements; ++xyz) {
        // Raw data pointer to write to.
        std::array<double, 3> components, normComponents;
        for (int n = 0; n < 3; ++n) {
            // Compute the percolation analysis scalar value.
            components[n] = dataBuffer[n][xyz];
            normComponents[n] = dataBuffer[n][xyz];
            normComponents[n] = std::isfinite(normComponents[n]) ? normComponents[n] : 0;
        }
        scalar[xyz] = RmsFunction(normComponents, components);
    }

    // Get rid of some memory
    for (int i = 0; i < 3; ++i) {
        delete[] dataBuffer[i];
    }

    return scalar;
}

double* PercolationLoader::loadDuct() {
    ind numElements = BlockSize.prod();
    double* scalar;

    const std::string& pathVelocity = Directory + "/VELOCITY/";
    const std::string& pathAverage = Directory + "/STAT/";
    const std::string& pathRms = Directory + "/ZEXPORT_STAT_wall_correction/" + RmsFilename;
    const std::string& pathVertex = Directory + "/STAT/";

    //    PerformanceTimer Timer;
    ind numZeros = 4 - (ind)std::log10(TimeStep);
    std::string zeros = std::string(numZeros, '0');
    zeros = zeros;

    std::array<double*, 3> dataBuffer = {
        loadBlock(pathVelocity + zeros + std::to_string(TimeStep) + ".vx"),
        loadBlock(pathVelocity + zeros + std::to_string(TimeStep) + ".vy"),
        loadBlock(pathVelocity + zeros + std::to_string(TimeStep) + ".vx")};

    std::array<double*, 3> avgBuffer = {loadBlock(pathAverage + "average_vx", true),
                                        loadBlock(pathAverage + "average_vy", true),
                                        loadBlock(pathAverage + "average_vz", true)};

    double* rmsBuffer = loadBlock(pathRms, true);

    // Cut out file name, match to rms types
    size_t nameBeginPos = pathRms.find_last_of('/');
    std::string filenameRms = pathRms.substr(nameBeginPos + 1, pathRms.length() - 2 - nameBeginPos);

    getRmsTypeFromFilename(filenameRms);

    if (!dataBuffer[0] || !dataBuffer[1] || !dataBuffer[2] || !avgBuffer[0] || !avgBuffer[1] ||
        !avgBuffer[2] || !rmsBuffer) {
        for (int i = 0; i < 3; ++i) {
            delete[] dataBuffer[i];
            delete[] avgBuffer[i];
        }
        delete[] rmsBuffer;
        std::cerr << "Could not load all needed files.\n";
        return nullptr;
    }

    //    std::cout << "Raw data loading took " << Timer.ElapsedTime() << " seconds.";

    //   Timer.Reset();

    scalar = normalizedFromComponents(dataBuffer, avgBuffer, rmsBuffer);

    //    std::cout << "\t\tGrid creation and normalization took " << Timer.ElapsedTime() <<
    //    "seconds.";

    // Get rid of some memory
    for (int i = 0; i < 3; ++i) {
        delete[] dataBuffer[i];
        delete[] avgBuffer[i];
    }
    return scalar;
}

double* PercolationLoader::loadWing() {
    ind numElements = BlockSize.prod();

    const std::string& pathVelocity = Directory + "/Fields/";
    const std::string& pathAverage = Directory + "/Stat/";
    const std::string& pathRms = Directory + "/Stat/" + RmsFilename;
    const std::string& pathVertex = Directory + "/Stat/";

    //    PerformanceTimer Timer;
    ind numZeros = 4 - (ind)std::log10(TimeStep);
    std::string zeros = std::string(numZeros, '0');
    zeros = zeros;

    std::array<float*, 3> dataBuffer = {
        loadBlockFloat(pathVelocity + zeros + std::to_string(TimeStep) + ".vx"),
        loadBlockFloat(pathVelocity + zeros + std::to_string(TimeStep) + ".vy"),
        loadBlockFloat(pathVelocity + zeros + std::to_string(TimeStep) + ".vx")};

    std::array<float*, 3> avgBuffer = {loadBlockFloat(pathAverage + "average_vx", true),
                                       loadBlockFloat(pathAverage + "average_vy", true),
                                       loadBlockFloat(pathAverage + "average_vz", true)};

    float* rmsBuffer = loadBlockFloat(pathRms, true);

    // For now: Hard-code uv
    getRmsTypeFromFilename("uv");

    if (!dataBuffer[0] || !dataBuffer[1] || !dataBuffer[2] || !avgBuffer[0] || !avgBuffer[1] ||
        !avgBuffer[2] || !rmsBuffer) {
        for (int i = 0; i < 3; ++i) {
            delete[] dataBuffer[i];
            delete[] avgBuffer[i];
        }
        delete[] rmsBuffer;
        std::cerr << "Could not load all needed files.\n";
        return nullptr;
    }

    //    std::cout << "Raw data loading took " << Timer.ElapsedTime() << " seconds.";

    //   Timer.Reset();

    // Copy from normalizedFromComponents
    if (!RmsFunction) {
        std::cerr << "No rms function defined.";
        return nullptr;
    }

    ind numStatElements = BlockSize.x * BlockSize.y;

    double* scalar = new double[numElements];

#pragma omp parallel for
    for (ind xyz = 0; xyz < numElements; ++xyz) {
        ind xy = xyz % numStatElements;
        // Raw data pointer to write to.
        std::array<double, 3> components, normComponents;
        for (int n = 0; n < 3; ++n) {
            // Compute the percolation analysis scalar value.
            components[n] = static_cast<double>(dataBuffer[n][xyz]);
            normComponents[n] = fabs(dataBuffer[n][xyz] - static_cast<double>(avgBuffer[n][xy]));
            normComponents[n] = std::isfinite(normComponents[n]) ? normComponents[n] : 0;
        }
        scalar[xyz] =
            rmsBuffer[xy] > 0 ? RmsFunction(normComponents, components) / rmsBuffer[xy] : 0;
    }

    //    std::cout << "\t\tGrid creation and normalization took " << Timer.ElapsedTime() <<
    //    "seconds.";

    // Get rid of some memory
    for (int i = 0; i < 3; ++i) {
        delete[] dataBuffer[i];
        delete[] avgBuffer[i];
    }
    return scalar;
}

double* PercolationLoader::loadRandom() const {
    ind blockNum = BlockSize.prod();
    double* scalar = new double[blockNum];

    for (ind i = 0; i < blockNum; ++i) {
        scalar[i] = RandomDistribution(RandomGenerator);
    }
    return scalar;
}

double* PercolationLoader::loadScalarData() {
    switch (Mode) {
        case InputMode::COMBINED_VELOCITY_AVG_RMS_FILE:
            return loadDuct();
        case InputMode::COMBINED_VELOCITY_AVG_RMS_FILE_WING:
            return loadWing();
        case InputMode::VELOCITY_FILE:
            return loadDuctWithoutNormalization();
        case InputMode::COMBINED_VELOCITY_AVG_RMS_VALUE:
            return loadIsotrop();
        case InputMode::RANDOM_UNIFORM:
            return loadRandom();
        default:
            return nullptr;
    }
}

void PercolationLoader::getRmsTypeFromFilename(const std::string& rmsName) {
    // Get type of scalar. Sizeing is important!
    size_t compLength = rmsName.find('_', 0);
    std::string component = rmsName.substr(0, compLength);

    // Find the corresponding function.
    bool found = false;

    for (auto& variant : ScalarVariants)
        if (!component.compare(variant.first)) {
            RmsFunction = variant.second;
            found = true;
            break;
        }

    // Could not find a match.
    if (!found) {
        std::cerr << "Unknown scalar component name \"" << component << '\"';
        RmsFunction = nullptr;
        return;
    }
}

}  // namespace perc
