#include "percolationloader.h"
#include <iostream>
#include <fstream>
#include <cmath>

namespace perc {

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
    vec3i totalSize = TotalSize;
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

double* PercolationLoader::normalizedFromComponents(const std::array<double*, 3> velocity,
                                                    const std::array<double*, 3> average,
                                                    const double* rms) {
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

double* PercolationLoader::loadScalarData(ind timeSlice, const std::string& pathVelocity,
                                          const std::string& pathAverage,
                                          const std::string& pathRms,
                                          const std::string& pathVertex) {
    ind numElements = BlockSize.prod();

    //    PerformanceTimer Timer;
    ind numZeros = 4 - (ind)std::log10(timeSlice);
    std::string zeros = std::string(numZeros, '0');
    zeros = zeros;

    std::array<double*, 3> dataBuffer = {
        loadBlock(pathVelocity + zeros + std::to_string(timeSlice) + ".vx"),
        loadBlock(pathVelocity + zeros + std::to_string(timeSlice) + ".vy"),
        loadBlock(pathVelocity + zeros + std::to_string(timeSlice) + ".vx")};

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

    double* scalar = normalizedFromComponents(dataBuffer, avgBuffer, rmsBuffer);

    //    std::cout << "\t\tGrid creation and normalization took " << Timer.ElapsedTime() <<
    //    "seconds.";

    // Get rid of some memory
    for (int i = 0; i < 3; ++i) {
        delete[] dataBuffer[i];
        delete[] avgBuffer[i];
    }

    // Assume uniform volume for now

    return scalar;
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
