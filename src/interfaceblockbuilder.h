#pragma once
#include <vector>
#include <functional>
#include "handles.h"
#include "unionfindblock.h"
#include "unionfindsubblock.h"

namespace perc {
namespace interfaceblockbuilder {
// Create the red blocks as red or grey subblocks.
// @returns size of white block within.
template <typename BlockType>
void buildRedBlocks(const vec3i& parentBlockSize, const vec3i& parentBlockOffset,
                    const vec3i& totalSize, std::vector<UnionFindSubBlock<BlockType>>& newBlocks,
                    ID*& newPointers, ind& pointerBlockSize, vec3i& whiteBlockSize,
                    vec3i& whiteBlockOffset, UnionFindBlock& parent,
                    std::function<BlockType()> blockConstructor) {
    newBlocks.reserve(6);

    // Find out white block size, offset and size to allocate for pointers.
    pointerBlockSize = 0;
    whiteBlockOffset = parentBlockOffset;
    whiteBlockSize = parentBlockSize;

    // Gather white size in all dimensions.
    for (ind dim = 0; dim < 3; ++dim) {

        // Low side.
        if (whiteBlockOffset[dim] > 0) {
            whiteBlockOffset[dim]++;
            whiteBlockSize[dim]--;
        }

        // High side.
        if (whiteBlockOffset[dim] + whiteBlockSize[dim] < totalSize[dim]) whiteBlockSize[dim] -= 2;
    }

    // Gather sizes in all dimensions.
    for (ind dim = 0; dim < 3; ++dim) {
        vec3i blockSize = whiteBlockSize;
        blockSize[dim] = 1;

        // Low side.
        if (whiteBlockOffset[dim] > 0) pointerBlockSize += blockSize.prod();

        // High side.
        if (whiteBlockOffset[dim] + whiteBlockSize[dim] < totalSize[dim])
            pointerBlockSize += blockSize.prod();
    }

    // Allocate memory.
    newPointers = new ID[pointerBlockSize];
    ID* pointerIterator = newPointers;

    // Build blocks.
    for (ind dim = 0; dim < 3; ++dim) {
        vec3i blockSize = whiteBlockSize;
        blockSize[dim] = 1;

        // Lower slice.
        if (whiteBlockOffset[dim] > 0) {

            vec3i blockOffset = whiteBlockOffset;
            blockOffset[dim] -= 1;

            // Create block.
            newBlocks.emplace_back(blockSize, blockOffset, totalSize, parent, blockConstructor(),
                                   pointerIterator);
            pointerIterator += blockSize.prod();
        }

        // Upper slice.
        if (whiteBlockOffset[dim] + whiteBlockSize[dim] < totalSize[dim]) {

            vec3i blockOffset = whiteBlockOffset;
            // WhiteBlockSize has already been adapted for the red block
            blockOffset[dim] += whiteBlockSize[dim];

            // Create block.
            newBlocks.emplace_back(blockSize, blockOffset, totalSize, parent, blockConstructor(),
                                   pointerIterator);
            pointerIterator += blockSize.prod();
        }
    }
}

// Create a green blocks as green or grey subblocks.
template <typename BlockType>
UnionFindSubBlock<BlockType> buildGreenBlock(const vec3i& direction, const vec3i& whiteBlockSize,
                                             const vec3i& whiteBlockOffset, const vec3i& totalSize,
                                             ID*& memoryPointer, UnionFindBlock& parent,
                                             std::function<BlockType()> blockConstructor) {

    vec3i blockOffset, blockSize;
    static const ind ALL_ZERO = -1;
    static const ind SEVERAL_NONZERO = -2;

    ind nonZeroDim = ALL_ZERO;
    for (ind dim = 0; dim < 3; ++dim) {
        if (direction[dim] != 0) {
            if (nonZeroDim == ALL_ZERO)
                nonZeroDim = dim;
            else {
                nonZeroDim = SEVERAL_NONZERO;
                break;
            }
        }
    }
    assert(nonZeroDim != ALL_ZERO && "Null vector input direction.");

    // Corner block?
    if (nonZeroDim == SEVERAL_NONZERO) {

        // Little 3x3x3 corner block or 3x3xS edge.
        for (ind dim = 0; dim < 3; ++dim) {
            switch (direction[dim]) {
                case 1:
                    blockOffset[dim] = whiteBlockOffset[dim] + whiteBlockSize[dim];
                    blockSize[dim] = 3;
                    break;

                case -1:
                    blockOffset[dim] = whiteBlockOffset[dim] - 3;
                    blockSize[dim] = 3;
                    break;

                case 0:
                    blockOffset[dim] = whiteBlockOffset[dim];
                    blockSize[dim] = whiteBlockSize[dim];
                    break;

                default:
                    assert(false && "Only allowing direction within {-1, 0, 1}³.");
            }
        }

        // Side or edge slice.
    } else {
        blockOffset = whiteBlockOffset;
        blockSize = whiteBlockSize;

        blockSize[nonZeroDim] = 1;
        switch (direction[nonZeroDim]) {
            case -1:
                assert(whiteBlockOffset[nonZeroDim] >= 4 && "No green block fits here.");
                blockOffset[nonZeroDim] = whiteBlockOffset[nonZeroDim] - 2;
                break;

            case 1:
                assert(whiteBlockOffset[nonZeroDim] + whiteBlockSize[nonZeroDim] <
                           totalSize[nonZeroDim] - 3 &&
                       "No green block here.");
                blockOffset[nonZeroDim] =
                    whiteBlockOffset[nonZeroDim] + whiteBlockSize[nonZeroDim] + 1;
                break;

            default:
                assert(false && "Only allowing direction within {-1, 0, 1}³.");
        }
    }

    auto block = UnionFindSubBlock<BlockType>(blockSize, blockOffset, totalSize, parent,
                                              blockConstructor(), memoryPointer);

    memoryPointer += blockSize.prod();
    return block;
}

}  // namespace interfaceblockbuilder

}  // namespace perc