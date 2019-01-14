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
std::pair<vec3i, vec3i> buildRedBlocks(const vec3i& parentBlockSize, const vec3i& parentBlockOffset,
                                       const vec3i& totalSize,
                                       std::vector<UnionFindSubBlock<BlockType>>& newBlocks,
                                       ID*& newPointers, UnionFindBlock& parent,
                                       std::function<BlockType()> blockConstructor) {
    newBlocks.reserve(6);

    // Find out white block size, offset and size to allocate for pointers.
    ind pointerBlockSize = 0;
    vec3i whiteBlockOffset = parentBlockOffset;
    vec3i whiteBlockSize = parentBlockSize;

    // Gather white size in all dimensions.
    for (ind dim = 0; dim < 3; ++dim) {

        // Low side.
        if (whiteBlockOffset[dim] > 0) whiteBlockOffset[dim]++;

        // High side.
        if (whiteBlockOffset[dim] + whiteBlockSize[dim] < totalSize[dim]) whiteBlockSize[dim]--;
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
            UnionFindSubBlock<BlockType> redBlock = UnionFindSubBlock<BlockType>(
                blockSize, blockOffset, totalSize, parent, blockConstructor(), pointerIterator);
            newBlocks.push_back(redBlock);
            pointerIterator += blockSize.prod();
        }

        // Upper slice.
        if (whiteBlockOffset[dim] + whiteBlockSize[dim] < totalSize[dim]) {

            vec3i blockOffset = whiteBlockOffset;
            blockOffset[dim] += whiteBlockSize[dim] - 1;

            // Create block.
            UnionFindSubBlock<BlockType> redBlock = UnionFindSubBlock<BlockType>(
                blockSize, blockOffset, totalSize, parent, blockConstructor(), pointerIterator);
            newBlocks.push_back(redBlock);
            pointerIterator += blockSize.prod();
        }
    }

    return std::make_pair(whiteBlockSize, whiteBlockOffset);
}
}  // namespace interfaceblockbuilder

}  // namespace perc