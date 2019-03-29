#include "localblock.h"
#include "interfaceblockbuilder.h"
#include "mpicommuncation.h"

namespace perc {

LocalBlock::LocalBlock(const vec3i& blockSize, const vec3i& blockOffset, const vec3i& totalSize,
                       const int rank)
    : UnionFindBlock(totalSize)
    , Rank(rank)
    , RefPLOGs(new RefPLOGtype(10000, &ClusterID::hash))
    , LOLs(new ClusterList(LOCAL_LIST))
    , LOGs(new ClusterListRecordingSingle(GLOBAL_LIST))
    , LOGSubBlocks() {

    auto constructor = [this]() -> RedProcessor { return RedProcessor(*LOLs, *LOGs, *RefPLOGs); };
    vec3i whiteBlockSize, whiteBlockOffset;
    interfaceblockbuilder::buildRedBlocks<RedProcessor>(
        blockSize, blockOffset, totalSize, LOGSubBlocks, MemoryLOG, MemoryLOGSize, whiteBlockSize,
        whiteBlockOffset, *this, constructor);

    LOLSubBlock =
        new UnionFindSubBlock<WhiteProcessor>(whiteBlockSize, whiteBlockOffset, totalSize, *this,
                                              WhiteProcessor(*LOLs, *LOGs, *RefPLOGs), nullptr);

    std::vector<vec3i> directions;
    vec3i potentialMax = blockOffset + blockSize;

    for (ind dim = 0; dim < 3; ++dim) {
        for (ind sign = -1; sign <= 1; sign += 2) {
            vec3i dir(0);
            // Green side in the lower or upper limits
            if ((sign == -1 && blockOffset[dim] > 0) ||
                (sign == 1 && potentialMax[dim] < totalSize[dim])) {
                // Direction the green side lies at.
                dir[dim] = sign;

                // Add to possible direction combinations.
                for (ind d = directions.size() - 1; d >= 0; --d)
                    if (directions[d][dim] == 0) {
                        directions.push_back(directions[d] + dir);
                    }
                directions.push_back(dir);
            }
        }
    }

    // Sort directions, such that the blocks will have the same order as in the global block
    // first come all that have -1 in z, then all that have -1 in y, then in x.
    // If both are -1 or not -1 on all dimension, order according to linear index
    std::sort(directions.begin(), directions.end(), [](auto& a, auto& b) {
        for (ind dim = 2; dim >= 0; --dim) {
            if (a[dim] == -1 && b[dim] != -1) return true;
            if (a[dim] != -1 && b[dim] == -1) return false;
        }
        return (a + 1).toIndexOfTotal(vec3i(3)) < (b + 1).toIndexOfTotal(vec3i(3));
    });

    // Setup global block vector.
    GOGSubBlocks.reserve(directions.size());
    ID* dummyMemory;

    for (auto& dir : directions) {
        dummyMemory = nullptr;

        // Make a green block.
        GOGSubBlocks.push_back(interfaceblockbuilder::buildGreenBlock<GrayProcessor>(
            dir, whiteBlockSize, whiteBlockOffset, totalSize, dummyMemory, *this,
            []() { return GrayProcessor(); }));
    }

#ifndef SINGLENODE
    Rank = 0;
#endif

    /*
        std::cout << Rank << ": ";
        for (auto& gogBlock : GOGSubBlocks) {
            std::cout << gogBlock.blockOffset() << "\t";
        }
        std::cout << std::endl;
        */
}

LocalBlock::LocalBlock(const vec3i& totalSize)
    : UnionFindBlock(totalSize)
    , Rank(1)
    , RefPLOGs(new RefPLOGtype(10000, &ClusterID::hash))
    , LOLs(new ClusterList(LOCAL_LIST))
    , LOGs(new ClusterListRecordingSingle(GLOBAL_LIST)) {}

LocalBlock::LocalBlock(LocalBlock&& other)
    : UnionFindBlock(other.TotalSize)
    , Rank(std::move(other.Rank))
    , LOGSubBlocks(std::move(other.LOGSubBlocks))
    , GOGSubBlocks(std::move(other.GOGSubBlocks)) {
    // Take over log memory.
    MemoryLOG = other.MemoryLOG;
    other.MemoryLOG = nullptr;
    MemoryLOGSize = other.MemoryLOGSize;
    other.MemoryLOGSize = 0;

    // Take over cluster lists.
    LOLs = other.LOLs;
    other.LOLs = nullptr;
    LOGs = other.LOGs;
    other.LOGs = nullptr;
    RefPLOGs = other.RefPLOGs;
    other.RefPLOGs = nullptr;

    // Take over white block on heap.
    LOLSubBlock = other.LOLSubBlock;
    other.LOLSubBlock = nullptr;

    // Setting parent pointers anew, should always point to this.
    LOLSubBlock->Parent = this;
    for (auto& log : LOGSubBlocks) log.Parent = this;
    for (auto& gog : GOGSubBlocks) gog.Parent = this;
}

// LocalBlock LocalBlock::operator=(LocalBlock&& other) {
//     TotalSize = other.TotalSize;
//     LOLSubBlock(std::move(other.LOLSubBlock)) LOGSubBlocks(std::move(other.LOGSubBlocks));
//     GOGSubBlocks(std::move(other.GOGSubBlocks)) LOLs(std::move(other.LOLs));
//     LOGs(std::move(other.LOGs)) RefPLOGs(std::move(other.RefPLOGs));
//     CommData(std::move(other.CommData));
//     MemoryLOG = other.MemoryLOG;
//     other.MemoryLOG = nullptr;
//     MemoryLOGSize = other.MemoryLOGSize;
//     other.MemoryLOGSize = 0;
// }

LocalBlock::~LocalBlock() {
    delete[] MemoryLOG;
    delete LOLSubBlock;
    delete LOLs;
    delete LOGs;
    delete RefPLOGs;
}

LocalBlock* LocalBlock::makeGroundtruth(const vec3i& blockSize, const vec3i& blockOffset,
                                        const vec3i& totalSize) {
    LocalBlock* block = new LocalBlock(totalSize);
    block->LOLSubBlock = new UnionFindSubBlock<WhiteProcessor>(
        blockSize, blockOffset, totalSize, *block,
        WhiteProcessor(*block->LOLs, *block->LOGs, *block->RefPLOGs));
    block->LOLSubBlock->loadData();
    block->sortData(false);
    // Id Block for Red is empty
    block->MemoryLOG = nullptr;
    block->MemoryLOGSize = 0;

    return block;
}

LocalBlock* LocalBlock::makeWhiteRedTest(const vec3i& blockSize, const vec3i& blockOffset,
                                         const vec3i& totalSize) {
    LocalBlock* block = new LocalBlock(totalSize);

    vec3i min = blockOffset;
    vec3i max = blockOffset + blockSize;

    // Create 3 red slices: one left, two right.
    vec3i sliceSize = vec3i(blockSize.x, blockSize.y, 10);
    min[2] += sliceSize.z;
    max[2] -= sliceSize.z * 2;

    block->LOLSubBlock = new UnionFindSubBlock<WhiteProcessor>(
        max - min, min, totalSize, *block,
        WhiteProcessor(*block->LOLs, *block->LOGs, *block->RefPLOGs));
    block->LOLSubBlock->loadData();
    block->LOLSubBlock->sortData(false);

    block->MemoryLOGSize = sliceSize.prod() * 3;
    block->MemoryLOG = new ID[block->MemoryLOGSize];

    block->LOGSubBlocks.reserve(3);
    // Left slice.
    block->LOGSubBlocks.emplace_back(sliceSize, blockOffset, totalSize, *block,
                                     RedProcessor(*block->LOLs, *block->LOGs, *block->RefPLOGs),
                                     block->MemoryLOG);

    // Two slices right.
    block->LOGSubBlocks.emplace_back(sliceSize, vec3i(blockOffset.x, blockOffset.y, max.z),
                                     totalSize, *block,
                                     RedProcessor(*block->LOLs, *block->LOGs, *block->RefPLOGs),
                                     block->MemoryLOG + sliceSize.prod());
    block->LOGSubBlocks.emplace_back(
        sliceSize, vec3i(blockOffset.x, blockOffset.y, max.z + sliceSize.z), totalSize, *block,
        RedProcessor(*block->LOLs, *block->LOGs, *block->RefPLOGs),
        block->MemoryLOG + 2 * sliceSize.prod());

    block->loadData();
    block->sortData(false);

    return block;
}

LocalBlock* LocalBlock::makeWhiteRedGreenTest(const vec3i& blockSize, const vec3i& blockOffset,
                                              const vec3i& totalSize) {
    LocalBlock* block = new LocalBlock(totalSize);

    vec3i min = blockOffset;
    vec3i max = blockOffset + blockSize;

    // Create 3 green slices: two left, one right.
    // Create 2 red slices: one left, one right,
    vec3i sliceSize = vec3i(blockSize.x, blockSize.y, 10);
    min[2] += sliceSize.z * 3;
    max[2] -= sliceSize.z * 2;

    block->LOLSubBlock = new UnionFindSubBlock<WhiteProcessor>(
        max - min, min, totalSize, *block,
        WhiteProcessor(*block->LOLs, *block->LOGs, *block->RefPLOGs));

    block->MemoryLOGSize = sliceSize.prod() * 2;
    block->MemoryLOG = new ID[block->MemoryLOGSize];

    // Red data.
    block->LOGSubBlocks.reserve(2);
    // Left slice.
    block->LOGSubBlocks.emplace_back(
        sliceSize, vec3i(blockOffset.x, blockOffset.y, min.z - sliceSize.z), totalSize, *block,
        RedProcessor(*block->LOLs, *block->LOGs, *block->RefPLOGs), block->MemoryLOG);

    // Right right.
    block->LOGSubBlocks.emplace_back(sliceSize, vec3i(blockOffset.x, blockOffset.y, max.z),
                                     totalSize, *block,
                                     RedProcessor(*block->LOLs, *block->LOGs, *block->RefPLOGs),
                                     block->MemoryLOG + sliceSize.prod());
    // Green data.
    block->GOGSubBlocks.reserve(2);

    // Left slice (only the right most one of the two left slices that the global block holds is
    // contained in this local block
    block->GOGSubBlocks.emplace_back(
        sliceSize, vec3i(blockOffset.x, blockOffset.y, blockOffset.z + sliceSize.z), totalSize,
        *block, GrayProcessor());

    // Right slice.
    block->GOGSubBlocks.emplace_back(sliceSize,
                                     vec3i(blockOffset.x, blockOffset.y, max.z + sliceSize.z),
                                     totalSize, *block, GrayProcessor());

    block->loadData();
    block->sortData(false);

    return block;
}

void LocalBlock::reset() {
    for (auto& green : GOGSubBlocks) green.reset();
    for (auto& red : LOGSubBlocks) red.reset();
    LOGs->reset();
    LOLs->reset();
}

void LocalBlock::loadData() {
    LOLSubBlock->loadData();
    for (auto& red : LOGSubBlocks) red.loadData();
    for (auto& green : GOGSubBlocks) green.loadData();
}

void LocalBlock::sortData(bool useBuckets) {
    LOLSubBlock->sortData(useBuckets);
    for (auto& red : LOGSubBlocks) red.sortData(useBuckets);
    for (auto& green : GOGSubBlocks) green.sortData(useBuckets);
}

void LocalBlock::doWatershed(const double minVal) {
    LOLSubBlock->doWatershed(minVal);
    for (auto& log : LOGSubBlocks) log.doWatershed(minVal);
}

ClusterID* LocalBlock::findClusterID(const vec3i& idx, vec3i& lastClusterID) {
    // One might want to do this more cleverly, especialy in the sheet tree.
    if (LOLSubBlock->contains(idx)) return LOLSubBlock->findClusterID(idx, lastClusterID);
    for (auto& log : LOGSubBlocks)
        if (log.contains(idx)) return log.findClusterID(idx, lastClusterID);

    for (auto& gog : GOGSubBlocks)
        if (gog.contains(idx)) return gog.findClusterID(idx, lastClusterID);

    assert(false && "Can not find block containing this idx.");
    return nullptr;
}

ID* LocalBlock::setID(const vec3i& idx, const ID& id) {
    // One might want to do this more cleverly, especialy in the sheet tree.
    ID* ptr = nullptr;
    if (LOLSubBlock->contains(idx))
        ptr = LOLSubBlock->PointerBlock.getPointer(idx);
    else {
        for (auto& log : LOGSubBlocks)
            if (log.contains(idx)) {
                ptr = log.PointerBlock.getPointer(idx);
                break;
            }
    }

    // No need to look into GOGBlocks, we should not be writing there anyways

    assert(ptr && "Can not find block containing this idx.");
    *ptr = id;

    return ptr;
}

double LocalBlock::getClusterVolume(ClusterID cluster) {
    return cluster.isGlobal() ? LOGs->getClusterVolume(cluster) : LOLs->getClusterVolume(cluster);
}

void LocalBlock::receiveData() {
    ind numNewLOGs;
    ind startOfLocalPlog;
    std::vector<ind> merges;
#ifdef COMMUNICATION
    // Reuse statuss in blocking receive
    MPI_Status status;

    // Tags are int
    int rankTag = Rank << MPICommunication::RANK_SHIFT;

    // Receive data common for all processes
#ifdef COLLECTIVES
    MPI_Bcast(&numNewLOGs, 1, MPI_IND, 0, MPI_COMM_WORLD);
    ind mergesSize;
    MPI_Bcast(&mergesSize, 1, MPI_IND, 0, MPI_COMM_WORLD);
    merges.resize(mergesSize);
    MPI_Bcast(merges.data(), mergesSize * sizeof(ind), MPI_BYTE, 0, MPI_COMM_WORLD);
#else
    MPICommunication::handleError(MPI_Recv(&numNewLOGs, 1, MPI_IND, 0,
                                           MPICommunication::NUMNEWCLUSTERS | rankTag,
                                           MPI_COMM_WORLD, &status));
    MPICommunication::handleError(MPICommunication::RecvVectorUknownSize(
        merges, 0, MPICommunication::MERGES | rankTag, MPI_COMM_WORLD, &status));
#endif  // COLLECTIVES

    MPICommunication::handleError(MPI_Recv(&startOfLocalPlog, 1, MPI_IND, 0,
                                           MPICommunication::STARTOFPLOG | rankTag, MPI_COMM_WORLD,
                                           &status));

    // Receive updated green blocks
    int counter = 0;
    for (int counter = 0; counter < GOGSubBlocks.size(); ++counter) {
        auto& gogBlock = GOGSubBlocks[counter];
        // Todo: Should this potentially be Non-Blocking?
        MPICommunication::handleError(MPI_Recv(
            gogBlock.PointerBlock.PointerBlock, gogBlock.blockSize().prod() * sizeof(ID), MPI_BYTE,
            0, MPICommunication::GREENPOINTERS | rankTag | counter, MPI_COMM_WORLD, &status));
    }

#else   // !COMMUNICATION
    numNewLOGs = CommData.PLOGs.size();
    startOfLocalPlog = 0;
    std::vector<std::vector<ClusterMerge>> mergesLocal = {LOGs->Merges};
    merges = ClusterMerge::mergeClusterAsList(ClusterMerge::mergeClustersFromLists(mergesLocal));
#endif  // COMMUNICATION

    // Add some incognito clusters of other compute nodes.
    LOGs->addClusters(startOfLocalPlog);

    // For each PLOG, add a new cluster and reference it.
    for (ClusterData& c : CommData.PLOGs) {
        vec3i cPos = vec3i::fromIndexOfTotal(c.Index.RawID, TotalSize);
        ClusterID newID = LOGs->addCluster(c.Volume);
        // Add representative and repointer PLOG (now LOG)

        setID(cPos, newID);
        LOGs->setRepresentative(newID, c.Index);
    }
    // Add Remaining new Clusters
    LOGs->addClusters(numNewLOGs - startOfLocalPlog - CommData.PLOGs.size());
    LOGs->clearVolumesAndMerges();
    CommData.PLOGs.clear();

    // First change pointers, then merge (cluster representative information is lost on merge).
    repointerMultipleMerges(merges);
    LOGs->mergeClusterFromList(merges);

    checkConsistency();
}

void LocalBlock::sendData() {
    checkConsistency();

    CommData.PLOGs.clear();
    CommData.PLOGs.reserve(RefPLOGs->size());

    for (ClusterID plog : *RefPLOGs) {
        Cluster c = LOLs->getCluster(plog);
        CommData.PLOGs.emplace_back(c.Index, c.Volume);
        LOLs->removeCluster(plog);
    }
    RefPLOGs->clear();
#ifdef COMMUNICATION
    ind rankTag = Rank << MPICommunication::RANK_SHIFT;

    ind numMessages = 8;
    ind messageId = 0;

#ifdef COLLECTIVES
    numMessages -= 4;
#endif

    MPI_Request* requests = new MPI_Request[numMessages];

    CommData.VectorSizes[0] = static_cast<ind>(CommData.PLOGs.size());
    CommData.VectorSizes[1] = static_cast<ind>(LOGs->Merges.size());

    MPICommunication::handleError(MPI_Isend(CommData.VectorSizes, 2, MPI_IND, 0,
                                            MPICommunication::VECTORSIZES | rankTag, MPI_COMM_WORLD,
                                            &requests[messageId++]));

    // New global clusters in red
    MPICommunication::handleError(
        MPICommunication::IsendVector(CommData.PLOGs, 0, MPICommunication::PLOGS | rankTag,
                                      MPI_COMM_WORLD, &requests[messageId++]));

    // Merges, nothing more to be done here with them
    MPICommunication::handleError(
        MPICommunication::IsendVector(LOGs->Merges, 0, MPICommunication::MERGES | rankTag,
                                      MPI_COMM_WORLD, &requests[messageId++]));

    // Send updated red blocks
    MPICommunication::handleError(MPI_Isend(MemoryLOG, MemoryLOGSize * sizeof(ID), MPI_BYTE, 0,
                                            MPICommunication::REDPOINTERS | rankTag, MPI_COMM_WORLD,
                                            &requests[messageId++]));

    CommData.NumClusters = numClusters();
    CommData.MaxVolume = maxVolume();
    CommData.TotalVolume = totalVolume();

// Number of local Clusters, maxVolume and totalVolume, updated volumes
#ifdef COLLECTIVES
    MPICommunication::handleError(
        MPI_Reduce(&CommData.NumClusters, nullptr, 1, MPI_IND, MPI_SUM, 0, MPI_COMM_WORLD));
    MPICommunication::handleError(
        MPI_Reduce(&CommData.MaxVolume, nullptr, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD));
    MPICommunication::handleError(
        MPI_Reduce(&CommData.TotalVolume, nullptr, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD));
    MPICommunication::handleError(MPI_Reduce(LOGs->volumes().data(), nullptr,
                                             LOGs->volumes().size(), MPI_DOUBLE, MPI_SUM, 0,
                                             MPI_COMM_WORLD));

#else
    MPICommunication::handleError(MPI_Isend(&CommData.NumClusters, 1, MPI_IND, 0,
                                            MPICommunication::NUMCLUSTERS | rankTag, MPI_COMM_WORLD,
                                            &requests[messageId++]));

    MPICommunication::handleError(MPI_Isend(&CommData.MaxVolume, 1, MPI_DOUBLE, 0,
                                            MPICommunication::MAXVOLUME | rankTag, MPI_COMM_WORLD,
                                            &requests[messageId++]));

    MPICommunication::handleError(MPI_Isend(&CommData.TotalVolume, 1, MPI_DOUBLE, 0,
                                            MPICommunication::TOTALVOLUME | rankTag, MPI_COMM_WORLD,
                                            &requests[messageId++]));

    // Volumes for LOGs (Additional volume)
    MPICommunication::handleError(
        MPICommunication::IsendVector(LOGs->volumes(), 0, MPICommunication::VOLUMES | rankTag,
                                      MPI_COMM_WORLD, &requests[messageId++]));
#endif  // COLLECTIVES

#ifndef SINGLENODE
    MPI_Waitall(numMessages, requests, MPI_STATUS_IGNORE);
#else
// delete[] requests;
#endif  // SINGLENODE

#endif  // COMMUNCATION
}  // namespace perc

void LocalBlock::repointerMultipleMerges(const std::vector<ind>& connComps) {
    for (auto it = connComps.begin(); it != connComps.end(); ++it) {
        ind compSize = *it;
        ClusterID ontoCluster(*(++it), false);
        for (ind c = 0; c < compSize - 1; ++c) {
            ClusterID fromCluster(*(++it), false);
            VertexID fromRep = LOGs->getCluster(fromCluster).Index;
            VertexID ontoRep = LOGs->getCluster(ontoCluster).Index;
            // The from rep may not be part of this node (-> Empty Rep)
            if (fromRep.RawID != -1) {
                // Onto Rep not part of this node, make fromRep the Representative
                if (ontoRep.RawID == -1) {
                    setID(fromRep, ontoCluster);
                    LOGs->setRepresentative(ontoCluster, fromRep);
                } else  // Both reps are part of this node
                {
                    assert(fromRep != ontoRep && "Self pointing rep");
                    setID(fromRep, ontoRep);
                }
            }
        }
    }
}

ind LocalBlock::memEstimate() const {
    ind memSize = 0;
    memSize += LOLSubBlock->memEstimate() + LOLs->memEstimate() + LOGs->memEstimate();

    for (auto& log : LOGSubBlocks) memSize += log.memEstimate();
    for (auto& gog : GOGSubBlocks) memSize += gog.memEstimate();

    return memSize;
}

void LocalBlock::checkConsistency() const {
#ifndef NDEBUG
    LOLSubBlock->checkConsistency();
    for (auto& log : LOGSubBlocks) log.checkConsistency();

    for (auto& merge : LOGs->Merges) {
        assert(LOGs->getClusterVolume(merge.From) >= 0 &&
               "Cluster recorded to merge from does not exist.");
        assert(LOGs->getClusterVolume(merge.Onto) > 0 &&
               "Cluster recorded to merge onto does not exist.");
    }
#endif
}

std::vector<std::pair<vec3i, double>> LocalBlock::getVoluminaForAddedVertices(double maxVal) {
    std::vector<std::pair<vec3i, double>> result;

    LOLSubBlock->getVoluminaForAddedVertices(maxVal, result);
    for (auto& log : LOGSubBlocks) {
        log.getVoluminaForAddedVertices(maxVal, result);
    }

    std::sort(result.begin(), result.end(), [this](auto& a, auto& b) {
        return a.first.toIndexOfTotal(TotalSize) > b.first.toIndexOfTotal(TotalSize);
    });

    return result;
}

void LocalBlock::outputFrontBlocks(std::vector<char>& field, ind slice, ind selection) {
    vec3i size = TotalSize;
    size.z = 1;
    assert(field.size() == size.prod() && "Field size not correct.");

    if (selection == 0 || selection == 1) {
        // Mark white block elements.
        if (LOLSubBlock->blockOffset().z <= slice &&
            LOLSubBlock->blockOffset().z + LOLSubBlock->blockSize().z > slice) {
            for (ind x = LOLSubBlock->blockOffset().x;
                 x < LOLSubBlock->blockOffset().x + LOLSubBlock->blockSize().x; ++x)
                for (ind y = LOLSubBlock->blockOffset().y;
                     y < LOLSubBlock->blockOffset().y + LOLSubBlock->blockSize().y; ++y) {

                    field[x + y * size.x] = '.';
                }
        }
    }

    // Mark red block elements.
    if (selection == 0 || selection == 2) {
        for (ind r = 0; r < LOGSubBlocks.size(); ++r) {
            auto& red = LOGSubBlocks[r];
            if (red.blockOffset().z <= slice && red.blockOffset().z + red.blockSize().z > slice) {
                for (ind x = red.blockOffset().x; x < red.blockOffset().x + red.blockSize().x; ++x)
                    for (ind y = red.blockOffset().y; y < red.blockOffset().y + red.blockSize().y;
                         ++y) {
                        assert(field[x + y * size.x] == ' ' && "Overlapping blocks.");
                        field[x + y * size.x] = '0' + r;
                    }
            }
        }
    }

    // Mark gray block elements.
    if (selection == 0 || selection == 3) {
        for (ind g = 0; g < GOGSubBlocks.size(); ++g) {
            auto& gray = GOGSubBlocks[g];
            if (gray.blockOffset().z <= slice &&
                gray.blockOffset().z + gray.blockSize().z > slice) {
                for (ind x = gray.blockOffset().x; x < gray.blockOffset().x + gray.blockSize().x;
                     ++x)
                    for (ind y = gray.blockOffset().y;
                         y < gray.blockOffset().y + gray.blockSize().y; ++y) {
                        //                    assert(field[x + y * size.x] == ' ' &&
                        //                    "Overlapping blocks.");
                        field[x + y * size.x] = 'a' + g;
                    }
            }
        }
    }
}

}  // namespace perc