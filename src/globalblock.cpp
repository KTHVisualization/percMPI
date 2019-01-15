#include "globalblock.h"
#include "performancetimer.h"
#include "mpicommuncation.h"
#include "interfaceblockbuilder.h"

namespace perc {

GlobalBlock::GlobalBlock(const vec3i& blockSize, const vec3i& blockOffset, const vec3i& totalSize,
                         const vec3i& numNodes)
    : NumNodes(numNodes), UnionFindBlock(totalSize), GOGs(GLOBAL_LIST) {

    // TODO: Actually load green blocks in correct locations, and set up number of red blocks
    int numBlocks = 0;  // Currently set to 0 for communications test
    float fraction = 1.0 / numBlocks;
    vec3i subBlocksize = vec3i(blockSize.x, blockSize.y, ind(blockSize.z * fraction));
    vec3i subOffset = vec3i(0, 0, ind(blockSize.z * fraction));

    GOGSubBlocks.reserve(numBlocks);

    for (int i = 0; i < numBlocks; i++) {
        int offsetScale = numBlocks - i;
        UnionFindSubBlock<GreenProcessor> newGOGBlock = UnionFindSubBlock<GreenProcessor>(
            subBlocksize, {subOffset.x * i, subOffset.y * i, subOffset.z * i}, totalSize, *this,
            GreenProcessor(GOGs));
        GOGSubBlocks.push_back(newGOGBlock);
    }

    for (auto gogBlock : GOGSubBlocks) {
        gogBlock.loadData();
    }

    // Set up data for sending and receiving as well as red blocks
    ReceivedMerges.resize(NumNodes.prod());

    for (ind p = 0; p < numNodes.prod(); ++p) {
        InfoPerProcess dataPerProcess;

        dataPerProcess.Merges = &ReceivedMerges[p];

        // TODO: Add correct greenIndices
        std::vector<ind> greenIndices;
        dataPerProcess.GreenIndices = greenIndices;

        auto constructor = [this]() -> GrayProcessor { return GrayProcessor(); };
        vec3i whiteBlockSize, whiteBlockOffset;

        interfaceblockbuilder::buildRedBlocks<GrayProcessor>(
            blockSize, blockOffset, totalSize, LOGSubBlocks, dataPerProcess.MemoryLOG,
            dataPerProcess.MemoryLOGSize, whiteBlockSize, whiteBlockOffset, *this, constructor);

        PerProcessData.push_back(dataPerProcess);
    }
}

void GlobalBlock::doWatershed(const double minVal) {
    ind numClusters = GOGs.numClusters();
    for (auto gogBlock : GOGSubBlocks) {
        gogBlock.doWatershed(minVal);
    }
    NumNewClusters += numClusters - GOGs.numClusters();

    // Merges have only been recorded -> Do merges now
    // Include
    ReceivedMerges.push_back(GOGs.Merges);
    // PerformanceTimer timer;
    // timer.Reset();
    // std::cout << "Merge ";
    Merges = ClusterMerge::mergeClusterAsList(ClusterMerge::mergeClustersFromLists(ReceivedMerges));
    // std::cout << "took " << timer.ElapsedTimeAndReset() << " seconds." << std::endl;

    // Merges are processed, clear for next step (remove from received again)
    GOGs.Merges.clear();
    ReceivedMerges.pop_back();

    // Do the representatives merge and repointering first
    repointerMultipleMerges(Merges);
    // Merge Clusters and Representatives in the list (this process will be recomputed by all other
    // blocks that receive the merges)
    GOGs.mergeClusterFromList(Merges);

    checkConsistency();
}

ClusterID* GlobalBlock::findClusterID(const vec3i& idx, vec3i& lastClusterID) {
    // TODO: Do this more cleverly, and for all types of subblocks:
    // One might want to do this more cleverly, especialy in the sheet tree.
    for (auto gogBlock : GOGSubBlocks)
        if (gogBlock.contains(idx)) return gogBlock.findClusterID(idx, lastClusterID);
    for (auto logBlock : LOGSubBlocks) {
        if (logBlock.contains(idx)) return logBlock.findClusterID(idx, lastClusterID);
    }
    assert(false && "Can not find block containing this idx.");
    return nullptr;
}

ID* GlobalBlock::setID(const vec3i& idx, const ID& id) {
    // One might want to do this more cleverly, especialy in the sheet tree.
    ID* ptr = nullptr;
    for (auto gogBlock : GOGSubBlocks)
        if (gogBlock.contains(idx)) {
            ptr = gogBlock.PointerBlock.getPointer(idx);
            break;
        }

    if (!ptr) {
        for (auto logBlock : LOGSubBlocks)
            if (logBlock.contains(idx)) {
                ptr = logBlock.PointerBlock.getPointer(idx);
                break;
            }
    }

    assert(ptr && "Can not find block containing this idx.");
    *ptr = id;

    return ptr;
}

double GlobalBlock::getClusterVolume(ClusterID cluster) { return GOGs.getClusterVolume(cluster); }

void GlobalBlock::repointerMultipleMerges(const std::vector<ind>& connComps) {
    for (auto it = connComps.begin(); it != connComps.end(); ++it) {
        ind compSize = *it;
        ind onto = *(++it);
        for (ind c = 0; c < compSize - 1; ++c) {
            int from = *(++it);
            const std::vector<GOG>& oldRepsFrom = GOGs.getRepresentatives(from);
            const std::vector<VertexID> newRepsFrom = GOGs.mergeRepresentatives(from, onto);
            assert(oldRepsFrom.size() == newRepsFrom.size() &&
                   "Size mismatch between old and new representatives");
            auto itOld = oldRepsFrom.begin();
            for (auto itNew = newRepsFrom.begin();
                 itNew != newRepsFrom.end() && itOld != oldRepsFrom.end(); itOld++, itNew++) {
                UnionFindSubBlock<GreenProcessor>* parent =
                    reinterpret_cast<UnionFindSubBlock<GreenProcessor>*>(itOld->ParentBlock);
                // Rep is a new rep for the merged cluster -> point to it
                if (*itNew == itOld->ID) {
                    parent->PointerBlock.setPointer(
                        vec3i::fromIndexOfTotal(itOld->ID.baseID(), parent->totalSize()), onto);
                }
                // There is already a rep in the same block -> point to the rep from onto
                else {
                    parent->PointerBlock.setPointer(
                        vec3i::fromIndexOfTotal(itOld->ID.baseID(), parent->totalSize()), *itNew);
                }
            }
        }
    }
}

void GlobalBlock::receiveData() {
    // Init / Reset datastructures
    NumNewClusters = 0;
    NumClustersLocal = 0;
    MaxVolumeLocal = 0;
    TotalVolumeLocal = 0;
    ind plogsAddedSoFar = 0;
    Merges.clear();

    // TODO: Use p + 1, since 0 will be the one green block
    for (int p = 0; p < NumNodes.prod(); p++) {
        ind numMessages = 7;
        MPI_Status status;
        int err;

        // Receive number of local Clusters and maxVolume and totalVolume (Might want to be gathers)
        int numClusters;
        err == MPI_Recv(&numClusters, 1, MPI_INT, p, MPICommunication::NUMCLUSTERS, MPI_COMM_WORLD,
                        &status);
        NumClustersLocal += numClusters;
        double maxVolume;
        err = MPI_Recv(&maxVolume, 1, MPI_DOUBLE, p, MPICommunication::MAXVOLUME, MPI_COMM_WORLD,
                       &status);
        MaxVolumeLocal = std::max(maxVolume, MaxVolumeLocal);
        double volume;
        err = MPI_Recv(&volume, 1, MPI_DOUBLE, p, MPICommunication::TOTALVOLUME, MPI_COMM_WORLD,
                       &status);
        TotalVolumeLocal += volume;

        // Receive volumes for LOGs (Additional volume) and update volumes (Should be some collected
        // for this too)
        std::vector<double> commVolumes(GOGs.volumes().size());
        err = MPI_Recv(commVolumes.data(), commVolumes.size(), MPI_DOUBLE, p,
                       MPICommunication::VOLUMES, MPI_COMM_WORLD, &status);

        ind counter = 0;
        for (double vol : commVolumes) {
            if (vol != 0) {
                GOGs.extendCluster(ClusterID(counter, false), vol);
            }
            counter++;
        }

        // Receive PLOGs and add a new cluster for each and initialize with their respective
        // volume
        std::vector<ClusterData> commPLOGs;
        err = MPICommunication::RecvVectorUknownSize(commPLOGs, p, MPICommunication::PLOGS,
                                                     MPI_COMM_WORLD, &status);
        // Record where PLOGS for this process start (as in: How many other PLOGs for other
        // processes where added before)
        PerProcessData[p].StartOfLocalPlog = plogsAddedSoFar;
        NumNewClusters += commPLOGs.size();
        plogsAddedSoFar += commPLOGs.size();
        for (auto cluster : commPLOGs) {
            UnionFindSubBlock<GrayProcessor>* parentBlock;
            vec3i cPos = vec3i::fromIndexOfTotal(cluster.Index.RawID, TotalSize);
            for (auto logBlock : LOGSubBlocks)
                if (logBlock.contains(cPos)) {
                    parentBlock = &logBlock;
                    break;
                }
            assert(parentBlock && "No Block for PLOG found");
            ClusterID newCluster = GOGs.addCluster(cluster.Index, cluster.Volume, parentBlock);
            parentBlock->PointerBlock.setPointer(cPos, newCluster);
        }

        // Receive merges, nothing more to be done here with them
        std::vector<ClusterMerge>& merges = ReceivedMerges[p];
        err = MPICommunication::RecvVectorUknownSize(merges, p, MPICommunication::MERGES,
                                                     MPI_COMM_WORLD, &status);

        // Receive updated red blocks
        err = MPI_Recv(PerProcessData[p].MemoryLOG, PerProcessData[p].MemoryLOGSize * sizeof(ID),
                       MPI_BYTE, p, MPICommunication::REDPOINTERS, MPI_COMM_WORLD, &status);
    }
}

void GlobalBlock::sendData() {
    // Broadcast number of new Clusters and updated Merges
    /* MPI_Bcast(&NumNewClusters, 1, MPI_INT, 0, MPI_COMM_WORLD);
    int mergesSize = Merges.size();
    // Probably Necessary to broadcast size of merges first
    MPI_Bcast(&mergesSize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&Merges[0], Merges.size() * sizeof(ind), MPI_BYTE, 0, MPI_COMM_WORLD);
    */

    // TODO: Use p + 1, since 0 will be the one green block
    for (ind p = 0; p < NumNodes.prod(); p++) {
        // NumNewClusters, Merges Vector, PLOG range,
        ind numMessages = 3 + PerProcessData[p].GreenIndices.size();
        MPI_Request* requests = new MPI_Request[numMessages];

        // TODO: Change this back to a Broadcast, see above
        MPI_Isend(&NumNewClusters, 1, MPI_INT, 0, MPICommunication::NUMNEWCLUSTERS, MPI_COMM_WORLD,
                  &requests[0]);
        MPICommunication::IsendVector(Merges, 0, MPICommunication::MERGES, MPI_COMM_WORLD,
                                      &requests[1]);

        // Plog Range (This might want to be a scatter operation)
        MPI_Isend(&PerProcessData[p].StartOfLocalPlog, 1, MPI_INT, 0, MPICommunication::STARTOFPLOG,
                  MPI_COMM_WORLD, &requests[2]);

        // Updated green blocks
        std::vector<ind>& greenIndices = PerProcessData[p].GreenIndices;
        int counter = 0;
        for (auto id : greenIndices) {
            auto gogBlock = GOGSubBlocks[id];
            MPI_Isend(gogBlock.PointerBlock.PointerBlock, gogBlock.blockSize().prod() * sizeof(ID),
                      MPI_BYTE, 0, MPICommunication::GREENPOINTERS & counter, MPI_COMM_WORLD,
                      &requests[3 + counter]);
        }

        // This should free our requests as well??? -> Would mean blocking after each process, might
        // not be what we want
        // TODO:
        // MPI_Waitall(numMessages, requests, MPI_STATUS_IGNORE);
    }
}

void GlobalBlock::checkConsistency() const {
#ifndef NDEBUG
    for (auto gog : GOGSubBlocks) gog.checkConsistency();

    for (auto& merge : GOGs.Merges) {
        assert(GOGs.getClusterVolume(merge.From) >= 0 &&
               "Cluster recorded to merge from does not exist.");
        assert(GOGs.getClusterVolume(merge.Onto) > 0 &&
               "Cluster recorded to merge onto does not exist.");
    }
#endif
}

std::vector<std::pair<vec3i, double>> GlobalBlock::getVoluminaForAddedVertices(double maxVal) {
    std::vector<std::pair<vec3i, double>> result;

    for (auto gog : GOGSubBlocks) {
        gog.getVoluminaForAddedVertices(maxVal, result);
    }

    std::sort(result.begin(), result.end(), [this](auto& a, auto& b) {
        return a.first.toIndexOfTotal(TotalSize) > b.first.toIndexOfTotal(TotalSize);
    });

    return result;
}

}  // namespace perc