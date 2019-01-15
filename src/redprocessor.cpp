#include "redprocessor.h"
#include "unionfindsubblock.h"

namespace perc {

void RedProcessor::setParent(UnionFindSubBlock<RedProcessor>* parent) {
    assert(!Parent && "Parent already set!");
    assert(&(parent->NeighborProcessor) == this && "This is not our parent.");

    Parent = parent;
}

ID RedProcessor::doWatershed(VertexID pos, double volume, std::vector<Neighbor>& neighClusters) {
    switch (neighClusters.size()) {
        // New LOL that is directly marked as PLOG.
        case 0: {
            ClusterID plog = LOLs.addCluster(pos, volume);
            PLOGs.insert(plog);
            return plog;
        }

        // Extend LOL/LOG, set pointer to representative.
        case 1: {
            Neighbor& neigh = neighClusters[0];
            if (neigh.Cluster.isGlobal()) {
                LOGs.extendCluster(neigh.Cluster, volume);
                // No LOG for GOG yet?
                if (!LOGs.getRepresentative(neigh.Cluster).isValid()) {
                    LOGs.setRepresentative(neigh.Cluster, pos);
                    return pos;
                }
                return neigh.Representative.toIndexOfTotal(Parent->totalSize());
            } else {
                LOLs.extendCluster(neigh.Cluster, volume);

                // Mark as PLOG iff it isn't already. Make this the new representative.
                if (PLOGs.insert(neigh.Cluster).second) {
                    // Make this position the cluster representative.
                    Parent->Parent.setID(neigh.Representative, pos);
                    LOLs.setRepresentative(neigh.Cluster, pos);
                    return neigh.Cluster;
                } else
                    return neigh.Representative.toIndexOfTotal(Parent->totalSize());
            }
        }

        // Merge onto LOG, if any.
        default:
            auto log = std::find_if(neighClusters.begin(), neighClusters.end(),
                                    [](const Neighbor& neigh) { return neigh.Cluster.isGlobal(); });

            // Merge all onto a random LOG.
            VertexID mergeDest;
            if (log != neighClusters.end()) {
                // In case of green: use this as representative.
                if (!LOGs.getRepresentative(log->Cluster).isValid()) {
                    LOGs.setRepresentative(log->Cluster, pos);
                    log->Representative = vec3i::fromIndexOfTotal(pos.RawID, Parent->totalSize());
                }
                mergeDest = log->Representative.toIndexOfTotal(Parent->totalSize());
                LOGs.extendCluster(log->Cluster, volume);
                for (Neighbor& neigh : neighClusters) {
                    // Don't merge log onto itself.
                    if (neigh.Cluster == log->Cluster) continue;
                    if (neigh.Cluster.isGlobal()) {  // Merge LOG->LOG.
                        LOGs.mergeClusters(neigh.Cluster, log->Cluster);
                    } else {  // Merge LOL->LOG.
                        LOGs.extendCluster(log->Cluster, LOLs.getClusterVolume(neigh.Cluster));
                        PLOGs.erase(neigh.Cluster);
                        LOLs.removeCluster(neigh.Cluster);
                        if (Parent->contains(neigh.Representative))
                            Parent->PointerBlock.setPointer(neigh.Representative, mergeDest);
                        else
                            Parent->Parent.setID(neigh.Representative, mergeDest);
                    }
                }
            }

            // No global cluster.
            else {
                Neighbor* dest = nullptr;

                // PLOG representative must lay in red block, so merge onto PLOG if any.
                for (auto& plog : neighClusters) {
                    if (PLOGs.count(plog.Cluster)) {
                        dest = &plog;
                        break;
                    }
                }

                assert(dest && "No LOG or PLOG neigbor found - impossible at a red merge.");
                mergeDest = dest->Representative.toIndexOfTotal(Parent->totalSize());

                // Extend by the volume of the voxel that has caused the merge
                LOLs.extendCluster(dest->Cluster, volume);
                for (auto neigh = neighClusters.begin(); neigh != neighClusters.end(); ++neigh) {
                    if (dest->Cluster == neigh->Cluster) continue;

                    PLOGs.erase(neigh->Cluster);
                    LOLs.mergeClusters(neigh->Cluster, dest->Cluster);
                    Parent->Parent.setID(neigh->Representative, mergeDest);
                }
            }
            return mergeDest;
    }
}

void RedProcessor::checkConsistency() const {
#ifndef NDEBUG
    double totVolume = 0;
    for (ind i = 0; i < Parent->blockSize().prod(); ++i) {
        ID curr = Parent->PointerBlock.PointerBlock[i];
        vec3i currPos = Parent->blockOffset() + vec3i::fromIndexOfTotal(i, Parent->blockSize());
        VertexID currPosID = currPos.toIndexOfTotal(Parent->totalSize());

        // Was processed yet?
        if (curr.baseID() >= 0) {
            totVolume += 1.0;
            // Iff pointing to cluster directly, assert that LOG/PLOG and correct back-pointer.
            if (curr.isCluster()) {
                ClusterID* currClust = curr.asCluster();
                if (currClust->isGlobal()) {
                    assert(LOGs.getRepresentative(*currClust).RawID >= 0 &&
                           "No valid representative set in pointed-to cluster..");
                    assert(LOGs.getRepresentative(*currClust) == currPosID &&
                           "Directly pointing non-representative.");
                } else {
                    assert(std::find(PLOGs.begin(), PLOGs.end(), *currClust) != PLOGs.end() &&
                           "Red can not point to LOLs that are not PLOGs.");
                    assert(LOLs.getCluster(*currClust).Index == currPosID &&
                           "Directly pointing non-representative.");
                }
            }
        }
    }
#endif
}

}  // namespace perc