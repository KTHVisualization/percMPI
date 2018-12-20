#include "localprocessor.h"
#include "unionfindsubblock.h"

namespace perc {

void LocalLocalProcessor::setParent(UnionFindSubBlock<LocalLocalProcessor>* parent) {
    assert(!Parent && "Parent already set!");
    assert(&(parent->NeighborProcessor) == this && "This is not our parent.");

    Parent = parent;
}

ID LocalLocalProcessor::doWatershed(VertexID pos, double volume,
                                    std::vector<Neighbor>& neighClusters) {
    switch (neighClusters.size()) {
        // New LOL.
        case 0:
            return LOLs.addCluster(pos, volume);

        // Extend LOL or LOG, set pointer to representative.
        case 1: {
            Neighbor& neigh = neighClusters[0];
            if (neigh.Cluster.isGlobal())
                LOGs.extendCluster(neigh.Cluster, volume);
            else
                LOLs.extendCluster(neigh.Cluster, volume);
            return neigh.Representative.toIndexOfTotal(Parent->totalSize());
        }

        // Merge onto LOG, if any.
        default:
            auto log = std::find_if(neighClusters.begin(), neighClusters.end(),
                                    [](const Neighbor& neigh) { return neigh.Cluster.isGlobal(); });

            // Merge all onto an random LOG.
            VertexID mergeDest;
            if (log != neighClusters.end()) {
                mergeDest = log->Representative.toIndexOfTotal(Parent->totalSize());
                LOGs.extendCluster(log->Cluster, volume);
                for (Neighbor& neigh : neighClusters) {
                    // Don't merge log onto itself.
                    if (neigh.Cluster == log->Cluster) continue;

                    // LOG->LOG merge.
                    if (neigh.Cluster.isGlobal()) {
                        LOGs.mergeClusters(neigh.Cluster, log->Cluster);
                    } else {  // LOL->LOG merge.
                        LOGs.extendCluster(log->Cluster, LOLs.getClusterVolume(neigh.Cluster));
                        LOLs.removeCluster(neigh.Cluster);

                        // Remove possible appearance in PLOG list.
                        PLOGs.erase(neigh.Cluster);

                        Parent->PointerBlock.setPointer(neigh.Representative, mergeDest);
                    }
                }
            }

            // No global cluster.
            else {
                auto lol = neighClusters[0];
                mergeDest = lol.Representative.toIndexOfTotal(Parent->totalSize());
                auto lolIsPlog = PLOGs.find(lol.Cluster);
                for (auto neigh = ++neighClusters.begin(); neigh != neighClusters.end(); ++neigh) {
                    LOLs.mergeClusters(neigh->Cluster, lol.Cluster);
                    Parent->PointerBlock.setPointer(neigh->Representative, mergeDest);
                    // Extend by the volume of the voxel that has caused the merge
                    LOLs.extendCluster(lol.Cluster, volume);

                    // Remove possible appearance in PLOG list.
                    if (PLOGs.erase(neigh->Cluster) > 0)  // Erased an element.
                        PLOGs.insert(lol.Cluster);
                }
            }
            vec3i posIdx = vec3i::fromIndexOfTotal(pos.baseID(), Parent->totalSize());
            return mergeDest;
    }
}

// Local Global Processor //

void LocalGlobalProcessor::setParent(UnionFindSubBlock<LocalGlobalProcessor>* parent) {
    assert(!Parent && "Parent already set!");
    assert(&(parent->NeighborProcessor) == this && "This is not our parent.");

    Parent = parent;
}

ID LocalGlobalProcessor::doWatershed(VertexID pos, double volume,
                                     std::vector<Neighbor>& neighClusters) {
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
                if (!LOGs.getRepresentative(neigh.Cluster).isValid())
                    LOGs.setRepresentative(neigh.Cluster, pos);
            } else {
                LOLs.extendCluster(neigh.Cluster, volume);

                // Mark as PLOG iff it isn't already. Make this the new representative.
                if (PLOGs.insert(neigh.Cluster).second) {
                    // Make this position the cluster representative.
                    ID* lolRep = Parent->PointerBlock.getPointer(neigh.Representative);
                    *lolRep = pos;
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

            // Merge all onto an random LOG.
            VertexID mergeDest;
            if (log != neighClusters.end()) {
                mergeDest = log->Representative.toIndexOfTotal(Parent->totalSize());
                LOGs.extendCluster(log->Cluster, volume);
                for (Neighbor& neigh : neighClusters) {
                    // Don't merge log onto itself.
                    if (neigh.Cluster == log->Cluster) continue;
                    if (neigh.Cluster.isGlobal()) {  // Merge LOG->LOG.
                        LOGs.mergeClusters(neigh.Cluster, log->Cluster);
                    } else {  // Merge LOL->LOG.
                        LOGs.extendCluster(log->Cluster, LOLs.getClusterVolume(neigh.Cluster));
                        LOLs.removeCluster(neigh.Cluster);
                        Parent->PointerBlock.setPointer(neigh.Representative, mergeDest);
                    }
                }
            }

            // No global cluster.
            else {
                auto lol = neighClusters[0];
                mergeDest = lol.Representative.toIndexOfTotal(Parent->totalSize());
                for (auto neigh = neighClusters.begin() + 1; neigh != neighClusters.end();
                     ++neigh) {

                    LOLs.mergeClusters(neigh->Cluster, lol.Cluster);
                    Parent->PointerBlock.setPointer(neigh->Representative, mergeDest);
                    // Extend by the volume of the voxel that has caused the merge
                    LOLs.extendCluster(lol.Cluster, volume);
                }
            }
            vec3i posIdx = vec3i::fromIndexOfTotal(pos.baseID(), Parent->totalSize());
            return mergeDest;
            break;
    }
}

}  // namespace perc