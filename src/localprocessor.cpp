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
                        LOGs.extendCluster(log->Cluster, LOLs.getCluster(neigh.Cluster).Volume);
                        LOLs.removeCluster(neigh.Cluster);

                        // Remove possible appearance in PLOG list.
                        auto found = std::find(PLOGs.begin(), PLOGs.end(), neigh.Cluster);
                        if (found != PLOGs.end()) {
                            *found = PLOGs[PLOGs.size() - 1];
                            PLOGs.pop_back();
                        }

                        Parent->PointerBlock.setPointer(neigh.Representative, mergeDest);
                    }
                }
            }

            // No global cluster.
            else {
                auto lol = neighClusters[0];
                mergeDest = lol.Representative.toIndexOfTotal(Parent->totalSize());
                bool lolIsPlog = std::find(PLOGs.begin(), PLOGs.end(), lol.Cluster) != PLOGs.end();
                for (auto neigh = ++neighClusters.begin(); neigh != neighClusters.end(); ++neigh) {
                    LOLs.mergeClusters(neigh->Cluster, lol.Cluster);
                    Parent->PointerBlock.setPointer(neigh->Representative, mergeDest);
                    // Extend by the volume of the voxel that has caused the merge
                    LOLs.extendCluster(lol.Cluster, volume);

                    // Remove possible appearance in PLOG list.
                    auto found = std::find(PLOGs.begin(), PLOGs.end(), neigh->Cluster);
                    if (found != PLOGs.end()) {
                        if (!lolIsPlog) {
                            *found = lol.Cluster;
                            lolIsPlog = true;
                        } else {
                            *found = PLOGs[PLOGs.size() - 1];
                            PLOGs.pop_back();
                        }
                    }
                }
            }
            vec3i posIdx = vec3i::fromIndexOfTotal(pos.baseID(), Parent->totalSize());
            return mergeDest;
            break;
    }
}

}  // namespace perc