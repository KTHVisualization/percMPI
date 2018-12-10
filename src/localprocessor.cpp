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
                for (Neighbor& neigh : neighClusters) {
                    if (neigh.Cluster == log->Cluster)
                        continue;  // TODO: Should that case not be already excluded?
                    if (neigh.Cluster.isGlobal()) {
                        LOGs.mergeClusters(neigh.Cluster, log->Cluster);
                        Parent->PointerBlock.setPointer(neigh.Representative, mergeDest);
                    } else {
                        LOGs.extendCluster(log->Cluster, LOLs.getCluster(neigh.Cluster).Volume);
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