#include "globalprocessor.h"
#include "unionfindsubblock.h"

namespace perc {

void GlobalProcessor::setParent(UnionFindSubBlock<GlobalProcessor>* parent) {
    assert(!Parent && "Parent already set!");
    assert(&(parent->NeighborProcessor) == this && "This is not our parent.");

    Parent = parent;
}

ID GlobalProcessor::doWatershed(VertexID pos, double volume, std::vector<Neighbor>& neighClusters) {
    switch (neighClusters.size()) {

        // Completely new GOG.
        case 0: {
            return GOGs.addCluster(pos, volume, Parent);
        }
        // Extend GOG
        case 1: {
            GOGs.extendCluster(neighClusters[0].Cluster, volume);
            // If we are extending a cluster that is not in our current block, we might be a new
            // representative
            if (!Parent->contains(neighClusters[0].Representative)) {
                // Will either be us or the representative in the parent block
                VertexID rep = GOGs.setRepresentative(neighClusters[0].Cluster, pos, false, Parent);
                // We are the new rep for the block, point to cluster directly
                if (rep == pos) {
                    return neighClusters[0].Cluster;
                }
                // Point to representative
                else {
                    return rep;
                }
            }
            return neighClusters[0].Representative.toIndexOfTotal(Parent->totalSize());
        }
        // Record merges
        default: {
            // All clusters are global, so we just merge on the first one
            Neighbor& gog = *neighClusters.begin();

            GOGs.extendCluster(gog.Cluster, volume);
            for (Neighbor& neigh : neighClusters) {
                // Don't merge gog onto itself.
                if (neigh.Cluster == gog.Cluster) continue;
                /*const std::vector<GOG>& reps = GOGs.getRepresentatives(neigh.Cluster);
                // Directly point representativevs to the cluster to avoid duplicate merges
                for (auto rep : reps) {
                    UnionFindSubBlock<GlobalProcessor>* blockParent =
                        reinterpret_cast<UnionFindSubBlock<GlobalProcessor>*>(rep.ParentBlock);
                    blockParent->PointerBlock.setPointer(
                        vec3i::fromIndexOfTotal(rep.ID.baseID(), Parent->totalSize()), gog.Cluster);
                }*/
                GOGs.mergeClusters(neigh.Cluster, gog.Cluster);
            }

            // Check if we need to be a representative of cluster we are merging onto, see extend
            // case
            // TODO: The cluster we are merging from might already have reps in this block?
            if (!(Parent->contains(gog.Representative))) {
                // Will either be us or the representative in the parent block
                VertexID rep = GOGs.setRepresentative(gog.Cluster, pos, false, Parent);
                // We are the new rep for the block, point to cluster directly
                if (rep == pos) {
                    return neighClusters[0].Cluster;
                }
                // Point to representative
                else {
                    return rep;
                }
            }

            return gog.Representative.toIndexOfTotal(Parent->totalSize());
        }
    }
}

}  // namespace perc