
# Load required packages
install.packages("TreeDist")
install.packages("ape")
install.packages("phangorn")
install.packages("remotes")

#below codes were modified by cloning the TreeDist github to enable instalation for M1:
#git clone https://github.com/ms609/TreeDist.git
#Edit src/information.h, and change: "constexpr double log_2 = log(2.0);" to "const double log_2 = std::log(2.0);" 

remotes::install_local("/Users/waodedwidaningrat/viridans/tree_comparison/TreeDist")
install.packages("/Users/waodedwidaningrat/viridans/tree_comparison/TreeDist", repos = NULL, type="source")
packageVersion("TreeDist")

system.file(package = "TreeDist")

library(TreeDist)
library(ape)
library(phangorn)
library(TreeTools)

# Increase max supported taxa
options(TreeDist.max.tips = 3000)

setwd("/Users/waodedwidaningrat/viridans/tree_comparison/")
#SUBTREE MITIS
# Read in the two trees from .treefile format
# Replace with your actual file paths
coregene <- read.tree("viridans_coregene_iqtree.treefile")
mlsa <- read.tree("viridans_housekeeping_iqtree.treefile")
recomb_coregene <- read.tree("coregene_prefix.labelled_tree.newick")
recomb_mlsa <- read.tree("mlsa_prefix.labelled_tree.newick")

# Ensure both trees have the same tips [coregene vs mlsa]
common_tips_subcoregene <- intersect(coregene$tip.label, mlsa$tip.label)
coregene <- drop.tip(coregene, setdiff(coregene$tip.label, common_tips_subcoregene))
mlsa <- drop.tip(mlsa, setdiff(mlsa$tip.label, common_tips_subcoregene))

length(coregene$tip.label)
length(mlsa$tip.label)
length(common_tips_subcoregene)

# Ensure both trees have the same tips (coregene vs adjusted-recombination coregene)
common_tips_subcoregenerecomb <- intersect(coregene$tip.label, recomb_coregene$tip.label)
coregene <- drop.tip(coregene, setdiff(coregene$tip.label, common_tips_subcoregenerecomb))
recomb_coregene <- drop.tip(recomb_coregene, setdiff(recomb_coregene$tip.label, common_tips_subcoregenerecomb))

length(coregene$tip.label)
length(recomb_coregene$tip.label)
length(common_tips_subcoregenerecomb)

# Ensure both trees have the same tips )mlsa vs adjusted-recombination mlsa)
common_tips_submlsa <- intersect(mlsa$tip.label, recomb_mlsa$tip.label)
mlsa <- drop.tip(mlsa, setdiff(mlsa$tip.label, common_tips_subcoregene))
recomb_mlsa <- drop.tip(recomb_mlsa, setdiff(recomb_mlsa$tip.label, common_tips_submlsa))

length(mlsa$tip.label)
length(recomb_mlsa$tip.label)
length(common_tips_subcoregene)

# Compute RF (used this one for the manuscript) :
# coregene vs mlsa
RF.dist(coregene, mlsa, normalize = TRUE)
# coregene vs coregene recombination-adjusted
RF.dist(coregene, recomb_coregene, normalize = TRUE)
# mlsa vs mlsa recombination-adjusted
RF.dist(mlsa, recomb_mlsa, normalize = TRUE)


# Compute Generalized RF
grf_subcoregene <- JaccardRobinsonFoulds(coregene, mlsa)

# alternative to InfoRobinsonFoulds
InfoDist(coregene, mlsa)



# Compute Generalized RF
grf_submlsa <- JaccardRobinsonFoulds(tree3, tree4)

# Compute Generalized RF
grf_subtree3 <- JaccardRobinsonFoulds(tree5, tree6)

cat("Generalized RF distance for core gene vs mlsa:", round(grf_subcoregene, 5), "\n")
cat("Generalized RF distance for mitis subcluster I:", round(grf_submlsa, 5), "\n")
cat("Generalized RF distance for mitis subtree I with pneumo collaped:", round(grf_subtree3, 5), "\n")

ClusteringEntropy(coregene)
ClusteringEntropy(mlsa)
SplitwiseInfo(coregene)

SharedPhylogeneticInfo(coregene, mlsa)
MutualClusteringInfo(coregene, mlsa)
InfoRobinsonFoulds(coregene, mlsa)

SharedPhylogeneticInfo(tree3, tree4)
MutualClusteringInfo(tree3, tree4)
InfoRobinsonFoulds(tree3, tree4)

SharedPhylogeneticInfo(tree5, tree6)
MutualClusteringInfo(tree5, tree6)
InfoRobinsonFoulds(tree5, tree6)

#Smith (2020) recommend MCI (mutual clustering info) as a metric to quantify tree comparison compare to RF
#SharePhylogeneticInfo can be used as the maximum phylogenetic info two trees can share
#so if the MCI was 50, while the SPI is 1,000, then 50/1,000 is 5% meaning that the two tree only share 5% similarity

mci <- MutualClusteringInfo(coregene, mlsa)  # Mutual clustering info between coregene and mlsa
similarity <- mci / SharedPhylogeneticInfo(coregene, coregene)  # Normalized similarity
mci
similarity

#read more here to confirm undestanding:
#https://cran.r-project.org/web/packages/TreeDist/vignettes/treespace.html
#https://cran.r-project.org/web/packages/TreeDist/vignettes/Generalized-RF.html

# Use VisualizeMatching with the custom plot function
TwoTreePlot()
VisualizeMatching(RobinsonFouldsMatching, coregene, mlsa)
VisualizeMatching(MutualClusteringInfo, coregene, mlsa)
VisualizeMatching(SharedPhylogeneticInfo, coregene, mlsa)
VisualizeMatching(InfoRobinsonFoulds, coregene, mlsa)

noLabelsPlot <- function(tree, ...) {
  tree$node.label <- NULL  # Remove support values (node labels)
  plot(tree, show.tip.label = FALSE, ...)  # Plot without tip labels
}

# Define a custom plot function to hide leaf labels and bootstrap values
noLabelPlot <- function(tree, ...) {
  plot(tree, show.tip.label = FALSE, edge.label = NULL, ...)
}

VisualizeMatching(RobinsonFouldsMatching, coregene, mlsa, Plot = noLabelPlot)
VisualizeMatching(MutualClusteringInfo, coregene, mlsa, Plot = noLabelPlot)
VisualizeMatching(SharedPhylogeneticInfo, coregene, mlsa, Plot = noLabelPlot)
VisualizeMatching(InfoRobinsonFoulds, coregene, mlsa, Plot = noLabelPlot)

VisualizeMatching(RobinsonFouldsMatching, tree3, tree4, Plot = noLabelPlot)
VisualizeMatching(MutualClusteringInfo, tree3, tree4, Plot = noLabelPlot)
VisualizeMatching(SharedPhylogeneticInfo, tree3, tree4, Plot = noLabelPlot)
VisualizeMatching(InfoRobinsonFoulds, tree3, tree4, Plot = noLabelPlot)

VisualizeMatching(RobinsonFouldsMatching, tree5, tree6, Plot = noLabelPlot)
VisualizeMatching(MutualClusteringInfo, tree5, tree6, Plot = noLabelPlot)
VisualizeMatching(SharedPhylogeneticInfo, tree5, tree6, Plot = noLabelPlot)
VisualizeMatching(InfoRobinsonFoulds, tree5, tree6, Plot = noLabelPlot)

ClusteringInfo(coregene)
ClusteringInfo(mlsa)

CIcoregene <- ClusteringInfo(coregene)
CImlsa <- ClusteringInfo(mlsa)
mutual_infoCI <- MutualClusteringInfo(coregene, mlsa)

normalized_similarityCI <- 2 * mutual_infoCI / (CIcoregene + CImlsa)

CIcoregene
CImlsa
mutual_infoCI
normalized_similarityCI

PIcoregene <- PhylogeneticInfo(coregene)
PImlsa <- PhylogeneticInfo(mlsa)
mutual_infoPI <- SharedPhylogeneticInfo(coregene, mlsa)

normalized_similarityPI <- 2 * mutual_infoPI / (PIcoregene + PImlsa)

PIcoregene
PImlsa
mutual_infoPI
normalized_similarityPI

TreeTools::NSplits(coregene)
NyeSimilarity(coregene, mlsa)
NyeSimilarity(coregene, mlsa, normalize = TreeTools::NSplits(coregene))

#The convenience function TreeDistance() returns the variation of clustering information between two trees,
#normalized against the total information content of all splits.

distance <- TreeDistance(coregene, mlsa)

distance


legend("topright", legend = c("Shared", "Unique to Tree 1", "Unique to Tree 2"),
       col = c("black", "red", "blue"), lty = 1, bty = "n", cex = 0.8)

matches <- InfoRobinsonFoulds(coregene, mlsa, reportMatching = TRUE)
summary(matches)
        

VisualizeMatching(SharedPhylogeneticInfo, coregene, mlsa, 
                  Plot = TreeDistPlot, matchZeros = FALSE)

VisualizeMatching(RobinsonFouldsMatching, coregene, mlsa, Plot = TreeDistPlot)

VisualizeMatching(MutualClusteringInfo, coregene, mlsa, 
                  Plot = TreeDistPlot, matchZeros = FALSE)


TreeDist::MapTrees()




#FULLTREE
# Read in the two trees from .treefile format
# Replace with your actual file paths
treecore <- read.tree("coregene_viridans_fasttree.newick")
treemlsa <- read.tree("housekeeping_viridans_fasttree.newick")

# Ensure both trees have the same tips
common_tips_fulltree <- intersect(treecore$tip.label, treemlsa$tip.label)
treecore <- drop.tip(treecore, setdiff(treecore$tip.label, common_tips_fulltree))
treemlsa <- drop.tip(treemlsa, setdiff(treemlsa$tip.label, common_tips_fulltree))

length(treecore$tip.label)
length(treemlsa$tip.label)
length(common_tips_fulltree)

# Compute Generalized RF
grf_fulltree <- JaccardRobinsonFoulds(treecore, treemlsa)

cat("Generalized RF distance for full tree:", round(grf_fulltree, 5), "\n")

