library(ggplot2)
library(phyloseq)
library(reshape2)
library(ape)
library(gridExtra)
library(plotly)

source("graphical_methods.R")
source("tree_methods.R")
source("plot_merged_trees.R")
source("specificity_methods.R")
source("ternary_plot.R")
source("richness.R")
source("edgePCA.R")
source("copy_number_correction.R")
source("import_frogs.R")
source("prevalence.R")
source("compute_niche.R")

## Import Chaillou data
food <- import_biom(BIOMfilename = "test/Chaillou/chaillou.biom", 
                    treefilename = "test/Chaillou/tree.nwk",
                    parseFunction = parse_taxonomy_greengenes)

sample_data(food)$EnvType <- factor(sample_data(food)$EnvType, 
                                    levels = c("DesLardons", "MerguezVolaille", "BoeufHache", "VeauHache", 
                                               "SaumonFume", "FiletSaumon", "FiletCabillaud", "Crevette"))

## Change behavior of plot composition
p <- plot_composition(food, "Kingdom", "Bacteria", "Phylum", fill = "Phylum", facet_grid = "~EnvType")
ggplotly(p, tooltip = c("fill", "y", "x"))

## Change behavior of plot ordination
dist.uf <- distance(food, method = "unifrac")
p <- plot_ordination(food, ordination = ordinate(food, distance = dist.uf, method = "MDS"), 
                     color = "EnvType")
plot(p)
ggplotly(p + geom_text(aes(label = EnvType), color = "transparent"), tooltip = c("label"))



## Test plot_composition

## test different gradient for plot_dist_as_heatmap
dist.uf <- distance(food, "unifrac")
order <- levels(reorder(sample_names(food), 
                        as.numeric(get_variable(food, "EnvType"))))
grid.arrange(plot_dist_as_heatmap(dist.uf, order, low = "#B1F756", high = "#132B13"), 
             plot_dist_as_heatmap(dist.uf, order, low = "#e4aad7", high = "#d43f77"), 
             plot_dist_as_heatmap(dist.uf, order, low = "#fdc9a1", high = "#e33c25"), 
             plot_dist_as_heatmap(dist.uf, order, low = "#457363", high = "#df3d77"),
             ncol = 2)

## test plot_clust
plot_clust(food, dist.uf, method = "ward.D2", color = "EnvType")
