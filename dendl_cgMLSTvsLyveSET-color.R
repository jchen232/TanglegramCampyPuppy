#!/usr/bin/env Rscript
# Authors: Beau Bruce and Weidong Gu
# Modified by Lee Katz and Jess Chen

library("phytools")
library("ape")
library("dendextend")

treefile1 <- read.tree(file = "labpaper2018_cgMLST_newprod.dnd")
treefile2 <- read.tree(file = "out.RAxML_bipartitions.dnd")

outbreak <- read.table('clade_numbers2.txt',
                       sep="\t", header=T, stringsAsFactors=F)

tree1 <- reorder(midpoint.root((treefile1)), order = "cladewise")
tree2 <- reorder(midpoint.root((treefile2)), order = "cladewise")

dend_tree1 <- force.ultrametric(tree1)
dend_tree2 <- force.ultrametric(tree2)

min_length <- 0.000000000000000000001

dend_tree1$edge.length[ dend_tree1$edge.length < min_length ] <- min_length
dend_tree2$edge.length[ dend_tree2$edge.length < min_length ] <- min_length

dend_tree1=(midpoint.root(dend_tree1))
dend_tree2=(midpoint.root(dend_tree2))

clade1 <- match(outbreak$WGS_id[outbreak$Clade == 1],dend_tree1$tip.label)
clade2 <- match(outbreak$WGS_id[outbreak$Clade == 2],dend_tree1$tip.label)

myColors <- c()
myColors[clade1] <- 'red'
myColors[clade2] <- 'blue'


myNewColors <- c(myColors[myColors=="blue"],myColors[myColors=="red"])

print("untangle")
dendl <- dendextend::untangle(as.dendrogram(dend_tree1), 
                     as.dendrogram(dend_tree2), 
                     method = "step2side") 



# Make the branches look nice
dendl %>% set("branches_lwd", 1) %>%
  set("labels_col", "white") -> dendl

print("entanglement...");
myEntanglement <- entanglement(dendl)
cophenetic <- cor.dendlist(dendl, method = "cophenetic")
baker      <- cor.dendlist(dendl, method = "baker")


# Start off the viz
tiff(file="tangle-dendl-mlst-hqsnp-color2.tiff", width = 7, height = 7, unit = "in",res =600 )
tanglegram(dendl,
           main_left='cgMLST',
           main_right='lyve-SET',
           lab.cex=0.3,
           highlight_distinct_edges = FALSE,
           color_lines=myNewColors, 
	     lwd=1,
	     common_subtrees_color_lines = FALSE
           )
myReturn <- dev.off();

