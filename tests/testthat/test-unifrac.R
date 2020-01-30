library(ape)

# Toy example
u_counts <- c(1, 0, 0, 4, 1, 2, 3, 0)
v_counts <- c(0, 1, 1, 6, 0, 1, 0, 0)
table <- cbind(u_counts, v_counts)
dimnames(table) <- list(paste0("OTU", 1:8), c("u", "v"))
tree = ape::read.tree(text="(((((OTU1:0.5,OTU2:0.5):0.5,OTU3:1.0):1.0):0.0,(OTU4:0.75,(OTU5:0.5,((OTU6:0.33,OTU7:0.62):0.5,OTU8:0.5):0.5):0.5):1.25):0.0)root;")
physeq = phyloseq::phyloseq(phyloseq::otu_table(table, taxa_are_rows = TRUE),
                            phyloseq::phy_tree(tree))

test_that("UniFrac fails when the tree is missing, unrooted or lacks branch lengthes", {
  ## Missing root
  physeq@phy_tree <- unroot(unroot(unroot(tree)))
  expect_error(UniFrac(physeq), "Tree is not rooted. Make sure your tree is rooted before attempting UniFrac calculation. See ?ape::root",
               fixed = TRUE)
  ## Missing branch lengthes
  physeq@phy_tree$edge.length <- NULL
  expect_error(UniFrac(physeq), "Tree has no branch lengths. See tree$edge.length. Cannot compute UniFrac without branch lengths",
               fixed = TRUE)
  ## Missing tree
  physeq@phy_tree <- NULL
  expect_error(UniFrac(physeq), "UniFrac distances require a phylogenetic tree.",
               fixed = TRUE)
})

test_that("UniFrac returns a distance matrix", {
  expect_true(inherits(UniFrac(physeq), "dist"))
})

test_that("Unweighted UniFrac works", {
  expect_equal(UniFrac(physeq, weighted = FALSE, normalized = FALSE)[1],
               3.12)
  expect_equal(UniFrac(physeq, weighted = FALSE, normalized = TRUE)[1],
               3.12 / 8.45)
})

test_that("Unweighted UniFrac works on non-binary trees", {
  tree <- "(((OTU1:0.5,OTU2:0.5):0.5,OTU3:1):1,(OTU4:0.75,(OTU5:0.5,((OTU6:0.33,OTU7:0.62,OTU8:0.5):0.5):0.5):0.5):1.25);"
  physeq@phy_tree <- read.tree(text = tree)
  expect_equal(UniFrac(physeq, weighted = FALSE, normalized = FALSE)[1],
               3.12)
  expect_equal(UniFrac(physeq, weighted = FALSE, normalized = TRUE)[1],
               3.12 / 8.45)
})

test_that("Weighted UniFrac works", {
  expect_equal(UniFrac(physeq, weighted = TRUE, normalized = FALSE)[1],
               1.54343434343434)
  expect_equal(UniFrac(physeq, weighted = TRUE, normalized = TRUE)[1],
               1.54343434343434 / 4.71272727272727)
})
