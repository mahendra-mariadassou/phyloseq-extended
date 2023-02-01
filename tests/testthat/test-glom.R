test_that("staggered throws an error for overlapping preserved taxa", {
  expect_error(staggered_tax_glom(food, atomic_taxa = c("BS11 gut group", "Bacteroidetes"), taxrank = "Phylum"))
})

test_that("staggered throws an error for when agglomerating rank is too deep compared to preserved taxa", {
  expect_error(staggered_tax_glom(food, atomic_taxa = c("BS11 gut group"), taxrank = "Genus"))
})

test_that("staggered gives expected number of taxa", {
  phy1 <- staggered_tax_glom(food, atomic_taxa = c("BS11 gut group"), taxrank = "Phylum")
  expect_equal(phy1 |> ntaxa(), 12)
  expect_equal(rank_names(phy1), c("Kingdom", "Phylum", "Class", "Order", "Family"))
  phy2 <- staggered_tax_glom(food, atomic_taxa = c("BS11 gut group", "Serratia"), taxrank = "Phylum")
  expect_equal(ntaxa(phy2), 13)
  expect_equal(rank_names(phy2), c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"))
})
