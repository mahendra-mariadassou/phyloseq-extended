library(phyloseq)
food <- import_biom(
  BIOMfilename = "data-raw/chaillou/chaillou.biom",
  treefilename = "data-raw/chaillou/tree.nwk",
  parseFunction = parse_taxonomy_greengenes
)
sample_data(food)$EnvType <- factor(sample_data(food)$EnvType,
  levels = c(
    "DesLardons", "MerguezVolaille", "BoeufHache", "VeauHache",
    "SaumonFume", "FiletSaumon", "FiletCabillaud", "Crevette"
  )
)
usethis::use_data(food)
