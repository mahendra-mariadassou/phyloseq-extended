#' @title Merge samples based on a grouping variable given either as a sample variable or a factor.
#'
#' @description The purpose of this method is to merge/agglomerate the sample indices of a phyloseq object according to a categorical variable contained in a sample_data or a provided factor. Unlike the original \code{\link{merge_samples}}, sample metadata are not coerced to factors. Instead, the metadata from the first sample in each group is kept and the samples are renamed according to the levels of the factors.
#'
#' @inheritParams phyloseq::merge_samples
#' @param fun Either "sum" (default) or "mean". Controls how abundances are merged within each factor level.
#' @param update.names Logical. If `TRUE` (default) merged samples are named according to the levels of the grouping variable. Else, keep the name of the first sample in each group.
#'
#' @return A phyloseq object that has had its sample indices merged according to the factor indicated by the group argument. The output class matches x.
#'
#' @seealso \code{\link{merge_samples}}
#'
#' @export
#'
#' @examples
#' data(food)
#' mg <- merge_group(food, "EnvType")
#' mg
#' ## Note that sample data are preserved
#' phyloseq::sample_data(mg)
#' @importFrom methods as
#' @importFrom phyloseq get_variable otu_table sample_data sample_names sample_variables taxa_are_rows
merge_group <- function(physeq, group, fun = c("sum", "mean"), update.names = TRUE) {
  fun <- match.arg(fun)

  # Build grouping factor
  if (!is.null(sample_data(physeq, FALSE))) {
    if (class(group) == "character" & length(group) == 1) {
      if (!group %in% sample_variables(physeq)) {
        stop("group not found among sample variable names.")
      }
      group <- get_variable(physeq, group)
    }
  }
  if (class(group) != "factor") {
    group <- factor(group)
  }
  # Group abundances
  cdf <- as(otu_table(physeq), "matrix")
  if (taxa_are_rows(physeq)) cdf <- t(cdf)
  cdf_merged <- rowsum(cdf, group, reorder = FALSE)
  if (fun == "mean") {
    scaling <- table(group)[rownames(cdf_merged)]
    cdf_merged <- cdf_merged / matrix(scaling[row(cdf_merged)], ncol = ncol(cdf_merged))
  }
  ## Rename merged samples according to first sample in group
  rownames(cdf_merged) <- sample_names(physeq)[!duplicated(group)]

  ## Update otu_table
  otu_table(physeq) <- otu_table(cdf_merged, taxa_are_rows = FALSE)

  ## Update sample names
  if (update.names) {
    sample_names(physeq) <- unique(group)
  }

  physeq
}

#' extract in order all ranks above above a given rank
#'
#' @param physeq Required. \code{\link{phyloseq-class}} object with a taxonomic table.
#' @param ranks  Required. One or more ranks whose ancestry are search (e.g. "Phylum", "Species", etc). Values not encountered in `physeq` are silently ignored
#'
#' @return a vector of rank names
#' @export
#'
#' @examples
#' data(food)
#' find_upper_ranks(food, "Phylum") ## c("Kingdom", "Phylum")
#' find_upper_ranks(food, c("Class", "Genus")) ## c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
#' @importFrom phyloseq rank_names
find_upper_ranks <- function(physeq, ranks) {
  rank_numbers <- match(ranks, rank_names(physeq))
  if (all(is.na(rank_numbers))) {
    stop("Bad ranks. Must be among the values of rank_names(physeq)")
  }
  rank_names(physeq)[1:max(rank_numbers, na.rm = TRUE)]
}

#' Fast alternative to [tax_glom()]
#'
#' @param physeq Required. \code{\link{phyloseq-class}} object
#' @param taxrank A character string specifying the taxonomic level that you want to agglomerate over. Should be among the results of \code{rank_names(physeq)}. The default value is \code{rank_names(physeq)[1]}, which may agglomerate too broadly for a given experiment. You are strongly encouraged to try different values for this argument.
#' @param bad_empty (Optional). Character vector. Default to empty white spaces and tabs. Defines the bad/empty values that should be replaced by "Unknown".
#'
#' @return A taxonomically-agglomerated \code{\link{phyloseq-class}} object.
#' @export
#'
#' @details fast_tax_glom differs from [tax_glom()] in three important ways:
#' * It is based on dplyr function and thus much faster
#' * It does not preserve the phy_tree slot
#' * It only preserves taxonomic ranks up to \code{taxrank} and discard all other ones.
#' The "archetype" (OTU name) and sequence for each group are chosen as those from the most abundant taxa
#' within that group (for compatibility with [tax_glom()])
#'
#' @seealso [tax_glom()], [merge_taxa()]
#' @importFrom dplyr across all_of arrange as_tibble cur_group_id desc group_by mutate select slice
#' @importFrom methods as
#' @importFrom phyloseq access otu_table phyloseq rank_names sample_data tax_table taxa_are_rows taxa_names taxa_sums
#' @importFrom tibble column_to_rownames
#' @examples
#' data(food)
#' fast_tax_glom(food, "Species")
fast_tax_glom <- function(physeq, taxrank = rank_names(physeq)[1], bad_empty = c("", " ", "\t")) {
  ranks <- find_upper_ranks(physeq, taxrank)
  ## change NA and empty to Unknown before merging
  tax <- as(tax_table(physeq), "matrix")
  tax[is.na(tax) | tax %in% bad_empty] <- "Unknown"
  ## create groups
  tax <- tax[ , ranks, drop = FALSE] %>%
    dplyr::as_tibble(tax, .name_repair = "minimal") %>%
    dplyr::mutate(Abundance = taxa_sums(physeq),
                  archetype = taxa_names(physeq)) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(ranks))) %>%
    dplyr::mutate(group = cur_group_id())
  ## create new_taxonomy
  new_tax <- tax %>%
    dplyr::arrange(desc(Abundance)) %>%
    dplyr::slice(1) %>% dplyr::arrange(group) %>%
    dplyr::select(-group, -Abundance) %>%
    tibble::column_to_rownames(var = "archetype") %>% as.matrix()
  ## create new count table
  otutab <- otu_table(physeq)
  if (!taxa_are_rows(physeq)) otutab  <- t(otutab)
  otutab <- rowsum(otutab, group = tax$group, reorder = TRUE)
  rownames(otutab) <- rownames(new_tax)
  ## create new refseq
  seqs <- access(physeq, "refseq")
  if (!is.null(seqs)) seqs <- seqs[rownames(new_tax)]
  ## return merged phyloseq
  phyloseq(sample_data(physeq),
           tax_table(new_tax),
           otu_table(otutab, taxa_are_rows = TRUE),
           seqs
  )
}

#' Staggered taxonomic agglomeration
#'
#' @param physeq Required. \code{\link{phyloseq-class}} object
#' @param atomic_taxa A character string specificying the taxa that won't be agglomerated. Note that taxa can be mixed (i.e. Species and Phylum) but are not allowed to overlap.
#' @param taxrank A character string specifying the taxonomic level that you want to agglomerate over (except for the taxa in `atomic_taxa`).
#'
#' @return A taxonomically-agglomerated \code{\link{phyloseq-class}} object.
#' @export
#'
#' @details staggered_tax_glom differs from [fast_tax_glom()] by preserving some taxa during the taxonomic agglomeration phase
#'
#' @importFrom dplyr if_else
#' @importFrom phyloseq merge_phyloseq ntaxa prune_taxa rank_names tax_table
#' @importFrom stringr str_remove
#'
#' @examples
#' data(food)
#' staggered_tax_glom(food, atomic_taxa = "BS11 gut group", taxrank = "Phylum") |> phyloseq::tax_table()
#' staggered_tax_glom(food, atomic_taxa = c("BS11 gut group", "Serratia"), taxrank = "Phylum")  |> phyloseq::tax_table()
staggered_tax_glom <- function(physeq, atomic_taxa, taxrank) {
  taxa_index <- matrix(FALSE, nrow = ntaxa(physeq), ncol = length(rank_names(physeq)))
  taxa_index[] <- tax_table(physeq) %in% atomic_taxa
  if (any(rowSums(taxa_index) > 1)) {
    stop("Overlap between preserved taxa detected. Make sure to specify which taxa to conserve.")
  }
  preserved_ranks <- which(colSums(taxa_index) > 0)
  preserved_taxa <- rowSums(taxa_index) == 1
  min_rank <- min(preserved_ranks); max_rank <- max(preserved_ranks)
  glom_rank <- which(rank_names(physeq) == taxrank)
  if (min_rank < glom_rank) {
    stop("The agglomeration rank conflicts with some of the taxa. Choose a higher rank.")
  }
  ## Find rightmost rank and set all further values to NA
  preserved_taxa_coord <- which(taxa_index, arr.ind = TRUE)
  rightmost_rank <- rep(glom_rank, ntaxa(physeq))
  rightmost_rank[preserved_taxa_coord[, 1]] <- preserved_taxa_coord[, 2]
  tax_table(physeq)[col(tax_table(physeq)) > rightmost_rank] <- NA_character_
  ## Extract most accurate taxonomy for each taxa and replace
  # rightmost_name <- as(tax_table(physeq), "matrix")[matrix(c(1:ntaxa(physeq), rightmost_rank), ncol = 2)]
  # replacement_mask <- col(tax_table(physeq)) > rightmost_rank
  # tax_table(physeq)[replacement_mask] <- rightmost_name[row(tax_table(physeq))][replacement_mask]
  other_taxa <- prune_taxa(!preserved_taxa, physeq) |>
    fast_tax_glom(taxrank)
  preserved_taxa <- prune_taxa(preserved_taxa, physeq) |>
    fast_tax_glom(rank_names(physeq)[max_rank])
  ## pad ranks for other taxonomic levels
  # if (min_rank == max_rank) {
  ## manage rank names in other_taxa
  preserved_glom <- tax_table(preserved_taxa)[, glom_rank] %>% unique() %>% as.character()
  srank <- tax_table(other_taxa)[, glom_rank] %>% as.character()
  tax_table(other_taxa)[, glom_rank] <- if_else(
    srank %in% preserved_glom,
    paste("Other", srank),
    srank
  )
  # }
  res <- merge_phyloseq(other_taxa, preserved_taxa) %>%
    tax_spread(explicit = FALSE)
  tax_table(res)[ , glom_rank] <- stringr::str_remove(tax_table(res)[ , glom_rank], "^Other ")
  res
}

#' Compute core microbiome
#'
#' @param physeq Required. \code{\link{phyloseq-class}} object
#' @param group Optional. Grouping variable, used to compute the core microbiome in a groupwise fashion. Either a single character string matching a variable name in the corresponding `sample_data` of `physeq`, or a factor with the same length as the number of samples in `physeq`. If `NULL` all samples belong to the same group.
#' @param ab_threshold Numeric. The minimum (relative) abundance for a taxa to be considered present. Can be a vector
#' @param prev_threshold Numeric. The minimum prevalence (across samples in the group) for a taxa to be considered core. Defaults to 0.5
#'
#' @return a tibble with components
#' * `OTU` taxa name
#' * `group` grouping level (only if group is speficied)
#' * `raw_prevalence` prevalence in the group (using no abundance threshold)
#' * `raw_abundance` average relative abundance in the group
#' * `ab_threshold` abundance threshold used to compute coreness
#' * `coreness` coreness value (prevalence above the abundance threshold) in the group
#' * `is_core` logical indicating whether the taxa is core in the group
#'
#' @export
#'
#' @seealso
#' [estimate_prevalence]
#'
#' @examples
#' extract_core(food, group = "EnvType")
#' @importFrom dplyr as_tibble filter group_by mutate select summarise ungroup
#' @importFrom methods as
#' @importFrom phyloseq get_variable nsamples otu_table sample_data sample_variables taxa_are_rows transform_sample_counts
#' @importFrom tidyr crossing pivot_longer
extract_core <- function(physeq, group = NULL, ab_threshold = 0, prev_threshold = 0.5) {

  # Build grouping factor
  if (is.null(group)) {
    group <- rep(1, phyloseq::nsamples(physeq))
  }
  if (!is.null(phyloseq::sample_data(physeq, FALSE))) {
    if (class(group) == "character" & length(group) == 1) {
      if (!group %in% phyloseq::sample_variables(physeq)) {
        stop("group not found among sample variable names.")
      }
      group <- phyloseq::get_variable(physeq, group)
    }
  }
  if (class(group) != "factor") {
    group <- factor(group)
  }

  # Melt count table and add grouping information
  cdf <- physeq %>%
    phyloseq::transform_sample_counts(function(x) { x / sum(x)}) %>%
    phyloseq::otu_table() %>%
    as("matrix")
  if (taxa_are_rows(physeq)) cdf <- t(cdf)
  cdf %>% dplyr::as_tibble(rownames = "Sample") %>%
    dplyr::mutate(group = group) %>%
    tidyr::pivot_longer(cols = -c(Sample, group), names_to = "OTU", values_to = "freq") %>%
    tidyr::crossing(ab_threshold = ab_threshold) %>%
    dplyr::group_by(OTU, group, ab_threshold) %>%
    dplyr::summarise(prevalence = mean(freq > 0),
                     abundance  = mean(freq),
                     coreness   = mean(freq > ab_threshold)
                     ) %>%
    dplyr::mutate(is_core = coreness >= prev_threshold) %>%
    dplyr::group_by(OTU) %>%
    dplyr::mutate(any_core = any(is_core)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(any_core) %>% dplyr::select(-any_core)
}

#' clr transformation
#'
#' @param data a integer vector or numeric vector of non-negative values
#'
#' @return the clr-transformed vector
#'
#'
#' @examples
#' clr(c(1, 2, 4))
#'
#' @export
#'
clr <- function(x) {
  if (any(x <= 0)) stop("clr is not defined for negative data")
  log_x <- log(x)
  log_x - mean(log_x)
}

#' mclr transformation
#'
#' @param data a integer vector or numeric vector of non-negative values
#' @param c positive constant used to separate the minimum value (after clr transform) and the zeros. The default value is the one recommended by the authors.
#'
#' @return the modified clr-transformed vector as defined in Yoon, Gaynanova and Müeller (2019)
#'
#' @references Yoon, Gaynanova and Müeller (2019), Frontiers in Genetics, Microbial Networks in SPRING - Semi-parametric Rank-Based Correlation and Partial Correlation Estimation for Quantitative Microbiome Data. doi:10.3389/fgene.2019.00516
#'
#' @details The mclr transform of a vector x is given by
#' \deqn{mclr(x) = \left(0, \frac{log(x_i)}{g(x)} + \varepsilon, \dots \right)}{%
#' mclr(x) = (0, log(x_i)/g(x) + epsilon, ...)}
#' where \eqn{g(x)} is the geometric mean of x (excluding 0s) and all non-zeros values are shifted by epsilon to preserve the original ordering of the data.
#' @examples
#' counts <- matrix(c(0, 1, 1, 1, 2, 0), nrow = 2)
#' mclr(counts)
#'
#' @export
#'
mclr <- function(data, c = 1) {
  zero_mask <- data == 0
  log_data <- log(data)
  log_data[zero_mask] <- NA_real_
  clr_data <- log_data - rowMeans(log_data, na.rm = TRUE)
  shift <- min(clr_data, na.rm = TRUE)
  mclr_data <- clr_data + abs(shift) + c
  mclr_data[zero_mask] <- 0
  mclr_data
}
