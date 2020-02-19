#' Compute local prevalence, abundance and specificity of an OTU in different levels of a group
#'
#' @param physeq phyloseq class object, otu abundances are extracted from this object
#' @param group Either the a single character string matching a variable name in
#'              the corresponding sample_data of `physeq`, or a factor with the same
#'              length as the number of samples in `physeq`.
#' @param format long (default) or wide. Should the result be organised as a long data.frame
#'               with columns "otu", "group", "prevalence", "specificity" and "abundance" or a numeric matrix
#'               with otus in rows and group specificty/group prevalence/group abundance in column.
#' @param rarefy (Logical). Default TRUE. Rarefy samples before computing prevalences.
#'
#' @details Specificity is defined here as local prevalence / global prevalence,
#'          and therefore only deals with presence/absence data. Rarefaction is
#'          therefore recommended.
#'
#' @return A tibble with components
#' * prevalence: observed prevalence in group
#' * specifity : specifity (prevalence / global prevalence) in group
#' * abundance : relative abundance in group
#' * group: factor level
#' * otu: otu
#' or a matrix m of size `ntaxa(physeq)` times 3*N where
#' N is the number of levels (groups) in `group` and
#' * `m[i, g]` is the prevalence of otu i in group g
#' * `m[i, N+g]` is the specificity of otu i to group g
#' * `m[i, 2N+g]` is the relative abundance of otu i in group g
#' @export
#'
#' @importFrom phyloseq sample_data otu_table rarefy_even_depth
#' @importFrom dplyr inner_join as_tibble
#' @importFrom tidyr gather
#'
#' @examples
#' data(food)
#' estimate_prevalence(food, group = "EnvType", rarefy = TRUE)
estimate_prevalence <- function(physeq, group, format = c("long", "wide"), rarefy = TRUE) {
    ## Get grouping factor
    if (!is.null(phyloseq::sample_data(physeq, FALSE))) {
        if (class(group) == "character" & length(group) == 1) {
            x1 <- data.frame(phyloseq::sample_data(physeq))
            if (!group %in% colnames(x1)) {
                stop("group not found among sample variable names.")
            }
            group <- x1[, group]
        }
    }
    if (class(group) != "factor") {
        group <- factor(group)
    }
    ## Rarefy samples
    if (rarefy) physeq <- phyloseq::rarefy_even_depth(physeq, rngseed = 1121983,
                                                      trimOTUs = FALSE, verbose = FALSE)
    ## Construct relative abundances by sample
    tdf <- as(phyloseq::otu_table(physeq), "matrix")
    if (!taxa_are_rows(physeq)) { tdf <- t(tdf) }
    ## abundance
    tdf[, ] <- apply(tdf, 2, function(x) { x / sum(x) })
    abund <- t(rowsum(t(tdf), group, reorder = TRUE)) / matrix(rep(table(group), each = nrow(tdf)),
                                                               nrow = nrow(tdf))
    tdf <- 0 + (tdf > 0)
    ## prevalence
    frac <- t(rowsum(t(tdf), group, reorder = TRUE)) / matrix(rep(table(group), each = nrow(tdf)),
                                                              nrow = nrow(tdf))
    ## specificity
    spec <- t(rowsum(t(tdf), group, reorder = TRUE)) / rowSums(tdf)
    spec[rowSums(tdf) == 0, ] <- 0
    format <- match.arg(format)
    if (format == "wide") {
        ## Pad colnames
        colnames(frac) <- paste("prev", colnames(frac), sep = "_")
        colnames(spec) <- paste("spec", colnames(spec), sep = "_")
        colnames(abund) <- paste("abund", colnames(abund), sep = "_")
        ## cbind prevalence, specificity and abundances
        return(cbind(frac, spec, abund))
    }
    ## Melt and join tables
    dplyr::as_tibble(frac, rownames = "otu") %>% tidyr::gather(key = "group", value = "prevalence", -otu) %>%
        dplyr::inner_join(dplyr::as_tibble(spec, rownames = "otu") %>% tidyr::gather(key = "group", value = "specificity", -otu),
                          by = c("otu", "group")) %>%
        dplyr::inner_join(dplyr::as_tibble(abund, rownames = "otu") %>% tidyr::gather(key = "group", value = "abundance", -otu),
                          by = c("otu", "group"))
}
