#' tax_spread
#'
#' Spread taxonomy according to last known taxa, to remove unknown and multi-affiliations
#'
#' @param physeq (Required). phyloseq-class.
#' @param pattern string. Pattern matching with taxon to rename. Can use regex, use `|` for separator
#'
#' @importFrom phyloseq tax_table
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom dplyr mutate group_by select last if_else
#' @importFrom stringr str_detect regex
#'
#' @return physeq object with correct name
#' @export
#'
#' @examples
#' physeq <- phyloseq(otu_table(matrix(1:4, 2, 2), taxa_are_rows = TRUE),
#' tax_table(matrix(c("Firmicutes", "Firmicutes", "Unknown", "Bacilli", NA, "Lactobacillales"), 2, 3)))
#' tax_spread(physeq)
tax_spread <- function(physeq,
                       pattern = "^NA$|Multi-affiliation|Unknown") {
  if (is.null(access(physeq, "tax_table"))) {
    stop("The tax_spread() function requires that physeq contain a taxonomyTable")
  }
  ranks <- rank_names(physeq)
  ## Spreading function
  .spread <- function(x) {
    last_aff <- x[!is.na(x) & str_detect(x, regex(pattern, ignore_case = TRUE), negate = TRUE)] %>% dplyr::last()
    dplyr::if_else(
      str_detect(x, regex(pattern, ignore_case = TRUE)) | is.na(x),
      paste("Unknown", last_aff, ranks),
      x
    )
  }
  tax_table(physeq) <- tax_table(physeq) %>% as('matrix') %>% apply(MARGIN = 1, .spread) %>% t() %>% `colnames<-`(ranks)
  return(physeq)
}
