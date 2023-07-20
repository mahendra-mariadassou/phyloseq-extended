#' tax_spread
#'
#' Spread taxonomy according to last known taxa, to remove unknown and multi-affiliations
#'
#' @param physeq (Required). phyloseq-class.
#' @param pattern string. Pattern matching with taxon to rename. Can use regex, use `|` for separator
#' @param explicit Logical. Should the spreading be explicit (e.g. "Unknown Firmicute species" for the species rank) (TRUE, default) or implicit (e.g. "Firmicutes" for all ranks below Phylum).
#'
#' @importFrom dplyr if_else last
#' @importFrom phyloseq access rank_names tax_table
#' @importFrom stringr regex str_detect
#'
#' @return physeq object with correct name
#' @export
#'
#' @examples
#' library(phyloseq)
#' physeq <- phyloseq(otu_table(matrix(1:4, 2, 2), taxa_are_rows = TRUE),
#' tax_table(matrix(c("Firmicutes", "Firmicutes", "Unknown", "Bacilli", NA, "Lactobacillales"), 2, 3)))
#' tax_spread(physeq) |> tax_table()
#' tax_spread(physeq, explicit = FALSE) |> tax_table()
tax_spread <- function(physeq,
                       pattern = "^NA$|Multi-affiliation|Unknown",
                       explicit = TRUE) {
  if (is.null(access(physeq, "tax_table"))) {
    stop("The tax_spread() function requires that physeq contain a taxonomyTable")
  }
  ranks <- rank_names(physeq)
  ## Spreading function
  .spread <- function(x) {
    last_aff <- x[!is.na(x) & str_detect(x, regex(pattern, ignore_case = TRUE), negate = TRUE)] %>% dplyr::last()
    dplyr::if_else(
      str_detect(x, regex(pattern, ignore_case = TRUE)) | is.na(x),
      if (explicit) paste("Unknown", last_aff, ranks) else last_aff,
      x
    )
  }
  tax_table(physeq) <- tax_table(physeq) %>% apply(MARGIN = 1, .spread) %>% t() %>% `colnames<-`(ranks)
  return(physeq)
}
