#' tax_spread
#'
#' Spread taxonomy according to last knowing taxa
#'
#' @param physeq (Required). phyloseq-class.
#' @param pattern string. Pattern matching with taxon to rename. Can use regex, use `|` for separator
#'
#' @importFrom phyloseq tax_table
#' @importFrom tibble rownames_to_column column_to_rownames
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom dplyr mutate group_by select
#' @importFrom stringr str_detect regex
#'
#' @return physeq object with correct name
#' @export
#'
#' @examples
#' data(food)
#' tax_spread(food)
tax_spread <- function(physeq,
                       pattern = "^NA$|Multi-affiliation|Unknown") {
  if (is.null(access(physeq, "tax_table"))) {
    stop("The tax_spread() function requires that physeq contain a taxonomyTable")
  }

  tax_table(physeq) <- tax_table(physeq) %>%
    as.data.frame() %>%
    rownames_to_column(var = "OTU") %>%
    pivot_longer(cols = -OTU, names_to = "rank", values_to = "taxon") %>%
    mutate(new_taxon = ifelse(str_detect(taxon, regex(pattern, ignore_case = TRUE)), NA, taxon)) %>%
    group_by(OTU) %>%
    mutate(new_taxon = ifelse(is.na(new_taxon),
      paste("Unknown", last(na.omit(new_taxon)), rank),
      new_taxon
    )) %>%
    select(-taxon) %>%
    pivot_wider(names_from = rank, values_from = new_taxon) %>%
    column_to_rownames("OTU") %>%
    as.matrix() %>%
    tax_table()
  return(physeq)
}
