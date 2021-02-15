#' Fit Sloan et al. (2006) Neutral Community Model and return results
#'
#' @param physeq phyloseq class object
#' @inheritParams ncm_loglik
#'
#' @return A list of two tibbles and a list
#' - a tibble called `ncm_prediction` containing the predicted NCM curve for the data, with 95% confidence interval for the occupancy estimate. Contains 4 columns:
#'   - `freq` (frequency in the regional pool)
#'   - `occ` (predicted occupancy)
#'   - `occ_[lower|upper]` (confidence interval for the occupancy based on the Wilson confidence interval around the predicted occupancy)
#' - a tibble called `observations` with one row per species containing its predicted frequency, occupancy (with confidence interval) and the observed occupancy under the fitted NCM. Contains 6 columns:
#'   - `otu` (otu name)
#'   - `freq` (estimated frequency in the regional pool)
#'   - `occ` (predicted occupancy)
#'   - `occ_[lower|upper]` (confidence interval for the occupancy based on the Wilson confidence interval around the predicted occupancy)
#'   - `occ_obs` (observed occupancy).
#' - a list called `fit` with the estimated value for the migration parameter m and the fitted ncm curve for the data.
#' @export
#'
#' @examples
#' data(food)
#' model <- ncm_fit(food)
ncm_fit <- function(physeq, threshold, absolute.threshold = 1) {
    abundances <- otu_table(physeq)
    ## Extract species count
    if (!taxa_are_rows(physeq)) {
        abundances <- t(abundances)
    }
    n_samples <- phyloseq::nsamples(physeq)
    depths <- phyloseq::sample_sums(physeq)
    if (!is.null(absolute.threshold)) {
     threshold <- absolute.threshold / depths
    }
    ## Fit NCM model to the data and extract estimated migration rate
    m_hat <- ncm_mle_fit(abundances, threshold, absolute.threshold)$maximum

    ## Compute predicted abundances in regional pool given migration rate m_hat
    freq <- local_community(abundances, m_hat)

    ncm_curve <- ncm_curve(threshold = threshold, N = depths, m = m_hat)
    ## Compute frequency occupancy curve on a regular grid
    freq_occ_ncm <- dplyr::tibble(
        freq     = 10^(-seq(0, 6, length.out = 100)),
        occ      = ncm_curve(freq),
        ci       = wilson_conf_int(p = occ, n = n_samples)
    ) %>%
        dplyr::mutate(occ_lower = ci[, "lower"], occ_upper = ci[, "upper"]) %>% dplyr::select(-ci)

    freq_occ_obs <- dplyr::tibble(
        otu     = names(freq),
        freq    = freq,
        occ     = ncm_curve(freq),
        occ_obs = apply(abundances, 1, function(x) { mean(x > threshold) }),
        ci      = wilson_conf_int(p = occ, n = n_samples)
    ) %>%
        dplyr::mutate(occ_lower = ci[, "lower"], occ_upper = ci[, "upper"],
                      neutral   = dplyr::if_else(occ_obs <= occ_upper & occ_obs >= occ_lower, TRUE, FALSE)) %>%
        dplyr::select(-ci)

    return(list(ncm_prediction = freq_occ_ncm,
                observations = freq_occ_obs,
                fit = list(m = m_hat, ncm_curve = ncm_curve)))

}

#' Plot a ncm_fit result
#'
#' @param ncm.fit scnm fit, as produced by [ncm_fit()]
#'
#' @return a ggplot object
#' @export
#'
#' @examples
#' data(food)
#' model <- ncm_fit(food)
#' plot_ncm(model)
plot_ncm <- function(ncm.fit) {
    ggplot(data = ncm.fit$observations, aes(x = freq, y = occ)) +
        geom_ribbon(data = ncm.fit$ncm_prediction, aes(ymin = occ_lower, ymax = occ_upper), alpha = 0.1) +
        geom_line(data = ncm.fit$ncm_prediction, group = 1) +
        geom_line(data = ncm.fit$ncm_prediction, aes(y = occ_lower), linetype = 2, group = 1) +
        geom_line(data = ncm.fit$ncm_prediction, aes(y = occ_upper), linetype = 2, group = 1) +
        geom_point(aes(y = occ_obs, color = neutral), alpha = 0.5) +
        scale_x_log10() +
        NULL
}

## fitting function for neutral community model of Sloan et al. (2006).
## Assumes equivalent species, migration rate m (to fit) and an array of relative abundances.

## log-likelihood of the occupancy data under a NCM for a given migration rate m and detection threshold. The local abundances are estimated from the data using a BLUE.
#' Compute the log-likelihood of community data under Sloan et al. (2006) Neutral Community Model
#'
#' @param abundances Required. Community data (count data matrix)
#' @param migration  Numeric. Migration rate (in [0, 1]) of the NCM
#' @param threshold  Numeric. Detection threshold (in [0, 1]): minimum relative abundance required to detect a species in a sample. Ignored if `absolute_threshold` is not `NULL`
#' @param absolute.threshold Integer. Minimum count required to detect a species in a sample.
#'
#' @return A numeric value, the log-likelihood of the observed data under the NCM model.
#' @export
#'
#' @references Sloan, W.T. and Lunn, M. and Woodcock, S. and Head, I.M. and Nee, S. and Curtis, T.P. (2006). Quantifying the roles of immigration and chance in shaping prokaryote community structure. Environmental Microbiology. 8:732–740. doi:10.1111/j.1462-2920.2005.00956.x.
#'
#' @examples
#' data(food)
#' ncm_loglik(otu_table(food), absolute.threshold = 1)
ncm_loglik <- function(abundances, migration = 0.1, threshold = 0.01, absolute.threshold = NULL) {
    Nt <- colSums(abundances)
    ## regional pool species composition
    local.abundances <- local_community(abundances, migration)
    ## compute the log-probability of appearance at a given threshold of an OTU given Nt,
    ## relative abundances in the local pool and the detection threshold
    .local <- function(pi, N, lower.tail = FALSE) {
        ## Parameters of the beta distributions
        alpha <-  migration*N*pi
        beta <-  migration*N*(1-pi)
        if (is.null(absolute.threshold)) {
            pbeta(threshold, shape1 = alpha, shape2 = beta,
                  lower.tail = lower.tail, log.p = TRUE)
        } else {
            pbeta(absolute.threshold/N, shape1 = alpha, shape2 = beta,
                  lower.tail = lower.tail, log.p = TRUE)
        }
    }
    ## Matrix of log-probability of appearance of each species in each community
    lprobs.presence <- outer(local.abundances, Nt, .local)
    lprobs.absence <- outer(local.abundances, Nt, .local, lower.tail = TRUE)
    ## Log-likelihood of the occupancy data
    loglik <- sum(ifelse(abundances > 0, lprobs.presence, lprobs.absence))
    return(loglik)
}

## MLE of migration rate m in NCM
ncm_mle_fit <- function(abundances, threshold, absolute.threshold = NULL) {
    optim.function <- function(m) { ncm_loglik(abundances, m, threshold, absolute.threshold) }
    solution <- optimize(optim.function, interval = c(1e-10, 1), maximum = TRUE)
    return(solution)
}



## Probability of occupancy of a site with total count `N` and detection threshold `threshold`
## by a species with abundance `pi` in the regional pool under a NCM with migration rate `m`. Vectorized over
## N and threshold if there are multiple sites (in which the function returns the mean across sites)
ncm_curve <- function(threshold, N, m) {
    if (length(threshold) == 1) {
        if (threshold >= 1) {
            threshold <- threshold / N
        } else {
            threshold <- rep(threshold, length(N))
        }
    }
    curve <- function(pi) {
        proba_sample <- matrix(NA, nrow = length(pi), ncol = length(N))
        for (i in seq_along(N)) {
            proba_sample[, i] <- pbeta(threshold[i],
                                     shape1 = N[i]*m*pi,
                                     shape2 = N[i]*m*(1-pi),
                                     lower.tail = FALSE)
        }
        return(rowMeans(proba_sample))
    }
    return(curve)
}

## Wilson confidence interval for mean of sum of independent bernoulli (with different parameters)
wilson_conf_int <- function(p, alpha  = 0.05, n) {
    z <- qnorm(1 - alpha/2)
    center <- (p + z^2 / (2*n)) / (1 + z^2 / n)
    dev    <- z / (1 + z^2 / n) * sqrt(p*(1-p)/n + z^2/(4*n^2))
    return(cbind(lower = center - dev, upper = center + dev))
}

## Estimate regional pool species composition using BLUE estimates under a NCM with migration rate m
local_community <- function(abundances, migration = 0.1) {
    Nt <- colSums(abundances)
    weights <- (Nt * migration + 1)
    weights <- weights / ( Nt * sum(weights) )
    ## estimates of local abundances, using BLUE estimates of the proportions
    weights.matrix <- matrix( rep(weights, each = nrow(abundances)) , ncol = ncol (abundances))
    local.abundances <- rowSums( abundances * weights.matrix ) ## pi
    return(local.abundances)
}
