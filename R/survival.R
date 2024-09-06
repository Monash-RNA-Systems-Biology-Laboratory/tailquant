
# df should have columns: tail, n_event, n_died (n_event includes n_died)
calc_km <- function(df) {
    # Kaplan-Meier calculation
    km <- df |>
        dplyr::arrange(tail) |>
        dplyr::collect() #Can't be arrow
    cum_n <- cumsum(km$n_event)
    total_n <- tail(cum_n,1)
    km$active_n <- total_n - dplyr::lag(cum_n,1,0)
    survivors_prop_survived <- 1 - km$n_died / km$active_n
    km$prop_after <- cumprod(survivors_prop_survived)
    km$prop_before <- dplyr::lag(km$prop_after,1,1)
    km$prop_died <- km$prop_before - km$prop_after
    
    #cat("Median",km$tail[ match(TRUE, km$prop_after <= 0.5) ],"\n")
    km
}

km_quantile <- function(km, prop) {
    km$tail[ match(TRUE, km$prop_after <= prop) ]
}


#' Join reads to sites
#'
#' @export
site_reads <- function(reads, sites, site_pad=10) {
    # Associate reads with sites
    hits <- tailquant:::find_overlaps(
        sites$chr, sites$pos-site_pad, sites$pos+site_pad, sites$strand,
        reads$chr, reads$pos, reads$pos, reads$strand)
    
    dplyr::tibble(site=sites$site[hits$index1], reads[hits$index2,]) |>
        arrow::as_arrow_table() |>
        dplyr::arrange(site, tail) |>
        arrow::as_arrow_table()
}

#' Count tail lengths
#'
#' Count tail lengths observed at each site, and how many are considered censored.
#'
#' @export
count_tails <- function(sited_reads, min_tail=19, length_trim=5) {
    cumulation <- sited_reads |>
        arrow::as_arrow_table() |> #Faster
        dplyr::transmute(
            site = site,
            length = length-.env$length_trim,
            died = genomic_bases+tail <= length,
            tail = ifelse(died, tail, length-genomic_bases)) |>
        dplyr::filter(tail >= .env$min_tail) |>
        dplyr::summarise(
            n_event = dplyr::n(),
            n_died = sum(died),
            .by = c(site, tail)) |>
        dplyr::arrange(site, tail) |>
        arrow::as_arrow_table()
    
    cumulation
}