
# df should have columns: tail, n_event, n_died (n_event includes n_died)
# If assume_all_died=TRUE, set n_died <- n_event.
calc_km <- function(df, assume_all_died=FALSE) {
    # Kaplan-Meier calculation
    km <- df |>
        dplyr::arrange(tail) |>
        dplyr::collect() #Can't be arrow
    
    if (assume_all_died) {
        km$n_died <- km$n_event
    }
    
    cum_n <- cumsum(km$n_event)
    total_n <- tail(cum_n,1)
    km$active_n <- total_n - dplyr::lag(cum_n,1,0)
    survivors_prop_survived <- 1 - km$n_died / km$active_n
    km$prop_after <- cumprod(survivors_prop_survived)
    km$prop_before <- dplyr::lag(km$prop_after,1,1)
    km$prop_died <- km$prop_before - km$prop_after
    km
}

km_quantile <- function(km, prop) {
    # Average if we exactly hit a proportion
    tail1 <- km$tail[ match(TRUE, km$prop_after <= prop) ]
    tail2 <- km$tail[ match(TRUE, km$prop_after < prop) ]
    (tail1 + tail2) / 2
}

km_complete <- function(km, min_tail=0, max_tail=NULL) {
    if (is.null(max_tail)) 
        max_tail <- max(km$tail)
    km |>
        dplyr::select(tail, n_event, n_died, prop_died) |>
        dplyr::filter(tail >= .env$min_tail, tail <= .env$max_tail) |>
        tidyr::complete(tail=seq(min_tail,max_tail),fill=list(n_event=0, n_died=0, prop_died=0))
}


#' Join reads to sites
#'
#' @param site_pad Reads ending within [position-site_pad, position+site_pad] will be counted to a site.
#'
#' @export
site_reads <- function(reads, sites, site_pad) {
    sites <- dplyr::collect(sites)
    reads <- dplyr::collect(reads)
    
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
#' @param min_tail Tail lengths shorter than this are discarded.
#'
#' @param length_trim This many bases are effectively trimmed from the ends of reads.
#'
#' @export
count_tails <- function(sited_reads, min_tail, length_trim) {
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
        arrow::as_arrow_table() |>
        set_attr("min_tail", min_tail) |>
        set_attr("length_trim", length_trim)
    
    cumulation
}


#' @export
combine_tail_counts <- function(tail_counts_list) {
    tail_counts_list <- ensure_list(tail_counts_list)
    
    result <- tail_counts_list |>
        purrr::map(dplyr::collect) |>
        dplyr::bind_rows() |>
        dplyr::summarise(n_event=sum(n_event), n_died=sum(n_died), .by=c(site, tail))
    
    # TODO: check consistency
    result |>
        set_attr("min_tail", get_attr(tail_counts_list[[1]], "min_tail")) |>
        set_attr("length_trim", get_attr(tail_counts_list[[1]], "length_trim")) |>
        set_attr("max_tail", max(0, result$tail))
}


#' @export
calc_site_stats <- function(tq) {
    sites <- tq$sites
    tail_counts <- combine_tail_counts(tq$samples$tail_counts)
    
    result <- tail_counts |>
        tidyr::nest(.by=site, .key="tail_counts") |>
        dplyr::mutate(
            km=map_progress(tail_counts, "Kaplan-Meier calculation", calc_km),
            n=purrr::map_dbl(tail_counts, \(df) sum(df$n_event)),
            n_died=purrr::map_dbl(tail_counts, \(df) sum(df$n_died)),
            tail10=purrr::map_dbl(km, km_quantile, 0.1),
            tail25=purrr::map_dbl(km, km_quantile, 0.25),
            tail50=purrr::map_dbl(km, km_quantile, 0.5),
            tail75=purrr::map_dbl(km, km_quantile, 0.75),
            tail90=purrr::map_dbl(km, km_quantile, 0.9))
    
    if (!is.null(sites))
        result <- dplyr::left_join(result, sites, by="site", copy=TRUE)
    
    result |>
        set_attr("min_tail", get_attr(tail_counts, "min_tail")) |>
        set_attr("length_trim", get_attr(tail_counts, "length_trim")) |>
        set_attr("max_tail", get_attr(tail_counts, "max_tail"))
}


#' @export
stats_density_mat <- function(stats, min_tail=0, max_tail=NULL, what="prop_died") {
    if (is.null(max_tail))
        max_tail <- max(purrr::map_dbl(stats$km,\(km) max(km$tail)))
    
    completes <- purrr::map(stats$km, \(km) km_complete(km, min_tail=min_tail, max_tail=max_tail))
    mat <- do.call(rbind, purrr::map(completes, what))
    colnames(mat) <- seq(min_tail, max_tail, by=1)
    rownames(mat) <- stats$site
    mat
}
