
#' Kaplan-Meier curve estimation
#'
#' @details
#' df should have columns tail, n_event, n_died (n_event includes n_died).
#'
#' If assume_all_died=TRUE, we set n_died = n_event.
#'
#' @export
calc_km <- function(df, assume_all_died=FALSE) {
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
    
    # Greenwood's formula for the standard error of prop_after
    # Ensure zero rather than NaN for zero survival
    #
    # Note: This turned out to be not terribly useful, since it is zero for proportions of 0 and 1, so I've removed the calculation for now.
    #
    #km$se <- sqrt(ifelse(
    #    km$prop_after == 0, 0, 
    #    km$prop_after * cumsum(km$n_died / (km$active_n*(km$active_n-km$n_died)))))
    
    # In the case of no censoring this matches the binomial variance estimate:
    #   p[i]*(1-p[i])/r[1]
    #
    # Let
    #   r[i] = active at timestep i
    #   d[i] = died at timestep i
    #
    # p[i] = r[i+1] / r[1]
    # 
    # With no censoring, d[i] = r[i]-r[i+1]
    #
    # The summation part is:
    #
    # sum( d[i]/(r[i]*r[i+1]) )
    # = sum( (r[i]-r[i+1]) / (r[i]*r[i+1]) )
    # = sum( 1/r[i+1] - 1/r[i] )
    # = 1/r[i+1] - 1/r[1]
    # = (r[1]/r[i+1)) * (1/r[1]) * (1 - r[i+1]/r[1])
    # = 1/p[i] * 1/r[1] * (1-p[i])
    #
    # So the final variance is
    #   p[i]*p[i]*sum = p[i] * (1-p[i]) / r[1]
    
    km
}

#' Quantile from a Kaplan-Meier curve
#'
#' @export
km_quantile <- function(km, prop) {
    # Average if we exactly hit a proportion
    tail1 <- km$tail[ match(TRUE, km$prop_after <= prop) ]
    tail2 <- km$tail[ match(TRUE, km$prop_after < prop) ]
    (tail1 + tail2) / 2
}

# Quantile on some number of standard errors from the Kaplan-Meier curve
#
# @export
#km_quantile_bound <- function(km, prop, z) {
#    km$prop_after <- km$prop_after + z * km$se
#    km_quantile(km, prop)
#}

#' Survival proportion for a specific tail length
#'
#' @export
km_at <- function(km, tail) {
    idx <- match(TRUE, km$tail >= tail) - 1
    if (is.na(idx))
        idx <- nrow(km)
    
    if (idx == 0)
        1
    else
        km$prop_after[idx]
    
    #if (idx == 0)
    #    tibble(prop=1,se=0)
    #else
    #    tibble(prop=km$prop_after[idx],se=km$se[idx])
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
#' @description
#' If a read is near multiple sites, the nearest is chosen.
#'
#' @param site_pad Reads ending within [position-site_pad, position+site_pad] will be counted to a site and used for tail estimation.
#'
#' @param site_upstrand Reads ending at most this far upstrand will also be included in UMI counts.
#'
#' @export
site_reads_into <- function(dest_filename, reads_filename, sites, site_pad, site_upstrand) {
    sites <- dplyr::collect(sites)
    
    has_umis <- "umi" %in% colnames(arrow::open_dataset(reads_filename))
    
    yield <- local_write_parquet(dest_filename)
    template <- dplyr::tibble(
        site=character(0),
        close_to_site=logical(0),
        chr=character(0),
        pos=integer(0),
        strand=integer(0),
        num_hits=integer(0),
        length=numeric(0),
        tail_start=numeric(0),
        tail=numeric(0))
    if (has_umis) {
        template$umi <- character(0)
    }
    yield(template)
    
    scan_parquet(reads_filename, \(reads) {
        reads <- dplyr::collect(reads)
        
        # Associate reads with sites
        left <- ifelse(sites$strand < 0,  site_pad,      site_upstrand)
        right <- ifelse(sites$strand < 0, site_upstrand, site_pad)
        hits <- tailquant:::find_overlaps(
            sites$chr, sites$pos-left, sites$pos+right, sites$strand,
            reads$chr, reads$pos, reads$pos, reads$strand)
        
        # If multiple hits, choose nearest to actual pos
        hits$offset <- abs(reads$pos[hits$index2] - sites$pos[hits$index1])
        #This was unexpectedly slow:
        #hits <- dplyr::slice_min(hits, offset, n=1, by=index2, with_ties=FALSE)
        #So:
        if (nrow(hits) > 0) { # Avoid spurious dplyr warning
            hits <- dplyr::summarize(hits, index1=index1[which.min(offset)], .by=c(index2))
        }
        
        hits$close_to_site <- abs(reads$pos[hits$index2] - sites$pos[hits$index1]) <= site_pad
        
        result <- dplyr::tibble(
                site=sites$site[hits$index1],
                close_to_site=hits$close_to_site,
                reads[hits$index2,]) |>
            dplyr::arrange(site, tail)
        # Unfortunately we can't sort the whole file this way.
        
        yield(result)
    })
}

#' Count tail lengths
#'
#' Count tail lengths observed at each site, and how many are considered censored.
#'
#' If umi column is present, weight each read by 1/n within UMI groups.
#'
#' UMI counts are n_event, n_died.
#'
#' Read counts are n_read_event, n_read_died.
#'
#' @param min_tail Tail lengths shorter than this are discarded.
#'
#' @param length_trim This many bases are effectively trimmed from the ends of reads.
#'
#' @param must_be_close_to_site Are tail lengths only valid if reads ended very close to the site position. This applies to tail lengths from read 1, but not tail lengths from read 2.
#'
#' @export
count_tails <- function(sited_reads, min_tail, length_trim, must_be_close_to_site) {
    sited_reads <- sited_reads |>
        arrow::as_arrow_table()
        
    if (must_be_close_to_site)
        sited_reads <- sited_reads |>
            dplyr::filter(close_to_site)    # Only count tail if read ends very near actual site
    
    if ("umi" %in% names(sited_reads)) {
        sited_reads <- dplyr::mutate(
            sited_reads,
            weight = 1 / dplyr::n(),
            .by = c(site, umi))
    } else {
        sited_reads <- dplyr::mutate(
            sited_reads, 
            weight = 1)
    }
    
    cumulation <- sited_reads |>
        dplyr::transmute(
            site = site,
            weight = weight,
            length = length-.env$length_trim,
            died = tail_start+tail-1 < length,
            tail = ifelse(died, tail, length-tail_start)) |>
        dplyr::filter(tail >= .env$min_tail) |>
        dplyr::summarise(
            n_event = sum(weight),
            n_died = sum(weight * ifelse(died,1,0)),
            n_read_event = dplyr::n(),
            n_read_died = sum(died),
            .by = c(site, tail)) |>
        dplyr::arrange(site, tail) |>
        arrow::as_arrow_table()
    
    cumulation
}


#' Count UMIs (or reads if UMIs not present)
#'
#' Reads are counted even if they don't have a tail. Use for examining site expression levels.
#'
#' @export
count_umis <- function(sited_reads) {
    sited_reads <- sited_reads |>
        arrow::as_arrow_table()
    
    if ("umi" %in% names(sited_reads)) {
        sited_reads <- dplyr::mutate(
            sited_reads,
            weight = 1 / dplyr::n(),
            .by = c(site, umi))
    } else {
        sited_reads <- dplyr::mutate(
            sited_reads, 
            weight = 1)
    }
    
    cumulation <- sited_reads |>
        dplyr::transmute(
            site = site,
            weight = weight,
            multimapper = num_hits > 1) |>
        dplyr::summarise(
            n = sum(weight),
            n_read = dplyr::n(),
            n_read_multimapper = sum(multimapper),
            .by = c(site)) |>
        dplyr::arrange(site) |>
        arrow::as_arrow_table()
    
    cumulation
}


#' @export
combine_tail_counts <- function(tail_counts_list) {
    tail_counts_list <- ensure_list(tail_counts_list)
    
    tail_counts_list |>
        purrr::map(dplyr::collect) |>
        dplyr::bind_rows() |>
        dplyr::summarise(
            n_event=sum(n_event), 
            n_died=sum(n_died), 
            n_read_event=sum(n_read_event),
            n_read_died=sum(n_read_died),
            .by=c(site, tail))
}

#' @export
combine_counts <- function(counts_list) {
    counts_list <- ensure_list(counts_list)
    
    counts_list |>
        purrr::map(dplyr::collect) |>
        dplyr::bind_rows() |>
        dplyr::summarise(
            n=sum(as.numeric(n)),
            n_read=sum(as.numeric(n_read)),
            n_read_multimapper=sum(as.numeric(n_read_multimapper)),
            .by=site)
}

#' @export
calc_site_stats <- function(tq) {
    sites <- tq@sites
    
    tail_counts <- combine_tail_counts(tq@samples$tail_counts)
    
    result <- tail_counts |>
        tidyr::nest(.by=site, .key="tail_counts") |>
        dplyr::mutate(
            km=purrr::map(tail_counts, calc_km, .progress=TRUE),
            tail_n_read=purrr::map_dbl(tail_counts, \(df) sum(df$n_read_event)),
            tail_n=purrr::map_dbl(tail_counts, \(df) sum(df$n_event)),
            tail_n_died=purrr::map_dbl(tail_counts, \(df) sum(df$n_died)),
            tail10=purrr::map_dbl(km, km_quantile, 0.1),
            tail25=purrr::map_dbl(km, km_quantile, 0.25),
            tail50=purrr::map_dbl(km, km_quantile, 0.5),
            tail75=purrr::map_dbl(km, km_quantile, 0.75),
            tail90=purrr::map_dbl(km, km_quantile, 0.9))
    
    counts <- combine_counts(tq@samples$counts) |>
        dplyr::select(site,all_n=n,all_n_read=n_read,all_n_read_multimapper=n_read_multimapper)
    
    result <- dplyr::full_join(result, counts, by="site")
    
    if (!is.null(sites)) {
        # Delete any existing results
        sites <- dplyr::collect(sites)
        sites <- sites[ c("site",setdiff(colnames(sites),colnames(result))) ]
        result <- sites |>
            dplyr::full_join(result, by="site") |>
            dplyr::mutate(
                all_n=tidyr::replace_na(all_n,0),
                all_n_read=tidyr::replace_na(all_n_read,0),
                all_n_read_multimapper=tidyr::replace_na(all_n_read_multimapper,0),
                tail_n=tidyr::replace_na(tail_n,0), 
                tail_n_read=tidyr::replace_na(tail_n_read,0), 
                tail_n_died=tidyr::replace_na(tail_n_died,0))
    }
    
    result
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
