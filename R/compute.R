
# Compute values from a tailquant directory

tq_get_cache <- function(tq, func, name=NA, version=NULL) {
    have_name <- !is.na(name)
    have <- FALSE
    
    if (have_name) {
        filename <- file.path(tq@dir, "cache", name)
    }
    
    if (have_name && file.exists(filename)) {
        #message("Have ", name)
        result <- qs2::qs_read(filename)
        have <- is.null(version) || identical(version, attr(result, "version"))
        if (!have) {
            warning("Version mismatch, updating ", name)
        }
    }
    
    if (!have) {
        #message("Computing ", name)
        result <- func()
        if (!is.null(version)) {
            attr(result, "version") <- version
        }
        
        if (have_name) {
            dir.create(file.path(tq@dir, "cache"), showWarnings=FALSE)
            qs2::qs_save(result, filename)
        }
    }
    
    result
}

tq_cached <- function(name, func, version=NULL) {
    function(tq) {
        tq_get_cache(tq, func=\() func(tq), name=name, version=version)
    }
}

#' @export
tq_site_ids <- function(tq) {
    tq@sites |> dplyr::select(site) |> dplyr::collect() |> dplyr::pull(site)
}

#' @export
tq_tail_range <- tq_cached("tail_range.qs2", \(tq) {
    ranges <- purrr::map_dfr(tq@samples$tail_counts, \(tail_counts) {
        tail_counts |>
            dplyr::summarize(min=min(tail), max=max(tail)) |> 
            dplyr::collect()
    })
    
    c(min(ranges$min), max(0, ranges$max))
})


counts_helper <- function(tq, colname) {
    sites <- tq_site_ids(tq)
    result <- matrix(0, nrow=length(sites), ncol=nrow(tq@samples))
    rownames(result) <- sites
    colnames(result) <- tq@samples$sample
    
    for(i in seq_len(ncol(result))) {
        df <- tq@samples$counts[[i]] |> dplyr::select(site, n=all_of(colname)) |> dplyr::collect()
        result[ match(df$site, sites), i ] <- df$n
    }
    
    result
}

#' @export
tq_counts <- tq_cached("counts.qs2",\(tq) {
    counts_helper(tq, "n")
})

#' @export
tq_counts_read <- tq_cached("counts_read.qs2",\(tq) {
    counts_helper(tq, "n_read")
})

#' @export
tq_counts_read_multimapper <- tq_cached("counts_read_multimapper.qs2",\(tq) {
    counts_helper(tq, "n_read_multimapper")
})


tail_count_helper <- function(tq, what, tail_min=-Inf, tail_max=Inf) {
    sites <- tq_site_ids(tq)
    result <- matrix(0, nrow=length(sites), ncol=nrow(tq@samples))
    rownames(result) <- sites
    colnames(result) <- tq@samples$sample
    
    for(i in seq_len(ncol(result))) {
        df <- tq@samples$tail_counts[[i]] |>
            dplyr::filter(tail >= tail_min, tail <= tail_max) |>
            dplyr::select(site, n=dplyr::all_of(what)) |>
            dplyr::summarize(n=sum(n), .by=site) |> 
            dplyr::collect()
        result[ match(df$site, sites), i ] <- df$n
    }
    
    result
}

#' @export
tq_counts_tail <- tq_cached("counts_tail.qs2",\(tq) {
    tail_count_helper(tq, "n_event")
})

#' @export
tq_counts_tail_at_least <- function(tq, tail_min=0) {
    name <- paste0("counts_tail_at_least_",tail_min,".qs2")
    tq_get_cache(tq, name=name, \() {
        tail_count_helper(tq, "n_event", tail_min=tail_min)
    })
}

#' @export
tq_counts_tail_ended <- tq_cached("counts_tail_ended.qs2",\(tq) {
    tail_count_helper(tq, "n_died")
})

#' @export
tq_counts_tail_ended_in_range <- function(tq, tail_min=0, tail_max=Inf) {
    name <- paste0("counts_tail_ended_in_",tail_min,"_",tail_max,".qs2")
    tq_get_cache(tq, name=name, \() {
        tail_count_helper(tq, "n_died", tail_min, tail_max)
    })
}


#' @export
tq_lib_sizes <- tq_cached("lib_sizes.qs2",\(tq) {
    counts <- tq_counts(tq)
    
    dplyr::tibble(sample=colnames(counts), lib_size=colSums(counts))
})


#' @export
tq_sample_stats <- tq_cached("sample_stats.qs2",version=3,\(tq) {
    n <- tq_counts(tq) |> colSums()
    n_tail <- tq_counts_tail(tq) |> colSums()
    n_tail_ended <- tq_counts_tail_ended(tq) |> colSums()
    n_read <- tq_counts_read(tq) |> colSums()
    
    # Backwards compatability
    n_read_multimapper <- rep(NA, length(n))
    try({
        n_read_multimapper <- tq_counts_read_multimapper(tq) |> colSums()
    }, silent=TRUE)
    
    mean_tail <- purrr::map_dbl(tq@samples$tail_counts, \(tail_counts)
        tail_counts |> 
            dplyr::summarize(mean_tail=sum(tail*n_event)/sum(n_event)) |> 
            dplyr::collect() |> 
            dplyr::pull(mean_tail))
    
    dplyr::tibble(
        sample=names(n),
        n=n,
        n_tail=n_tail,
        n_tail_ended=n_tail_ended,
        n_read=n_read,
        n_read_multimapper=n_read_multimapper,
        mean_tail=mean_tail,
        reads_per_umi=n_read/n,
        multimapping=n_read_multimapper/n_read)
})

#' Get matrix of proprtion of tails as long or longer than a certain length
#'
#' @export
#
# This could be made much faster!
tq_proportions <- function(tq, tail) {
    name <- paste0("proportions_",tail,".qs2")
    tq_get_cache(tq, name=name, \() {
        sites <- tq_site_ids(tq)
        prop <- matrix(NA_real_, nrow=length(sites), ncol=nrow(tq@samples))
        rownames(prop) <- sites
        colnames(prop) <- tq@samples$sample
        #se <- prop
        
        for(i in seq_len(ncol(prop))) {
            tail_counts <- tq@samples$tail_counts[[i]] |>
                dplyr::select(site, tail, n_event, n_died) |>
                dplyr::collect() |>
                tidyr::nest(.by=site, .key="tail_counts")
            idx <- match(tail_counts$site, sites)
            for(j in seq_along(idx)) {
                km <- calc_km(tail_counts$tail_counts[[j]])
                prop[ idx[j], i ] <- km_at(km, tail)
                #result <- km_at(km, tail)
                #prop[ idx[j], i ] <- result$prop
                #se[ idx[j], i ] <- result$se
            }
        }
        
        #list(prop=prop, se=se)
        prop
    })
}

#' Get matrix of quantile tail lengths
#'
#' @export
#
# This could be made much faster!
tq_quantiles <- function(tq, prop) {
    name <- paste0("quantiles_",prop,".qs2")
    tq_get_cache(tq, name=name, \() {
        sites <- tq_site_ids(tq)
        tail <- matrix(NA_real_, nrow=length(sites), ncol=nrow(tq@samples))
        rownames(tail) <- sites
        colnames(tail) <- tq@samples$sample
        #low <- tail
        #high <- tail
        
        for(i in seq_len(ncol(tail))) {
            tail_counts <- tq@samples$tail_counts[[i]] |>
                dplyr::select(site, tail, n_event, n_died) |>
                dplyr::collect() |>
                tidyr::nest(.by=site, .key="tail_counts")
            idx <- match(tail_counts$site, sites)
            for(j in seq_along(idx)) {
                km <- calc_km(tail_counts$tail_counts[[j]])
                tail[ idx[j], i ] <- km_quantile(km, prop)
                #low [ idx[j], i ] <- km_quantile_bound(km, prop, -z)
                #high[ idx[j], i ] <- km_quantile_bound(km, prop, z)
            }
        }
        
        #list(tail=tail, low=low, high=high)
        tail
    })
}
