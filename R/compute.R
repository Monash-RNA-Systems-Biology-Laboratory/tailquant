
# Compute values from a tailquant directory

tq_get_cache <- function(tq, name, func) {
    filename <- file.path(tq@dir, "cache", name)
    if (file.exists(filename)) {
        qs2::qs_read(filename)
    } else {
        result <- func(tq)
        dir.create(file.path(tq@dir, "cache"), showWarnings=FALSE)
        qs2::qs_save(result, filename)
        result
    }
}

tq_cached <- function(name, func) {
    function(tq) {
        tq_get_cache(tq, name, func)
    }
}

#' @export
tq_site_ids <- function(tq) {
    tq@sites |> dplyr::select(site) |> dplyr::collect() |> dplyr::pull(site)
}

#' @export
tq_tail_range <- tq_cached("tail_range.qs2", \(tq) {
    ranges <- map_dfr(tq@samples$tail_counts, \(tail_counts) {
        tail_counts |>
            dplyr::summarize(min=min(tail), max=max(tail)) |> 
            dplyr::collect()
    })
    
    c(min(ranges$min), max(0, ranges$max))
})

#' @export
tq_counts <- tq_cached("counts.qs2",\(tq) {
    sites <- tq_site_ids(tq)
    result <- matrix(0, nrow=length(sites), ncol=nrow(tq@samples))
    rownames(result) <- sites
    colnames(result) <- tq@samples$sample
    
    for(i in seq_len(ncol(result))) {
        df <- tq@samples$counts[[i]] |> dplyr::select(site, n) |> dplyr::collect()
        result[ match(df$site, sites), i ] <- df$n
    }
    
    result
})


tail_count_helper <- function(tq, what) {
    sites <- tq_site_ids(tq)
    result <- matrix(0, nrow=length(sites), ncol=nrow(tq@samples))
    rownames(result) <- sites
    colnames(result) <- tq@samples$sample
    
    for(i in seq_len(ncol(result))) {
        df <- tq@samples$tail_counts[[i]] |> 
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
tq_counts_tail_ended <- tq_cached("counts_tail.qs2",\(tq) {
    tail_count_helper(tq, "n_died")
})


#' @export
tq_lib_sizes <- tq_cached("lib_sizes.qs2",\(tq) {
    counts <- tq_counts(tq)
    
    dplyr::tibble(sample=colnames(counts), lib_size=colSums(counts))
})


#' Get matrix of proprtion of tails as long or longer than a certain length
#'
#' @export
#
# This could be made much faster!
tq_proportions <- function(tq, tail) {
    name <- paste0("proportions_",tail,".qs2")
    tq_get_cache(tq, name, \(tq) {
        sites <- tq_site_ids(tq)
        prop <- matrix(NA_real_, nrow=length(sites), ncol=nrow(tq@samples))
        rownames(prop) <- sites
        colnames(prop) <- tq@samples$sample
        se <- prop
        
        for(i in seq_len(ncol(prop))) {
            tail_counts <- tq@samples$tail_counts[[i]] |>
                dplyr::select(site, tail, n_event, n_died) |>
                dplyr::collect() |>
                tidyr::nest(.by=site, .key="tail_counts")
            idx <- match(tail_counts$site, sites)
            for(j in seq_along(idx)) {
                km <- calc_km(tail_counts$tail_counts[[j]])
                result <- km_at(km, tail)
                prop[ idx[j], i ] <- result$prop
                se[ idx[j], i ] <- result$se
            }
        }
        
        list(prop=prop, se=se)
    })
}

#' Get matrix of quantile tail lengths
#'
#' @export
#
# This could be made much faster!
tq_quantiles <- function(tq, prop, z=1) {
    name <- paste0("quantiles_",prop,"_",z,".qs2")
    tq_get_cache(tq, name, \(tq) {
        sites <- tq_site_ids(tq)
        tail <- matrix(NA_real_, nrow=length(sites), ncol=nrow(tq@samples))
        rownames(tail) <- sites
        colnames(tail) <- tq@samples$sample
        low <- tail
        high <- tail
        
        for(i in seq_len(ncol(tail))) {
            tail_counts <- tq@samples$tail_counts[[i]] |>
                dplyr::select(site, tail, n_event, n_died) |>
                dplyr::collect() |>
                tidyr::nest(.by=site, .key="tail_counts")
            idx <- match(tail_counts$site, sites)
            for(j in seq_along(idx)) {
                km <- calc_km(tail_counts$tail_counts[[j]])
                tail[ idx[j], i ] <- km_quantile(km, prop)
                low [ idx[j], i ] <- km_quantile_bound(km, prop, -z)
                high[ idx[j], i ] <- km_quantile_bound(km, prop, z)
            }
        }
        
        list(tail=tail, low=low, high=high)
    })
}
