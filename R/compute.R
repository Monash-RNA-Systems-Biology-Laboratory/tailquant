
# Place to keep track of number of checks on background jobs
cache_state <- rlang::new_environment(list(tries=list()))

# Compute values from a tailquant directory
# - Use lock to prevent multiply simultaneous computations, writes.
# - File is moved into place once written.
tq_get_cache <- function(tq, func, name, version=NULL, blocking=TRUE, timeout=Inf) {
    dir.create(file.path(tq@dir, "cache", "lock"), showWarnings=FALSE)
    filename <- file.path(tq@dir, "cache", name)
    lockname <- file.path(tq@dir, "cache", "lock", paste0("lock_", name))
    writename <- file.path(tq@dir, "cache", "lock", paste0("write_", name))
    
    have <- FALSE
    result <- NULL
    
    get_it <- function() {
        if (file.exists(filename)) {
            result <<- qs2::qs_read(filename)
            have <<- TRUE
        }
        
        if (have) {
            have <<- identical(version, attr(result, "version"))
            if (!have) {
                warning("Version mismatch, updating ", name)
            }
        }
    }
    
    get_it()
    if (have) return(result)
    
    if (!blocking && future::nbrOfFreeWorkers() > 0) {
        # Launch worker and forget about it.
        # Conversion to promise ensures the future is cleaned up properly.
        future::future(seed=NULL, {
                tailquant:::tq_get_cache(tq=tq, func=func, name=name, version=version, timeout=0)
                NULL
            }) |> 
            promises::as.promise() |>
            promises::catch(\(e) {
                if (!inherits(e, "error_tailquant_locked"))
                    message("Background worker error: ", e$message) 
            })
        
        # Future may already have resolved, eg if using sequential plan
        get_it()
        if (have) return(result)
        
        tries <- cache_state$tries[[ lockname ]]
        tries <- if (is.null(tries)) 1 else tries + 1
        cache_state$tries[[ lockname ]] <- tries
        
        rlang::abort("Background worker launched", class="error_tailquant_running", tailquant_tries=tries)
    }
    
    # Time to acquire the write lock
    lock <- NULL
    withr::defer({ if (!is.null(lock)) filelock::unlock(lock) })
    lock <- filelock::lock(lockname, timeout=timeout)
    if (is.null(lock))
        rlang::abort("Lock acquisition timed out.", class="error_tailquant_locked")
    
    # Try to load again, as we may have waited for a computation to finish.
    get_it()
    if (have) return(result)
    
    result <- func()
    attr(result, "version") <- version
    
    qs2::qs_save(result, writename)
    file.rename(writename, filename)
    
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



#' @export
tq_genes <- function(tq) {
    unify <- function(vec) paste(unique(vec), collapse="/")
    
    tq@sites |>
        dplyr::select(site, gene_id, name, biotype, product) |>
        dplyr::collect() |>
        dplyr::summarize(
            sites=list(site), 
            name=unify(name), 
            biotype=unify(biotype), 
            product=unify(product), 
            .by=gene_id) |>
        dplyr::arrange(gene_id)
}

#' @export
counts_genesums <- function(counts, tq) {
    df <- tq@sites |> dplyr::select(site, gene_id) |> dplyr::collect()
    groups <- df$gene_id[ match(rownames(counts), df$site) ]
    good <- !is.na(groups)
    rowsum(counts[good,,drop=FALSE], groups[good])
}


#' General statistics for a group of sites
#'
#' (Used in site table of shiny app.)
#'
#' @export
tq_site_stats <- function(tq, samples=NULL, blocking=TRUE) {
    if (is.null(samples)) 
        samples <- tq@samples$sample
    
    samples <- sort(samples)
    key <- rlang::hash(samples)
    name <- paste0("site_stats_",key,".qs2")
    
    tq_get_cache(tq, name=name, version=samples, blocking=blocking, \() {
        keep <- tq@samples$sample %in% samples
        tail_counts <- combine_tail_counts(tq@samples$tail_counts[keep])
        
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
        
        # These are big, but are used in the multi-site heatmap.
        #result$tail_counts <- NULL
        
        # These are easily recreated
        result$km <- NULL
        
        counts <- combine_counts(tq@samples$counts[keep]) |>
            dplyr::select(site,all_n=n,all_n_read=n_read,all_n_read_multimapper=n_read_multimapper)
        
        result <- dplyr::full_join(result, counts, by="site")
        
        # Ensure complete
        result <- tq@sites |> 
            dplyr::select(site) |>
            dplyr::collect() |>
            dplyr::full_join(result, by="site") |>
            dplyr::mutate(
                all_n=tidyr::replace_na(all_n,0),
                all_n_read=tidyr::replace_na(all_n_read,0),
                all_n_read_multimapper=tidyr::replace_na(all_n_read_multimapper,0),
                tail_n=tidyr::replace_na(tail_n,0), 
                tail_n_read=tidyr::replace_na(tail_n_read,0), 
                tail_n_died=tidyr::replace_na(tail_n_died,0))
        
        result
    })
}


#' Ensure some computations are cached.
#'
#' Ensure some computations are cached, specifically those likely to be used by the Shiny app.
#'
#' @param tq A TailQuant object.
#'
#' @export
tq_warmup <- function(tq) {
    tq_tail_range(tq)
    tq_sample_stats(tq)
    tq_lib_sizes(tq)
    tq_site_stats(tq)
    parallel_walk(c(0.5, 0.25, 0.75, 0.1, 0.9), \(prop) {
        tq_quantiles(tq, prop)
    })
}

