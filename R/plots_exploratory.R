
#
# These functions primarily support the QC & Exploration part of the shiny interface.
#
# However they may be useful to use from R.
#
# Some of these functions use a grouping, which is a data frame with two columns:
# - sample - a character vector of sample names
# - group  - a factor of groups
#
#


#
# Convert a matrix to a form ready to be plotted as an RLE plot.
#
rle_frame <- function(mat) {
    mat |> 
        reshape2::melt(varnames=c("feature","sample")) |> 
        dplyr::as_tibble() |>
        dplyr::mutate(median=median(value), .by=feature)
}

#
# Use a grouping to average columns in a matrix.
#
group_mean_columns <- function(mat, grouping) {
    if (is.null(grouping))
        return(mat)
    
    mat <- mat[, grouping$sample, drop=FALSE]
    groups <- factor(grouping$group)
    t(rowsum(t(mat), groups) / as.vector(table(groups)))
}

#
# Choose a set of interesting patterns.
#
# If diversity is NA, choose patterns with the greatest range.
#
# If diversity is given, choose patterns using a moderated sphering method.
#
choose_interesting <- function(patterns, n=50, diversity=NA) {
    if (is.na(diversity)) {
        scores <- apply(patterns,1,max) - apply(patterns,1,min)
    } else {
        patterns <- sweep(patterns,1,rowMeans(patterns),"-")
        decomp <- svd(sweep(patterns,2,colMeans(patterns),"-"))
        d <- decomp$d
        u <- sweep(decomp$u,2,sqrt(d**2/((d*diversity)**2+d[1]**2)),"*")
        scores <- rowSums(u*u)
    }
    
    order(scores, decreasing=TRUE) |> head(n)
}

#
# For use with tail_bins_plot
#
tq_tail_bins <- function(tq, breaks, genesums=FALSE, grouping=NULL) {
    n <- length(breaks)
    i <- seq_len(n-1)
    bins <- dplyr::tibble(
        start = breaks[i],
        end = breaks[i+1]-((i+1)<n),
        bases = end-start+1)
        
    bins$counts <- purrr::pmap(bins,\(start,end,...) tq_counts_tail_ended_in_range(tq, start, end))
    
    if (genesums) {
        bins$counts <- purrr::map(bins$counts, counts_genesums, tq)
    }
    
    lib_sizes <- Reduce(`+`, purrr::map(bins$counts, colSums))
    
    # Counts per million UMIs, and apply grouping
    bins$cpms <- purrr::pmap(bins,\(start,end,counts,...) { 
        result <- sweep(counts,2,1e6/lib_sizes,"*") |>
            group_mean_columns(grouping)
        colnames(result) <- paste0(colnames(result), " ", start, "-", end)
        result
    })
    
    bins$cpmbs <- purrr::pmap(bins,\(bases,cpms,...) cpms/bases)
    
    bins
}

#
# Heatmap, breaking down samples/groups into tail length bins
#
tail_bins_plot <- function(bins, n=50, moderation=10, diversity=100, samples_outer=TRUE, naming=NULL) {
    row_cpms <- Reduce(`+`, purrr::map(bins$cpms, rowSums)) / nrow(bins)
    
    patterns <- do.call(cbind, bins$cpmbs)
    
    # Each pattern is centered
    # Each pattern is scaled (with moderation)
    patterns <- sweep(patterns,1,rowMeans(patterns),"-")
    patterns <- sweep(patterns,1,sqrt(rowSums(patterns**2)+moderation**2),"/")
    
    keep <- choose_interesting(patterns, n, diversity)
    
    patterns <- patterns[keep,,drop=FALSE]
    
    if (samples_outer) {
        remap <- matrix(seq_len(ncol(patterns)), ncol=nrow(bins)) |> t() |> as.vector()
        patterns <- patterns[, remap, drop=FALSE]
    }
    
    if (!is.null(naming)) {
        rownames(patterns) <- naming[ rownames(patterns) ]
    }
    
    varistran::plot_heatmap(patterns, baseline=0, show_baseline=FALSE, show_tree=FALSE, scale_label="z-score")
}


#' Plot a distribution of tail lengths per sample (or group)
#'
#' @export
tq_plot_sample_tails <- function(tq, grouping=NULL) {
    
    #TODO: survival analysis!
    
    if (is.null(grouping)) {
        samples <- unique(df$sample)
        grouping <- dplyr::tibble(sample=samples, group=forcats::fct_inorder(samples))
    }
    
    df <- tq@samples |>
        dplyr::filter(sample %in% grouping$sample) |>
        purrr::pmap(\(sample, tail_counts, ...) {
            tail_counts |>
                dplyr::summarize(n_event=sum(n_event), n_died=sum(n_died), .by=tail) |>
                dplyr::collect() |>
                dplyr::mutate(sample = sample)
        }) |> 
        dplyr::bind_rows() |>
        dplyr::left_join(grouping, by="sample") |>
        dplyr::summarize(n_event=sum(n_event), n_died=sum(n_died), .by=c(group, tail))
    
    ggplot2::ggplot(df) + 
        ggplot2::aes(x=tail) + 
        ggplot2::facet_wrap(~ group, scales="free_y") +
        ggplot2::geom_col(ggplot2::aes(y=n_event), width=1, fill="#888888") +
        ggplot2::geom_col(ggplot2::aes(y=n_died), width=1, fill="#000000") +
        ggplot2::labs(x="Tail length", y="Count") +
        ggplot2::theme_bw()
}

