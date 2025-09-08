
check_design <- function(tq, design) {
    assertthat::assert_that(
        !is.null(rownames(design)), 
        msg="Design matrix must have row names.")
    assertthat::assert_that(
        all(rownames(design) %in% tq@samples$sample), 
        msg="Design matrix row names must match sample names.")
}

# Which rows are estimable?
check_full_rank <- function(present, design, tol=1e-10) {
    full_rank <- function(selection) {
        if (sum(selection) < ncol(design)) return(FALSE)
        sum(abs(svd(design[selection,,drop=F])$d)>=tol) == ncol(design)
    }
    apply(present, 1, full_rank)
}

# Which rows should be kept?
# **DO NOT provide design for differential expression, as check_full_rank should not be applied here!**
get_keep <- function(counts, min_count, min_count_in, design=NULL) {
    lib_sizes <- colSums(counts)
    cutoffs <- min_count / median(lib_sizes) * lib_sizes
    present <- counts >= cutoffs
    keep <- rowSums(present) >= min_count_in
    
    if (!is.null(design)) {
        keep[keep] <- check_full_rank(present[keep,,drop=FALSE], design)
    }
    
    keep
}

tq_test_expression <- function(tq, design, contrasts, fdr=0.05, min_count=10, min_count_in=1, min_tail=NA, title="a test") {
    check_design(tq, design)
    
    if (is.na(min_tail)) {
        counts <- tq_counts(tq)
    } else {
        counts <- tq_counts_tail_at_least(tq, min_tail)
    }
    
    counts <- counts[,rownames(design),drop=FALSE]
    
    keep <- get_keep(counts, min_count, min_count_in)
    
    counts <- counts[keep,,drop=FALSE]
    dge <- edgeR::DGEList(counts) |>
        edgeR::calcNormFactors()
    voomed <- limma::voom(dge, design=design)
    
    result <- weitrix::weitrix_confects(voomed, design=design, contrasts=contrasts, fdr=fdr, full=TRUE)
    result$title <- paste0(
        "Differential expression",
        if (!is.na(min_tail)) paste0(", tail at least ", min_tail) else "",
        ": ", title)
    result$what <- "sites"
    
    sites <- tq@sites |>
        dplyr::select(site, gene_id, name, biotype, product) |>
        dplyr::collect()
    i <- match(result$table$name, sites$site)
    result$table$gene_name <- sites$name[i]
    result$table$biotype <- sites$biotype[i]
    
    result$plots <- list()
    result$plots[["Calibration vs sample"]] <- 
        compact_plot(weitrix::weitrix_calplot(voomed, design, cat=col))
    result$plots[["Calibration vs expression"]] <- 
        compact_plot(weitrix::weitrix_calplot(voomed, design, covar=mu))
    
    result
}

tq_test_quantile <- function(tq, prop, design, contrasts, fdr=0.05, min_count=10, min_count_in=1, title="a test") {
    check_design(tq, design)
    weights <- tq_counts_tail_ended(tq)[,rownames(design),drop=FALSE]
    tails <- tq_quantiles(tq, prop)[,rownames(design),drop=FALSE]
    
    keep <- get_keep(weights, min_count, min_count_in, design)
    
    weights <- weights[keep,,drop=FALSE]
    tails <- tails[keep,,drop=FALSE]
    
    # Tail may be missing if a lot tails weren't observed to end
    # Encode as missing data appropriately for weitrix
    bad <- is.na(tails)
    tails[bad] <- 0
    weights[bad] <- 0
    
    wei <- weitrix::as_weitrix(tails, weights)
    cal <- weitrix::weitrix_calibrate_all(
        wei, design, 
        ~ weitrix::well_knotted_spline(log(weight), 3) + weitrix::well_knotted_spline(mu, 3))
    
    result <- weitrix::weitrix_confects(cal, design=design, contrasts=contrasts, fdr=fdr, full=TRUE)
    result$title <- paste0("Differential tail length, quantile ",prop*100,"%: ", title)
    result$what <- "sites"
    
    sites <- tq@sites |>
        dplyr::select(site, gene_id, name, biotype, product) |>
        dplyr::collect()
    i <- match(result$table$name, sites$site)
    result$table$gene_name <- sites$name[i]
    result$table$biotype <- sites$biotype[i]
    
    result$plots <- list()
    result$plots[["Calibration vs sample"]] <- 
        compact_plot(weitrix::weitrix_calplot(cal, design, cat=col))
    result$plots[["Calibration vs prediction"]] <- 
        compact_plot(weitrix::weitrix_calplot(cal, design, covar=mu))
    
    result
}


tq_test_shift <- function(tq, design, contrasts, fdr=0.05, min_count=10, min_count_in=1, title="a test") {
    check_design(tq, design)
    
    counts <- tq_counts(tq)
    counts <- counts[,rownames(design),drop=FALSE]
    
    sites <- tq@sites |>
        dplyr::select(site, chr, pos, strand, gene_id, name, biotype, product) |>
        dplyr::collect() |>
        dplyr::filter(dplyr::n() >= 2, .by=gene_id) |>
        dplyr::arrange(gene_id, pos*strand)
    
    wei <- weitrix::counts_shift(
        counts[sites$site,,drop=FALSE], 
        dplyr::select(sites, group=gene_id, name=site))
    
    keep <- get_keep(weitrix::weitrix_weights(wei), min_count, min_count_in, design)
    wei <- wei[keep,]
    
    cal <- weitrix::weitrix_calibrate_all(wei, design)
    
    result <- weitrix::weitrix_confects(cal, design=design, contrasts=contrasts, fdr=fdr, full=TRUE)
    result$title <- paste0("End shift: ", title)
    result$what <- "genes"
    
    i <- match(result$table$name, sites$gene_id)
    result$table$gene_name <- sites$name[i]
    result$table$biotype <- sites$biotype[i]
    
    result$plots <- list()
    result$plots[["Calibration vs sample"]] <- 
        compact_plot(weitrix::weitrix_calplot(cal, design, cat=col))
    result$plots[["Calibration vs prediction"]] <- 
        compact_plot(weitrix::weitrix_calplot(cal, design, covar=mu))
    
    result
}


test_types <- list(
    "expression" = list(
        title="Site expression",
        func=tq_test_expression, 
        version=4),
    "expression13" = list(
        title="Site expression, tail of at least 13",
        func=\(tq, ...) tq_test_expression(tq, min_tail=13, ...), 
        version=4),
    "expression20" = list(
        title="Site expression, tail of at least 20",
        func=\(tq, ...) tq_test_expression(tq, min_tail=20, ...), 
        version=4),
    "expression30" = list(
        title="Site expression, tail of at least 30",
        func=\(tq, ...) tq_test_expression(tq, min_tail=30, ...), 
        version=4),
    "shift" = list(
        title="End shift",
        func=tq_test_shift,
        version=4),
    "quantile90" = list(
        title="Tail length, 90% are longer",
        func=\(tq, ...) tq_test_quantile(tq, 0.9, ...), 
        version=4),
    "quantile75" = list(
        title="Tail length, 75% are longer",
        func=\(tq, ...) tq_test_quantile(tq, 0.75, ...), 
        version=4),
    "quantile50" = list(
        title="Tail length, 50% are longer (median)",
        func=\(tq, ...) tq_test_quantile(tq, 0.5, ...), 
        version=4),
    "quantile25" = list(
        title="Tail length, 25% are longer",
        func=\(tq, ...) tq_test_quantile(tq, 0.25, ...), 
        version=4),
    "quantile10" = list(
        title="Tail length, 10% are longer",
        func=\(tq, ...) tq_test_quantile(tq, 0.1, ...), 
        version=4))


#' @export
tq_test <- function(
        tq,
        test,
        cache_key=NA,
        spec=NULL,
        ...) {
    name <- NA
    if (!is.na(cache_key)) 
        name <- paste0("test_", test, "_", cache_key, ".qs2")
    
    test_spec <- test_types[[test]]
    assertthat::assert_that(!is.null(test_spec), msg="Unknown test.")
    
    spec <- c(list(...), spec)
    version <- spec
    version$version <- test_spec$version
    
    func <- function() {
        do.call(test_spec$func, c(list(tq=tq), spec))
    }
    
    tq_get_cache(tq, func=func, name=name, version=version)
}

