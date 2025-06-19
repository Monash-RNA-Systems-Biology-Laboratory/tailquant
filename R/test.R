
check_design <- function(tq, design) {
    assertthat::assert_that(
        !is.null(rownames(design)), 
        msg="Design matrix must have row names.")
    assertthat::assert_that(
        all(rownames(design) %in% tq@samples$sample), 
        msg="Design matrix row names must match sample names.")
}

get_keep <- function(counts, min_count, min_count_in) {
    lib_sizes <- colSums(counts)
    cutoffs <- min_count / median(lib_sizes) * lib_sizes
    colSums(t(counts) >= cutoffs) >= min_count_in
}

tq_test_expression <- function(tq, design, contrasts, fdr=0.05, min_count=10, min_count_in=1, title="a test") {
    check_design(tq, design)
    
    counts <- tq_counts(tq)
    counts <- counts[,rownames(design),drop=FALSE]
    
    keep <- get_keep(counts, min_count, min_count_in)
    
    counts <- counts[keep,,drop=FALSE]
    dge <- edgeR::DGEList(counts) |>
        edgeR::calcNormFactors()
    voomed <- limma::voom(dge, design=design)
    
    result <- weitrix::weitrix_confects(voomed, design=design, contrasts=contrasts, fdr=fdr, full=TRUE)
    result$title <- paste0("Differential expression: ", title)
    result$what <- "sites"
    result
}

tq_test_quantile <- function(tq, prop, design, contrasts, fdr=0.05, min_count=10, min_count_in=1, title="a test") {
    check_design(tq, design)
    weights <- tq_counts_tail_ended(tq)
    tails <- tq_quantiles(tq, prop)
    
    keep <- get_keep(weights, min_count, min_count_in)
    
    weights <- weights[keep,,drop=FALSE]
    tails <- tails[keep,,drop=FALSE]
    
    # Tail may be missing if a lot tails weren't observed to end
    # Encode as missing data appropriately for weitrix
    bad <- is.na(tails)
    tails[bad] <- 0
    weights[bad] <- 0
    
    wei <- weitrix::as_weitrix(tails, weights)
    cal <- weitrix::weitrix_calibrate_all(
        wei, design, ~ weitrix::well_knotted_spline(log(weight), 3) + weitrix::well_knotted_spline(mu, 3))
    
    result <- weitrix::weitrix_confects(cal, design=design, contrasts=contrasts, fdr=fdr, full=TRUE)
    result$title <- paste0("Differential tail length, quantile ",prop*100,"%: ", title)
    result$what <- "sites"
    result
}


test_types <- list(
    "expression" = list(
        title="Site expression",
        func=tq_test_expression, 
        version=1),
    "quantile90" = list(
        title="Tail length, 90% are longer",
        func=\(tq, ...) tq_test_quantile(tq, 0.9, ...), 
        version=1),
    "quantile75" = list(
        title="Tail length, 75% are longer",
        func=\(tq, ...) tq_test_quantile(tq, 0.75, ...), 
        version=1),
    "quantile50" = list(
        title="Tail length, 50% are longer (median)",
        func=\(tq, ...) tq_test_quantile(tq, 0.5, ...), 
        version=1),
    "quantile25" = list(
        title="Tail length, 25% are longer",
        func=\(tq, ...) tq_test_quantile(tq, 0.25, ...), 
        version=1),
    "quantile10" = list(
        title="Tail length, 10% are longer",
        func=\(tq, ...) tq_test_quantile(tq, 0.1, ...), 
        version=1))


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

