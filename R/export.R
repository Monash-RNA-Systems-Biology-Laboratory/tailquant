
#' Tailquant data export as a set of parquet files
#'
#' @export
tq_export <- function(tq, out_dir, samples=NULL) {
    
    if (is.null(samples)) {
        samples <- tq@samples |>
            dplyr::select(sample) |>
            dplyr::collect()
    }
    
    sites <- tq@sites |>
        dplyr::select(site, location, chr, pos, strand, gene=gene_id, symbol=name, biotype, product) |>
        dplyr::collect() |>
        # ( Some processing specific to the yeast reference we are using, may be removed in future
        dplyr::mutate(biotype = gsub("_gene","",biotype)) |>
        dplyr::mutate(biotype = gsub("(^|/)gene($|/)","protein_coding",biotype)) |>
        # )
        dplyr::arrange(site)
    
    sample_site_tail <- purrr::pmap_dfr(tq@samples, \(sample, tail_counts, ...) {
            tail_counts |> dplyr::collect() |> dplyr::mutate(sample=sample)
        }) |>
        dplyr::select(sample, site, tail, n_tail=n_event, n_tail_ended=n_died) |>
        dplyr::arrange(match(sample,samples$sample), site, tail)
    
    sample_site <- purrr::pmap_dfr(tq@samples, \(sample, tail_counts, counts, ...) {
            cat(sample,"\n")
            tc <- tail_counts |> dplyr::collect() |>
                tidyr::nest(.by=site, .key="tail_counts") |>
                dplyr::mutate(
                    km=purrr::map(tail_counts, calc_km, .progress=TRUE),
                    n_tail=map_dbl(tail_counts, \(tc) sum(tc$n_event)),
                    n_tail_ended=map_dbl(tail_counts, \(tc) sum(tc$n_died)),
                    tail10=purrr::map_dbl(km, km_quantile, 0.1),
                    tail25=purrr::map_dbl(km, km_quantile, 0.25),
                    tail50=purrr::map_dbl(km, km_quantile, 0.5),
                    tail75=purrr::map_dbl(km, km_quantile, 0.75),
                    tail90=purrr::map_dbl(km, km_quantile, 0.9))
            
            counts |> dplyr::collect() |> 
                dplyr::full_join(tc, by="site") |>
                dplyr::mutate(sample=sample)
        }) |>
        dplyr::select(
            sample, site, 
            n, n_tail, n_tail_ended, 
            tail10, tail25, tail50, tail75, tail90, 
            n_read, n_read_multimapper) |>
        tidyr::complete(sample=samples$sample, site=sites$site) |>
        tidyr::replace_na(list(n=0, n_tail=0, n_tail_ended=0, n_read=0, n_read_multimapper=0)) |>
        dplyr::arrange(match(sample,samples$sample), site)
    
    ensure_dir(out_dir)
    arrow::write_parquet(samples, file.path(out_dir,"sample.parquet"))
    arrow::write_parquet(sites, file.path(out_dir,"site.parquet"))
    arrow::write_parquet(sample_site_tail, file.path(out_dir,"sample_site_tail.parquet"))
    arrow::write_parquet(sample_site, file.path(out_dir,"sample_site.parquet"))
}