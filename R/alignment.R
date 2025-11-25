
#
# Functions to run STAR aligner, calculate depth of coverage
#

run_aligner <- function(out_dir, in_dir, star_dir, sample_names, threads=parallel::detectCores()) {
    sample_filenames <- file.path(in_dir, paste0(sample_names,".fastq.gz"))
    assertthat::assert_that(all(file.exists(sample_filenames)))
    assertthat::assert_that(dir.exists(star_dir))
    
    ensure_dir(out_dir)
    
    for(i in seq_len(length(sample_names))) {
        command <- paste0(
            "STAR",
            " --runThreadN ", threads,
            " --genomeDir ", shQuote(star_dir),
            " --outFileNamePrefix ", shQuote(file.path(out_dir, paste0(sample_names[i],"."))),
            " --outSAMtype BAM SortedByCoordinate",
            " --readFilesIn ", shQuote(sample_filenames[i]),
            " --readFilesCommand zcat",
            " --outMultimapperOrder Random",
            " --runRNGseed 563",
            # No de-novo introns, annotated introns will still be used
            " --alignIntronMax 20")
        # By default multimappers will be output if there are less than 10 locations.
        message("Running: ", command)
        
        code <- system(command)
        assertthat::assert_that(code == 0, msg="Running STAR failed")
        
        # STAR gives the BAM file a verbose name. Also, index it.
        old_bam_filename <- file.path(out_dir, paste0(sample_names[i], ".Aligned.sortedByCoord.out.bam"))
        new_bam_filename <- file.path(out_dir, paste0(sample_names[i], ".bam"))
        file.rename(old_bam_filename, new_bam_filename)
        Rsamtools::indexBam(new_bam_filename)
    }
}


#' @export
# Note: keep_multimappers only affects "end" bigwigs!
make_sample_bigwigs <- function(out_prefix, bam_filename, parquet_filename, min_tail, keep_multimappers=TRUE) {
    coverage_fwd <- NULL
    coverage_rev <- NULL
    coverage_fwd_end <- NULL
    coverage_rev_end <- NULL
    coverage_fwd_unique <- NULL
    coverage_rev_unique <- NULL
    
    combine <- function(a,b) if (is.null(a)) b else a+b
    
    pq <- arrow::open_dataset(parquet_filename)
    
    scan_bam_chunks(
            bam_filename,
            Rsamtools::ScanBamParam(
                what="qname", 
                tag="NH",
                flag=Rsamtools::scanBamFlag(isSecondaryAlignment=FALSE)),
            \(alignments) {
        coverage_fwd <<- combine(coverage_fwd, 
            GenomicAlignments::coverage(alignments[ GenomicAlignments::strand(alignments) == "+" ]))
        coverage_rev <<- combine(coverage_rev, 
            GenomicAlignments::coverage(alignments[ GenomicAlignments::strand(alignments) == "-" ]))
        
        ranges <- methods::as(alignments, "GRanges") |> GenomicRanges::resize(1, fix="end")
        
        if (!keep_multimappers) {
            ranges <- ranges[ S4Vectors::mcols(ranges)$NH == 1 ]
        }
        
        ended_reads <- pq |>
            dplyr::inner_join(dplyr::tibble(readname=ranges$qname, i=seq_len(length(ranges))), by="readname") |>
            dplyr::filter(has_end==TRUE, tail >= .env$min_tail) |>
            dplyr::select(i) |>
            dplyr::collect()
        ranges <- ranges[ ended_reads$i ]
        coverage_fwd_end <<- combine(coverage_fwd_end, 
            GenomicRanges::coverage(ranges[ GenomicRanges::strand(ranges) == "+" ]))
        coverage_rev_end <<- combine(coverage_rev_end, 
            GenomicRanges::coverage(ranges[ GenomicRanges::strand(ranges) == "-" ]))
        
        unique <- alignments[ S4Vectors::mcols(alignments)$NH == 1 ]
        
        coverage_fwd_unique <<- combine(coverage_fwd_unique, 
            GenomicAlignments::coverage(unique[ GenomicAlignments::strand(unique) == "+" ]))
        coverage_rev_unique <<- combine(coverage_rev_unique, 
            GenomicAlignments::coverage(unique[ GenomicAlignments::strand(unique) == "-" ]))
    })
    
    rtracklayer::export(coverage_fwd, paste0(out_prefix,".fwd.all.bigwig"), "bigwig")
    rtracklayer::export(coverage_rev, paste0(out_prefix,".rev.all.bigwig"), "bigwig")
    rtracklayer::export(coverage_fwd_unique, paste0(out_prefix,".fwd.unique.bigwig"), "bigwig")
    rtracklayer::export(coverage_rev_unique, paste0(out_prefix,".rev.unique.bigwig"), "bigwig")
    rtracklayer::export(coverage_fwd_end, paste0(out_prefix,".fwd.end.bigwig"), "bigwig")
    rtracklayer::export(coverage_rev_end, paste0(out_prefix,".rev.end.bigwig"), "bigwig")
}


total_bigwigs <- function(out_filename, in_filenames) {
    result <- NULL
    combine <- function(a,b) if (is.null(a)) b else a+b
    for(filename in in_filenames) {
        result <- combine(result, rtracklayer::import(filename, "bigwig", as="RleList"))
    }
    
    rtracklayer::export(result, out_filename, "bigwig")
}

#' @export
# Note: keep_multimappers only affects "end" bigwigs!
make_bigwigs <- function(out_dir, in_dir_alignments, in_dir_parquets, sample_names, min_tail, keep_multimappers=TRUE) {
    assertthat::assert_that(length(intersect(sample_names, c("total","multimapping"))) == 0)
    
    ensure_dir(out_dir)
    
    parallel_walk(sample_names, \(sample) {
        make_sample_bigwigs(
            file.path(out_dir, sample),
            file.path(in_dir_alignments, paste0(sample, ".bam")),
            file.path(in_dir_parquets, paste0(sample, ".parquet")),
            min_tail=min_tail, keep_multimappers=keep_multimappers)
    })
    
    total_bigwigs(
        file.path(out_dir,"total.fwd.all.bigwig"),
        file.path(out_dir,paste0(sample_names,".fwd.all.bigwig")))
    
    total_bigwigs(
        file.path(out_dir,"total.rev.all.bigwig"),
        file.path(out_dir,paste0(sample_names,".rev.all.bigwig")))
    
    total_bigwigs(
        file.path(out_dir,"total.fwd.unique.bigwig"),
        file.path(out_dir,paste0(sample_names,".fwd.unique.bigwig")))
    
    total_bigwigs(
        file.path(out_dir,"total.rev.unique.bigwig"),
        file.path(out_dir,paste0(sample_names,".rev.unique.bigwig")))
    
    total_bigwigs(
        file.path(out_dir,"total.fwd.end.bigwig"),
        file.path(out_dir,paste0(sample_names,".fwd.end.bigwig")))
    
    total_bigwigs(
        file.path(out_dir,"total.rev.end.bigwig"),
        file.path(out_dir,paste0(sample_names,".rev.end.bigwig")))
    
    all <- rtracklayer::import(file.path(out_dir,"total.fwd.all.bigwig"), "bigwig", as="RleList") +
           rtracklayer::import(file.path(out_dir,"total.rev.all.bigwig"), "bigwig", as="RleList")
    unique <- rtracklayer::import(file.path(out_dir,"total.fwd.unique.bigwig"), "bigwig", as="RleList") +
              rtracklayer::import(file.path(out_dir,"total.rev.unique.bigwig"), "bigwig", as="RleList")
    
    pmax_rle <- function(x,y) S4Vectors::endoapply(x,BiocGenerics::pmax,y)
    
    multimapping <- (all-unique) / pmax_rle(all,1)
    rtracklayer::export(multimapping, file.path(out_dir,"multimapping.bigwig"), "bigwig")
}