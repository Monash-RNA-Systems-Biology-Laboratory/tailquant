
#
# Functions to run STAR aligner, calculate depth of coverage
#

run_aligner <- function(out_dir, in_dir, star_dir, samples, threads=parallel::detectCores()) {
    sample_filenames <- file.path(in_dir, paste0(samples,".fastq.gz"))
    assertthat::assert_that(all(file.exists(sample_filenames)))
    assertthat::assert_that(dir.exists(star_dir))
    
    ensure_dir(out_dir)
    
    for(i in seq_len(length(samples))) {
        command <- paste0(
            "STAR",
            " --runThreadN ", threads,
            " --genomeDir ", shQuote(star_dir),
            " --outFileNamePrefix ", shQuote(file.path(out_dir, paste0(samples[i],"."))),
            " --outSAMtype BAM SortedByCoordinate",
            " --outStd SAM",
            " --readFilesIn ", shQuote(sample_filenames[i]),
            " --readFilesCommand zcat",
            " --outMultimapperOrder Random",
            " --runRNGseed 563",
            # No de novo introns, annotated introns will still be used
            " --alignIntronMax 20")
        # By default multimappers will be output if there are less than 10 locations.
        
        code <- system(command)
        assertthat::assert_that(code == 0, msg="Running STAR failed")
    }
}


#' @export
make_sample_bigwigs <- function(out_prefix, bam_filename, parquet_filename) {
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
            Rsamtools::ScanBamParam(what="qname", tag="NH"),
            \(alignments) {
        coverage_fwd <<- combine(coverage_fwd, 
            GenomicAlignments::coverage(alignments[ GenomicAlignments::strand(alignments) == "+" ]))
        coverage_rev <<- combine(coverage_rev, 
            GenomicAlignments::coverage(alignments[ GenomicAlignments::strand(alignments) == "-" ]))
        
        ranges <- methods::as(alignments, "GRanges") |> GenomicRanges::resize(1, fix="end")
        ended_reads <- pq |> 
            dplyr::filter(readname %in% ranges$qname, has_end==TRUE) |> 
            dplyr::select(readname) |> 
            dplyr::collect()
        ranges <- ranges[ ranges$qname %in% ended_reads$readname ]
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
make_bigwigs <- function(out_dir, in_dir_alignments, in_dir_parquets, samples) {
    assertthat::assert_that(length(intersect(samples, c("total","multimapping"))) == 0)
    
    ensure_dir(out_dir)
    
    parallel_walk(samples, \(sample) {
        make_sample_bigwigs(
            file.path(out_dir, sample),
            file.path(in_dir_alignments, paste0(sample, ".Aligned.sortedByCoord.out.bam")),
            file.path(in_dir_parquets, paste0(sample, ".parquet")))
    })
    
    total_bigwigs(
        file.path(out_dir,"total.fwd.all.bigwig"),
        file.path(out_dir,paste0(samples,".fwd.all.bigwig")))
    
    total_bigwigs(
        file.path(out_dir,"total.rev.all.bigwig"),
        file.path(out_dir,paste0(samples,".rev.all.bigwig")))
    
    total_bigwigs(
        file.path(out_dir,"total.fwd.unique.bigwig"),
        file.path(out_dir,paste0(samples,".fwd.unique.bigwig")))
    
    total_bigwigs(
        file.path(out_dir,"total.rev.unique.bigwig"),
        file.path(out_dir,paste0(samples,".rev.unique.bigwig")))
    
    total_bigwigs(
        file.path(out_dir,"total.fwd.end.bigwig"),
        file.path(out_dir,paste0(samples,".fwd.end.bigwig")))
    
    total_bigwigs(
        file.path(out_dir,"total.rev.end.bigwig"),
        file.path(out_dir,paste0(samples,".rev.end.bigwig")))
    
    all <- rtracklayer::import(file.path(out_dir,"total.fwd.all.bigwig"), "bigwig", as="RleList") +
           rtracklayer::import(file.path(out_dir,"total.rev.all.bigwig"), "bigwig", as="RleList")
    unique <- rtracklayer::import(file.path(out_dir,"total.fwd.unique.bigwig"), "bigwig", as="RleList") +
              rtracklayer::import(file.path(out_dir,"total.rev.unique.bigwig"), "bigwig", as="RleList")
    
    pmax_rle <- function(x,y) S4Vectors::endoapply(x,BiocGenerics::pmax,y)
    
    multimapping <- (all-unique) / pmax_rle(all,1)
    rtracklayer::export(multimapping, file.path(out_dir,"multimapping.bigwig"), "bigwig")
}