
#' @export
scan_paired_chunks <- function(
        chunk1, chunk2, samples, max_mismatch, 
        clip_quality_char, clip_penalty,
        poly_penalty, suffix_penalty) {
    n <- length(chunk1) %/% 4
    i <- seq_len(n)
    readnames1 <- stringr::str_match(chunk1[(i-1)*4+1], "@(\\S*)")[2]
    readnames2 <- stringr::str_match(chunk2[(i-1)*4+1], "@(\\S*)")[2]
    assertthat::assert_that(identical(readnames1, readnames2))
    
    seq1 <- chunk1[(i-1)*4+2]
    qual1 <- chunk1[(i-1)*4+4]
    seq2 <- chunk2[(i-1)*4+2]
    qual2 <- chunk2[(i-1)*4+4]
    barcode <- substr(seq2, 1, 8)
    umi <- substr(seq2, 9, 18)
    
    # Identify sample
    dists <- stringdist::stringdistmatrix(barcode, samples$barcode, method="hamming")
    sample_index <- apply(dists, 1, \(x) {
        i <- which(x <= max_mismatch)
        if (length(i) != 1) NA else i
    })
    sample <- samples$sample[sample_index]
    barcode_mismatches <- purrr::imap_dbl(sample_index, \(index,i) dists[i,index])
    
    # Identify poly(A) and poly(T) sequence
    suffix <- paste0(barcode,umi) |>
        Biostrings::DNAStringSet() |>
        Biostrings::reverseComplement() |>
        as.character()
        
    
    n1 <- purrr::map2_dbl(seq1, qual1, 
        \(seq,qual) quality_clip(seq, qual, clip_quality_char, clip_penalty))
    sr1 <- purrr::imap(seq1, 
        \(seq,i) scan_suffix(seq, "A", suffix[i], n1[i], poly_penalty, suffix_penalty))
    
    n2 <- purrr::map2_dbl(seq2, qual2, 
        \(seq,qual) quality_clip(seq,qual, clip_quality_char, clip_penalty))
    sr2 <- purrr::map2(seq2, n2, 
        \(seq,n) scan_from(seq, "T", 19, n, poly_penalty))
    
    tibble::tibble(
        readname=readnames1, 
        sample=sample,
        barcode=barcode, 
        barcode_mismatches=barcode_mismatches,
        umi=umi,
        read_1_seq=seq1,
        read_1_qual=qual1,
        read_1_length=purrr::map_dbl(seq1, nchar),
        read_1_clip=n1,
        poly_a_start=purrr::map_dbl(sr1, 1), 
        poly_a_length=purrr::map_dbl(sr1,\(sr) sr[2]-sr[1]+1),
        poly_a_suffix=purrr::map_dbl(sr1, 3),
        read_2_seq=seq2, #Expensive if not needed!
        read_2_qual=qual2, #Expensive if not needed!
        read_2_length=purrr::map_dbl(seq2, nchar),
        read_2_clip=n2, 
        poly_t_start=purrr::map_dbl(sr2, 1), 
        poly_t_length=purrr::map_dbl(sr2,\(sr) sr[2]-sr[1]+1))
}

#' @export
process_paired <- function(
        output_filename, reads1, reads2, 
        samples, max_mismatch=1,
        clip_quality_char="I", clip_penalty=4,
        poly_penalty=1000, suffix_penalty=4,
        limit=Inf) {
    
    # Check/convert arguments
    
    samples <- samples |>
        tibble::as_tibble() |>
        dplyr::select(sample, barcode)
    
    # Begin
    
    cli::cli_progress_bar("Scanning read pairs")
    r1 <- gzfile(reads1, open="r")
    #r1 <- pipe(paste0("unpigz - <", shQuote(reads1)), open="r")
    withr::defer(close(r1))
    r2 <- gzfile(reads2, open="r")
    #r2 <- pipe(paste0("unpigz - <", shQuote(reads2)), open="r")
    withr::defer(close(r2))
    
    yield <- local_write_parquet(output_filename)
    queue <- local_queue()
    
    withr::local_options(future.globals.maxSize=Inf)
    
    # Configure schema
    yield(tibble::tibble(
        readname=character(0),
        sample=character(0),
        barcode=character(0),
        barcode_mismatches=numeric(0), 
        umi=character(0),
        read_1_seq=character(0),
        read_1_qual=character(0),
        read_1_length=numeric(0),
        read_1_clip=numeric(0),
        poly_a_start=numeric(0),
        poly_a_length=numeric(0),
        poly_a_suffix=numeric(0),
        read_2_seq=character(0),
        read_2_qual=character(0),
        read_2_length=numeric(0),
        read_2_clip=numeric(0),
        poly_t_start=numeric(0),
        poly_t_length=numeric(0)))
    
    chunk <- 1e5
    total <- 0
    repeat {
        lines_to_read <- 4*max(0, min(chunk, limit-total))
        chunk1 <- readLines(r1, n=lines_to_read)
        chunk2 <- readLines(r2, n=lines_to_read)
        stopifnot(length(chunk1) == length(chunk2))
        n <- length(chunk1) %/% 4
        total <- total+n
        if (n == 0) break
        
        cli::cli_progress_update(n)
        
        future_result <- future::future(
            scan_paired_chunks(
                chunk1, chunk2, 
                samples=samples, max_mismatch=max_mismatch,
                clip_quality_char=clip_quality_char, clip_penalty=clip_penalty,
                poly_penalty=poly_penalty, suffix_penalty=suffix_penalty), 
            seed=NULL)
        queue(\(item) yield(future::value(item)), future_result)
    }
}
