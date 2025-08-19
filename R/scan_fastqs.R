
ensure_samples <- function(samples) {
    samples |>
        tibble::as_tibble() |>
        dplyr::select(sample, barcode)
}

# Helper function for process_paired
scan_paired_chunks <- function(
        chunk1, chunk2, samples, max_mismatch, 
        clip_quality_char, clip_penalty,
        poly_penalty, suffix_penalty) {
    n <- length(chunk1) %/% 4
    i <- seq_len(n)
    readnames1 <- stringr::str_match(chunk1[(i-1)*4+1], "@(\\S*)")[,2]
    readnames2 <- stringr::str_match(chunk2[(i-1)*4+1], "@(\\S*)")[,2]
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
        best <- min(x, max_mismatch)
        i <- which(x <= best)
        if (length(i) != 1) NA else i
    })
    sample <- samples$sample[sample_index]
    barcode_mismatches <- purrr::imap_dbl(sample_index, \(index,i) 
        if (is.na(index)) NA else dists[i,index])
    
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


#' Initial loading and processing of Pooled PAT-Seq data into a parquet file
#'
#' @export
ingest_read_pairs <- function(
        output_filename, reads1, reads2, 
        samples, max_mismatch=1,
        clip_quality_char="I", clip_penalty=4,
        poly_penalty=1000, suffix_penalty=4,
        limit=Inf) {
    
    # Check/convert arguments
    samples <- ensure_samples(samples)
    
    # Begin
    
    cli::cli_progress_bar(format="Scanning read pairs | {scales::comma(cli::pb_current)} done | {scales::comma(cli::pb_rate_raw*60)}/m | {cli::pb_elapsed}")
    r1 <- file(reads1, open="r")
    withr::defer(close(r1))
    r2 <- file(reads2, open="r")
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
    
    chunk <- 2e5
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


#' Demultiplex into FASTQ files from parquet file
#'
#' Output read-1s into one FASTQ file per sample. Optionally clip off at poly(A) tail or low quality region. Read names have UMI and barcode appended.
#'
#' Various total read counts are provided in the output stats.parquet file. If clipping is disabled, statistics will reflect what would have been done if clipping was used.
#'
#' @param out_dir Directory to save fastq files.
#'
#' @param in_file Parquet filename, for file produced by ingest_read_pairs.
#'
#' @param clip Should reads be clipped. Don't use this if reads are intended for Tail Tools.
#'
#' @param clip_min_untemplated If clipping and a poly(A) tail plus expected suffix bases (UMI, barcode) was detected with at least this length, clip to this length. Otherwise clip at end of good quality region.
#'
#' @param clip_min_length If clipped read is shorter than this, discard it.
#'
#' @export
demux_reads <- function(
        out_dir, 
        in_file, 
        sample_names=NULL,
        clip=FALSE,
        clip_min_untemplated=12,
        clip_min_length=20,
        verbose=FALSE) {
    assertthat::assert_that(file.exists(in_file), msg="Input file doesn't exist.")
    
    if (is.null(sample_names)) {
        df <- arrow::open_dataset(in_file) |>
            dplyr::distinct(sample) |>
            dplyr::collect()
        sample_names <- sort(na.omit(df$sample))
    }
    
    ensure_dir(out_dir)
    
    # Do each sample separately.
    # A little inefficient, since it scans the parquet file for each sample.
    
    result <- parallel_map(sample_names, \(sample) {
        filename <- file.path(out_dir, paste0(sample, ".fastq.gz"))
        con <- gzfile(filename, open="w")
        withr::defer(close(con))
        
        # Need to record which reads had tails for site calling
        if (clip) {
            yield <- local_write_parquet(file.path(out_dir,paste0(sample,".parquet")))
            yield(dplyr::tibble(readname=character(0), has_end=logical(0), has_tail=logical(0), tail=numeric(0)))
        }
        
        total_reads <- 0
        total_reads_kept <- 0
        total_reads_ended <- 0
        total_reads_tailed <- 0
        
        scan_parquet(in_file, 
            columns=c(
                "readname", "sample", "barcode", "umi", "read_1_seq", "read_1_qual",
                "poly_a_length", "poly_a_suffix", "poly_a_start", "read_1_clip"), 
            callback=\(df) {
                df <- df |>
                    dplyr::filter(sample == .env$sample) |>
                    dplyr::mutate(
                        readname = paste0(readname, "_", barcode,"_", umi),
                        has_end  = poly_a_length+poly_a_suffix >= clip_min_untemplated,
                        has_tail = poly_a_length >= clip_min_untemplated,
                        clip = ifelse(has_end, poly_a_start-1, read_1_clip),
                        keep = clip >= clip_min_length)
                
                total_reads <<- total_reads + nrow(df)
                total_reads_kept <<- total_reads_kept +sum(df$keep)
                total_reads_ended <<- total_reads_tailed + sum(df$keep & df$has_end)
                total_reads_tailed <<- total_reads_tailed + sum(df$keep & df$has_tail)
                
                if (clip) {
                    df <- df |>
                        dplyr::filter(keep) |>
                        dplyr::mutate(
                            read_1_seq = substr(read_1_seq, 1, clip),
                            read_1_qual = substr(read_1_qual, 1, clip))
                    
                    yield(dplyr::select(df, readname, has_end, has_tail, tail=poly_a_length))
                }
                
                
                lines <- rbind(
                    paste0("@",df$readname),
                    df$read_1_seq,
                    rep("+",nrow(df)),
                    df$read_1_qual)
                writeLines(lines, con=con)
            })
        
        if (verbose) {
            message(sample, " done")
        }
        
        dplyr::tibble(
            sample=sample, 
            reads=total_reads, 
            reads_kept=total_reads_kept, 
            reads_ended=total_reads_ended,
            reads_tailed=total_reads_tailed)
    })
    
    result <- result |> dplyr::bind_rows()
    arrow::write_parquet(result, file.path(out_dir, "stats.parquet"))
    result
}

#Code to output files all at once. Would need to be vectorized at least to be performant.
#
#    dir.create(out_dir, showWarnings=FALSE)
#    n_samples <- length(sample_names)
#    connections <- vector("list", n_samples)
#    withr::defer(purrr::map(connections, \(con) if (!is.null(con)) close (con)))
#    for(i in seq_len(n_samples)) {
#        filename <- file.path(out_dir, paste0(sample_names[i], ".fastq.gz"))
#        connections[[i]] <- gzfile(filename, open="w")
#    }
#    
#    pq |>
#    dplyr::select(readname, sample, barcode, umi, read_1_seq, read_1_qual) |>
#    dplyr::filter(sample %in% .env$sample_names) |>
#    scan_query(\(df) {
#        matches <- match(df$sample, sample_names)
#        for(i in seq_len(nrow(df))) {
#            j <- matches[i]
#            seq <- df$read_1_seq[i]
#            qual <- df$read_1_qual[i]
#            # TODO: clipping logic
#            
#            writeLines(
#                c(
#                    paste0(">", df$readname[i], "_", df$barcode[i],"_",df$umi[i]),
#                    seq,
#                    "@",
#                    qual),
#                con=connections[[j]])
#        }
#    })
#}


#' Examine a random selection of reads
#' @export
reads_peek <- function(in_file, n=100, line_width=1000, seed=563) {
    withr::local_seed(seed)
    
    pq <- arrow::ParquetFileReader$create(in_file, mmap=FALSE)
    
    cat("
This is a random selection of read pairs. For each read, the three lines are:

1. Quality scores
2. DNA Sequence
3. tailquant base interpretation
     b = barcode
     u = umi
     ^ = tail
     - = any other sequence that passed quality clip

")
    
    # Minor niggles:
    # Samples without replacement.
    # Samples a row group, and then samples from the row group.
    for(nth in seq_len(n)) {
        df <- local({
            gr <- pq$ReadRowGroup( sample.int(pq$num_row_groups,1)-1 )
            dplyr::collect(gr[ sample.int(nrow(gr),1), ])
        })
        gc()
        
        cat(paste0(
            "\n\n\n",
            "Sample ", df$sample, ", ",
            df$barcode_mismatches, " barcode mismatches, ",
            "A*",df$poly_a_length, " T*",df$poly_t_length,"\n\n"))
        
        # Examine read 1
        
        seq <- stringr::str_split_1(df$read_1_seq,"")
        qual <- stringr::str_split_1(df$read_1_qual,"")
        ind <- seq_along(seq)
        anno <- rep("-", length(seq))
        anno[ind > df$read_1_clip] <- " "
        anno[ind >= df$poly_a_start & ind <= df$poly_a_start+df$poly_a_length-1] <- "^"
        
        cat("Read 1:\n")
        j <- 1
        while(j <= length(seq)) {
            k <- min(length(seq),j+line_width-1)
            ind <- seq(j,k)
            cat(paste0(
                paste(qual[ind],collapse=""),"\n",
                paste(seq[ind],collapse=""),"\n",
                paste(anno[ind],collapse=""),"\n\n"))
            j <- k+1
        }
        
        # Examine read 2
        
        seq <- stringr::str_split_1(df$read_2_seq,"")
        qual <- stringr::str_split_1(df$read_2_qual,"")
        ind <- seq_along(seq)
        anno <- rep("-", length(seq))
        anno[ind <= 8] <- "b"
        anno[ind >= 9 & ind <= 18] <- "u"
        anno[ind > df$read_2_clip] <- " "
        anno[ind >= df$poly_t_start & ind <= df$poly_t_start+df$poly_t_length-1] <- "^"
        
        cat("Read 2:\n")
        j <- 1
        while(j <= length(seq)) {
            k <- min(length(seq),j+line_width-1)
            ind <- seq(j,k)
            cat(paste0(
                paste(qual[ind],collapse=""),"\n",
                paste(seq[ind],collapse=""),"\n",
                paste(anno[ind],collapse=""),"\n\n"))
            j <- k+1
        }
    }
}