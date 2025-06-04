
# Filename convention
fn <- function(dir, subdir, ...) file.path(dir, subdir, paste0(...))

save_parquet <- function(obj, dir, subdir, ...) {
    dir.create(dir, showWarnings=FALSE)
    dir.create(file.path(dir, subdir), showWarnings=FALSE)
    arrow::write_parquet(obj, fn(dir, subdir, ...))
}

load_parquet <- function(dir, subdir, ...) {
    arrow::open_dataset(fn(dir, subdir, ...))
}

#' Load site information from Tail Tools output
#'
#' @return
#'
#' A tibble.
#'
#' pos is the final base of the transcript in the genome (1-based coordinates).
#'
#' @export
load_tt_sites <- function(path, reference_gff=NULL, reference_fasta=NULL) {
    if (dir.exists(path)) {
        peaks_filename <- file.path(path, "peaks/relation-child.gff")
    } else if (file.exists(path)) {
        peaks_filename <- path
    } else {
        stop("Expected a Tail Tools pipeline output directory or a gff file.")
    }
    
    peaks <- rtracklayer::import(peaks_filename, "gff3")
    peaks$Parent <- purrr::map_chr(as.list(peaks$Parent), \(v) if (length(v)==1) v else NA)
    names(peaks) <- peaks$id
    
    sites <- peaks |>
        GenomicRanges::resize(1,fix="end") |>
        dplyr::as_tibble() |>
        dplyr::transmute(
            site=id, 
            location=paste0(seqnames,":",start," ",strand), 
            chr=as.character(seqnames),
            pos=as.integer(start), 
            strand=strand_to_int(strand), 
            relation=Relation, 
            gene_id=Parent, 
            name=Name, 
            biotype=Biotype, 
            product=Product)
    
    # Discard antisense relations
    is_antisense <- (sites$relation == "Antisense") |> tidyr::replace_na(FALSE)
    sites$relation[is_antisense] <- NA
    sites$gene_id [is_antisense] <- NA
    sites$name    [is_antisense] <- NA
    sites$biotype [is_antisense] <- NA
    sites$product [is_antisense] <- NA
    
    sites
}


# Helper function to load Tail Tools clipping information
load_tt_clips <- function(path) {
    clips_filename <- file.path(path, "clipped_reads.clips.gz")
    
    # Clips file is TSV, commented header line
    # name, length, a_start, a_end, a_start_ignoring_adaptor, a_end_ignoring_adaptor, adaptor_bases
    # Zero-based!
    clips <- readr::read_tsv(clips_filename, col_types="ciiiiii", skip=1, col_names=FALSE, lazy=FALSE, progress=FALSE)
    colnames(clips) <- c(
        "name", "length", 
        "a_start", "a_end", 
        "a_start_ignoring_adaptor", "a_end_ignoring_adaptor", 
        "adaptor_bases")
    
    clips <- clips |> 
        dplyr::transmute(
            read=stringr::str_split_fixed(name,"\\s",2)[,1], 
            length=length, 
            #genomic_bases=a_start,
            tail_start=a_start,
            tail=a_end-a_start)
    
    # Extract UMI if always present in name
    umis <- stringr::str_match(clips$read, ".*_.*_(.*)")[,2]
    if (!any(is.na(umis)))
        clips$umi <- umis
    
    clips
}

# Helper function to load tailquant clipping information for poly(T)
load_read2_clips <- function(read_names, read_pairs_file) {
    # Reconcile naming. Tail Tools will be dealing with names containing barcode and umi
    read_names <- unique(read_names)
    parts <- stringr::str_match(read_names, "(.*)_.*_.*")
    base_read_names <- parts[,2]
    reads_wanted <- tibble::tibble(readname=base_read_names, read=read_names)
    
    arrow::open_dataset(read_pairs_file) |>
        dplyr::select(readname, length=read_2_length, tail_start=poly_t_start, tail=poly_t_length, umi) |>
        dplyr::inner_join(reads_wanted, by="readname") |>
        dplyr::select(!readname) |>
        dplyr::collect()
}

# Helper function to load tailquant clipping information for poly(A)
load_read1_clips <- function(read_names, read_pairs_file) {
    # Reconcile naming. Tail Tools will be dealing with names containing barcode and umi
    read_names <- unique(read_names)
    parts <- stringr::str_match(read_names, "(.*)_.*_.*")
    base_read_names <- parts[,2]
    reads_wanted <- tibble::tibble(readname=base_read_names, read=read_names)
    
    arrow::open_dataset(read_pairs_file) |>
        dplyr::select(readname, length=read_1_length, tail_start=poly_a_start, tail=poly_a_length, umi) |>
        dplyr::inner_join(reads_wanted, by="readname") |>
        dplyr::select(!readname) |>
        dplyr::collect()
}

#' Load data from a Tail Tools sample directory 
#'
#' Load data from a Tail Tools sample directory. Information about genomic length and tail length is included. Multimappers are discarded. Only the 3' end position is kept.
#'
#' @export
load_tt_sample_into <- function(dest_filename, path, tail_source, read_pairs_file=NULL, limit=NA) {
    bam_filename <- file.path(path, "alignments_filtered_sorted.bam")
    
    param <- Rsamtools::ScanBamParam(what=c("qname"), tag=c("NH"))
    #bam_file <- Rsamtools::BamFile(bam_filename, yieldSize=limit)
    #alignments <- GenomicAlignments::readGAlignments(bam_file, param=param) |>
    
    yield <- local_write_parquet(dest_filename)
    yield(dplyr::tibble(
        chr=character(0),
        pos=integer(0),
        strand=integer(0),
        length=numeric(0),
        tail_start=numeric(0),
        tail=numeric(0),
        umi=character(0)
    ))
    
    scan_bam_chunks(
        bam_filename, 
        param=param, 
        limit=limit, 
        callback=\(alignments) {
            alignments <- alignments |>
                as("GRanges") |>
                GenomicRanges::resize(1, fix="end") |>
                as.data.frame() |> 
                dplyr::as_tibble() |>
                dplyr::filter(NH == 1) |>
                dplyr::transmute(
                    read=qname, 
                    chr=as.character(seqnames), 
                    pos=as.integer(start), 
                    strand=strand_to_int(strand))
            
            if (tail_source == "tt") {
                clips <- load_tt_clips(path)
            } else if (tail_source == "read2") {
                clips <- load_read2_clips(alignments$read, read_pairs_file)
            } else if (tail_source == "read1") {
                clips <- load_read1_clips(alignments$read, read_pairs_file)
            } else {
                stop("Unknown tail source.")
            }
            
            alignments <- alignments |>
                dplyr::inner_join(clips, by="read") |>
                dplyr::select(!read)
            
            yield(alignments)
    })
    
    #alignments |>
    #    dplyr::bind_rows() |>
    #    dplyr::arrange(chr,strand,pos) #Improves compressability
}


#' Ingest Tail Tools output
#'
#' @param tail_source Should be one of "tt", "read1", or "read2". If "tt", Tail Tools tail lengths are used. If "read1", poly(A) tail lengths from read 1 are used. If "read2", poly(T) tail lengths from read 2 are used. If using "read1" or "read2", read_pairs_file must be given.
#'
#' @param read_pairs_file Parquet file or parquet files produced by ingest_read_pairs().
#'
#' @export
ingest_tt <- function(
        out_dir, in_dir,
        tail_source,
        read_pairs_file=NULL,
        site_file=NULL,
        site_pad=10,
        site_upstrand=300,
        min_tail=19,
        length_trim=10,
        limit=NA, # Max alignments to read per BAM file
        steps=1:7) {
    
    assertthat::assert_that(dir.exists(in_dir), msg="Input directory doesn't exist.")
    
    dir.create(out_dir, showWarnings=FALSE)
    
    # We can use more tail lengths if getting them from read 2
    must_be_close_to_site <- tail_source != "read2"
    if (!must_be_close_to_site)
        message("Tail lengths from read 2, more read-pairs provide tail lengths.")
    
    if (1 %in% steps) {
        message("Step 1: samples")
        meta <- jsonlite::fromJSON(file.path(in_dir, "plotter-config.json"))
        sample_names <- meta$samples$name
        
        dplyr::tibble(sample=sample_names) |>
            save_parquet(out_dir,".","samples.parquet")
    }
    
    sample_names <- load_parquet(out_dir,".","samples.parquet") |> 
        dplyr::collect() |> 
        dplyr::pull(sample)
    
    if (2 %in% steps) {
        message("Step 2: sites")
        if (is.null(site_file)) 
            site_file <- in_dir
        load_tt_sites(site_file) |>
            save_parquet(out_dir,".","sites.parquet")
    }
    
    if (3 %in% steps) {
        message("Step 3: reads")
        dir.create(file.path(out_dir, "reads"), showWarnings=FALSE)
        parallel_walk(sample_names, \(sample) {
            #message("Ingesting ", sample)
            #file.path(in_dir,"samples",sample) |>
            #    load_tt_sample(tail_source=tail_source, read_pairs_file=read_pairs_file, limit=limit) |>
            #    save_parquet(out_dir,"reads",sample,".reads.parquet")
            load_tt_sample_into(
                file.path(out_dir,"reads",paste0(sample,".reads.parquet")),
                file.path(in_dir,"samples",sample),
                tail_source=tail_source, read_pairs_file=read_pairs_file, limit=limit)
            NULL
        })
    }
    
    if (4 %in% steps) {
        message("Step 4: sited_reads")
        
        parallel_walk(sample_names, \(sample) {
            #message("Siting ", sample)
            sites <- load_parquet(out_dir,".","sites.parquet") |> dplyr::collect()
            load_parquet(out_dir,"reads",sample,".reads.parquet") |>
                site_reads(sites, site_pad=site_pad, site_upstrand=site_upstrand) |>
                save_parquet(out_dir,"sited_reads",sample,".sited_reads.parquet")
            NULL
        })
    }
    
    if (5 %in% steps) {
        message("Step 5: tail_counts")
        #for(sample in sample_names) {
        parallel_walk(sample_names, \(sample) {
            #message("Tail counting ", sample)
            load_parquet(out_dir,"sited_reads",sample,".sited_reads.parquet") |>
                count_tails(min_tail=min_tail, length_trim=length_trim, must_be_close_to_site=must_be_close_to_site) |>
                save_parquet(out_dir,"tail_counts",sample,".tail_counts.parquet")
            NULL
        })
    }
    
    if (6 %in% steps) {
        message("Step 6: counts")
        parallel_walk(sample_names, \(sample) {
            #message("Counting ", sample)
            load_parquet(out_dir,"sited_reads",sample,".sited_reads.parquet") |>
                count_umis() |>
                save_parquet(out_dir,"counts",sample,".counts.parquet")
            NULL
        })
    }
    
    if (7 %in% steps) {
        message("Step 7: stats")
        tq <- load_tq(out_dir)
        calc_site_stats(tq) |>
            save_parquet(out_dir,".","sites.parquet")
        rm(tq)
    }
    
    # Delete any cached files, as they may be out of date
    cache_dir <- file.path(out_dir, "cache")
    cached_files <- list.files(cache_dir, full.names=TRUE)
    for(filename in cached_files) {
        unlink(filename)
    }
}


# Backwards compatability
fix_column_names <- function(pq) {
    if ("genomic_bases" %in% names(pq))
        pq <- dplyr::rename(pq, tail_start = genomic_bases)
    pq
}


#' @export
load_tq <- function(in_dir) {
    sites <- load_parquet(in_dir, ".", "sites.parquet")
    
    samples <- load_parquet(in_dir, ".", "samples.parquet") |> 
        dplyr::collect()
    
    samples$reads <- purrr::map(samples$sample, \(sample) {
        load_parquet(in_dir,"reads",sample,".reads.parquet") |>
            fix_column_names()
    })
    
    samples$sited_reads <- purrr::map(samples$sample, \(sample) {
        load_parquet(in_dir,"sited_reads",sample,".sited_reads.parquet") |>
            fix_column_names()
    })
    
    samples$tail_counts <- purrr::map(samples$sample, \(sample) {
        load_parquet(in_dir,"tail_counts",sample,".tail_counts.parquet")
    })
    
    samples$counts <- purrr::map(samples$sample, \(sample) {
        load_parquet(in_dir,"counts",sample,".counts.parquet")
    })
    
    new("TailQuant", dir=in_dir, sites=sites, samples=samples)
}

