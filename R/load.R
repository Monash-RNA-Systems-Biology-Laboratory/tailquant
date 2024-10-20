
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
            location=paste0(seqnames,":",start), 
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


#' Load data from a Tail Tools sample directory 
#'
#' Load data from a Tail Tools sample directory. Information about genomic length and tail length is included Multimappers are discarded. Only the 3' end position is kept.
#'
#' @export
load_tt_sample <- function(path) {
    bam_filename <- file.path(path, "alignments_filtered_sorted.bam")
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
            genomic_bases=a_start, 
            tail=a_end-a_start)
    
    param <- Rsamtools::ScanBamParam(what=c("qname"), tag=c("NH"))
    alignments <- GenomicAlignments::readGAlignments(bam_filename, param=param) |>
        as("GRanges") |>
        GenomicRanges::resize(1, fix="end") |>
        as.data.frame() |> 
        dplyr::as_tibble() |>
        dplyr::filter(NH == 1) |>
        dplyr::transmute(
            read=qname, 
            chr=as.character(seqnames), 
            pos=as.integer(start), 
            strand=strand_to_int(strand)) |>
        dplyr::left_join(clips, by="read") |>
        dplyr::select(!read) |>
        dplyr::arrange(chr,strand,pos) #Improves compressability
    
    alignments
}

#' @export
ingest_tt <- function(
        out_dir, in_dir, 
        site_pad=10,
        min_tail=19,
        length_trim=10,
        steps=c("samples","sites","reads","sited_reads","tail_counts","stats")) {
    
    assertthat::assert_that(dir.exists(in_dir), msg="Input directory doesn't exist.")
    
    if ("samples" %in% steps) {
        message("Step: samples")
        meta <- jsonlite::fromJSON(file.path(in_dir, "plotter-config.json"))
        sample_names <- meta$samples$name
        
        dplyr::tibble(sample=sample_names) |>
            save_parquet(out_dir,".","samples.parquet")
    }
    
    sample_names <- load_parquet(out_dir,".","samples.parquet") |> 
        dplyr::collect() |> 
        dplyr::pull(sample)
    
    if ("sites" %in% steps) {
        message("Step: sites")
        load_tt_sites(in_dir) |>
            save_parquet(out_dir,".","sites.parquet")
    }
    
    if ("reads" %in% steps) {
        message("Step: reads")
        purrr::walk(sample_names, .progress=TRUE, \(sample) {
            #message("Ingesting ", sample)
            file.path(in_dir,"samples",sample) |>
                load_tt_sample() |>
                save_parquet(out_dir,"reads",sample,".reads.parquet")
            NULL
        })
    }
    
    if ("sited_reads" %in% steps) {
        message("Step: sited_reads")
        sites <- load_parquet(out_dir,".","sites.parquet") |> dplyr::collect()
        
        purrr::walk(sample_names, .progress=TRUE, \(sample) {
            #message("Siting ", sample)
            load_parquet(out_dir,"reads",sample,".reads.parquet") |>
                site_reads(sites, site_pad=site_pad) |>
                save_parquet(out_dir,"sited_reads",sample,".sited_reads.parquet")
            NULL
        })
        
        rm(sites)
    }
    
    if ("tail_counts" %in% steps) {
        message("Step: tail_counts")
        #for(sample in sample_names) {
        purrr::walk(sample_names, .progress=TRUE, \(sample) {
            #message("Counting ", sample)
            load_parquet(out_dir,"sited_reads",sample,".sited_reads.parquet") |>
                count_tails(min_tail=min_tail, length_trim=length_trim) |>
                save_parquet(out_dir,"tail_counts",sample,".tail_counts.parquet")
            NULL
        })
    }
    
    if ("stats" %in% steps) {
        message("Step: stats")
        tq <- load_tq(out_dir)
        calc_site_stats(tq) |>
            save_parquet(out_dir,".","sites.parquet")
        rm(tq)
    }
}


#' @export
load_tq <- function(in_dir) {
    sites <- load_parquet(in_dir, ".", "sites.parquet")
    samples <- load_parquet(in_dir, ".", "samples.parquet") |> 
        dplyr::collect()
    
    samples$reads <- purrr::map(samples$sample, \(sample) {
        load_parquet(in_dir,"reads",sample,".reads.parquet")
    })
    
    samples$sited_reads <- purrr::map(samples$sample, \(sample) {
        load_parquet(in_dir,"sited_reads",sample,".sited_reads.parquet")
    })
    
    samples$tail_counts <- purrr::map(samples$sample, \(sample) {
        load_parquet(in_dir,"tail_counts",sample,".tail_counts.parquet")
    })
    
    list(sites=sites, samples=samples)
}

