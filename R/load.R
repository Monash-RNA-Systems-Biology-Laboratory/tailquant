
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

