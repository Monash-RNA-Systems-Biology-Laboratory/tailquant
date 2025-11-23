
# The focus in this package is on arrow and tidyverse data structures.
# These functions ease use of Bioconductor algorithms.

strand_to_int <- function(strand) ifelse(strand=="-",-1L,1L)
strand_to_char <- function(strand) ifelse(strand<0,"-","+")

find_overlaps <- function(chr1,start1,end1,strand1, chr2,start2,end2,strand2) {
    r1 <- GenomicRanges::GRanges(
        seqnames=chr1, 
        ranges=IRanges::IRanges(start=start1,end=end1),
        strand=strand_to_char(strand1))
    
    r2 <- GenomicRanges::GRanges(
        seqnames=chr2, 
        ranges=IRanges::IRanges(start=start2,end=end2),
        strand=strand_to_char(strand2))
    
    suppressWarnings(GenomicRanges::findOverlaps(r1, r2)) |> 
        as.data.frame() |> 
        dplyr::as_tibble() |>
        dplyr::transmute(index1=queryHits, index2=subjectHits)
}



# Sometimes easier to process GRanges as a dataframe with signed positions 

granges_to_sranges <- function(gr) {
    df <- as.data.frame(gr) |> dplyr::as_tibble()
    df$width <- NULL
    df$strand <- strand_to_int(df$strand)
    start <- df$start
    end <- df$end
    df$start <- ifelse(df$strand<0, -end, start)
    df$end <- ifelse(df$strand<0, -start, end)
    df
}

# Any negative start positions are trimmed, to keep IGV happy.
sranges_to_granges <- function(df, sort=FALSE) {
    start <- df$start
    end <- df$end
    df$start <- ifelse(df$strand<0, -end, start) |> pmax(1)
    df$end <- ifelse(df$strand<0, -start, end)
    df$strand <- strand_to_char(df$strand)
    
    if (sort) {
        df <- dplyr::arrange(df, seqnames, start)
    }
    
    GenomicRanges::GRanges(df)
}

sranges_find_overlaps <- function(r1, r2) {
    r1 <- sranges_to_granges(r1)
    r2 <- sranges_to_granges(r2)
    
    suppressWarnings(GenomicRanges::findOverlaps(r1, r2)) |> 
        as.data.frame() |> 
        dplyr::as_tibble() |>
        dplyr::transmute(index1=queryHits, index2=subjectHits)
}

sranges_disjoin <- function(r) {
    dplyr::select(r, seqnames, start, end, strand) |>
        sranges_to_granges() |>
        GenomicRanges::disjoin() |>
        granges_to_sranges()
}