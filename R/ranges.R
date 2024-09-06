
# Spend as little time as possible with GRanges objects

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
