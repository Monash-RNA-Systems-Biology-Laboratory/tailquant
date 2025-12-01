#
# Organization of reference sequence and annotation, and STAR indexes
#

ensure_dna <- function(dna) {
    if (is.character(dna)) {
        dna <- Biostrings::readDNAStringSet(dna)
    }
    
    assertthat::assert_that(inherits(dna, "DNAStringSet"))
    
    dna
}

ensure_granges <- function(ranges) {
    if (is.character(ranges)) {
        ranges <- rtracklayer::import(ranges)
    }
    
    assertthat::assert_that(inherits(ranges, "GRanges"))
    
    ranges
}

ensure_reference <- function(ref) {
    if (is.character(ref)) {
        ref <- list(
            genome = file.path(ref, "genome.fasta"),
            annotation = file.path(ref, "annotation.gff3.bgz"),
            star = file.path(ref, "star"))
    }
    
    assertthat::assert_that(
        is.list(ref),
        all(c("genome","annotation","star") %in% names(ref)),
        inherits(ref$genome, "DNAStringSet") || file.exists(ref$genome),
        inherits(ref$genome, "GRanges") || file.exists(ref$annotation),
        dir.exists(ref$star))
    
    ref
}


#' @export
make_reference <- function(out_dir, genome, annotation) {
    genome <- ensure_dna(genome)
    annotation <- ensure_granges(annotation)
    
    ensure_dir(out_dir)
    Biostrings::writeXStringSet(genome, file.path(out_dir, "genome.fasta"))
    rtracklayer::export(annotation, file.path(out_dir,"annotation.gff3"), format="gff3", index=TRUE)
    
    # Needed by STAR
    rtracklayer::export(annotation, file.path(out_dir,"annotation.gff3"), format="gff3")
    
    rm(genome, annotation)
    
    command <- paste0("cd ", shQuote(out_dir), " && samtools faidx genome.fasta")
    message("Running: ", command)
    code <- system(command)
    assertthat::assert_that(code == 0, msg="Running samtools faidx failed.")
    
    command <- paste0("cd ", shQuote(out_dir), " && ",
        "STAR --runMode genomeGenerate",
        " --outFileNamePrefix star/",
        " --genomeDir star/",
        " --genomeFastaFiles genome.fasta",
        " --sjdbGTFfile annotation.gff3",
        " --sjdbGTFtagExonParentTranscript Parent",
        " --sjdbOverhang 299")
    message("Running: ", command)
    code <- system(command)
    assertthat::assert_that(code == 0, msg="STAR index generation failed.")
    
    unlink(file.path(out_dir,"annotation.gff3"))
}