
# Sites are called using a triangular convolution.
#
# Sites suppress nearby sites with lower scores.
# Even suppressed sites perform suppression.
#
call_sites_inner <- function(vec, min_reads=50, site_pad=10, suppression_pad=50) {
    vec <- as.numeric(vec)
    n <- length(vec)
    
    # Uniform convolution
    padder <- rep(0, site_pad)
    
    #uniform <- rep(1,site_pad*2+1)
    triangle <- c(seq_len(site_pad), site_pad+1, rev(seq_len(site_pad))) / (site_pad+1)
    
    spread <- stats::filter(c(padder,vec,padder), triangle)[(site_pad+1):(site_pad+n)]
    
    #triangle <- c(seq_len(site_pad), site_pad+1, rev(seq_len(site_pad)))
    #tiebreaker <- stats::filter(c(padder,vec,padder), triangle)[(site_pad+1):(site_pad+n)]
    #tiebreaker <- vec
    
    # Consider positions in order of convolved depth.
    ordering <- order(spread, decreasing=TRUE)
    
    is_suppressed <- rep(FALSE, n)
    is_site <- rep(FALSE, n)
    
    for(i in ordering) {
        if (spread[i] < min_reads) break
        is_site[i] <- !is_suppressed[i]
        is_suppressed[max(1,i-suppression_pad) : min(n,i+suppression_pad)] <- TRUE
    }
    
    which(is_site)
}

call_sites <- function(
        bigwig_dir, genome, 
        min_reads=50, site_pad=10, suppression_pad=50,
        downstrand_length=12) {
    depth_fwd <- rtracklayer::import(file.path(bigwig_dir, "total.fwd.end.bigwig"), "bigwig", as="RleList")
    depth_rev <- rtracklayer::import(file.path(bigwig_dir, "total.rev.end.bigwig"), "bigwig", as="RleList")
    
    sites <- list(dplyr::tibble(chr=character(0), pos=numeric(0), strand=numeric(0)))
    
    for(name in names(depth_fwd)) {
        this_sites <- call_sites_inner(
            depth_fwd[[name]], min_reads=min_reads, site_pad=site_pad, suppression_pad=suppression_pad)
            
        sites[[length(sites)+1]] <- dplyr::tibble(chr=name, pos=this_sites, strand=1)
    }
    
    for(name in names(depth_rev)) {
        this_sites <- call_sites_inner(
            rev(depth_rev[[name]]), min_reads=min_reads, site_pad=site_pad, suppression_pad=suppression_pad)
        this_sites <- length(depth_rev[[name]])+1 - this_sites
            
        sites[[length(sites)+1]] <- dplyr::tibble(chr=name, pos=this_sites, strand=-1)
    }
    
    sites <- dplyr::bind_rows(sites)
    
    ranges <- sites |>
        dplyr::transmute(seqnames=chr, start=pos, end=pos, strand=ifelse(strand>=0,"+","-")) |>
        GenomicRanges::GRanges(seqinfo=Biostrings::seqinfo(genome))
    
    downstrands <- GenomicRanges::flank(ranges, downstrand_length, start=FALSE) |>
        GenomicRanges::trim()
    
    sites$downstrand <- BSgenome::getSeq(genome, downstrands) |> as.character()
    
    #penalty <- downstrand_prop / (1-downstrand_prop)
    
    #sites$downstrand_a_length <- purrr::map_dbl(sites$downstrand, \(seq) {
    #    scan_from(seq,"A",1,nchar(seq),penalty)[2]
    #})
    
    sites$a_count <- stringr::str_count(sites$downstrand, "A")
    
    sites
}



# For all exons overlapping the extension, set extension$end to exon$start-1
# Filter any negative-length extensions 

last <- function(vec) vec[length(vec)]

transcript_to_ranges <- function(gene_id, transcript_id, exons, cds, pad, extension, noncoding_prop) {
    if (!is.null(cds)) {
        trim_to <- last(cds$end)+1
    }  else {
        widths <- exons$end-exons$start+1
        total <- sum(widths)
        remainder <- ceiling(total*(1-noncoding_prop))
        i <- 1
        while(widths[i] < remainder) {
            remainder <- remainder - widths[i]
            i <- i+1
        }
        trim_to <- exons$start[i] + remainder - 1
    }
    
    assertthat::assert_that(last(exons$end)+1 >= trim_to) 
    
    seqnames <- exons$seqnames[1]
    strand <- exons$strand[1]
    
    exons <- dplyr::filter(exons, end+1 >= trim_to)
    exons$start <- pmax(exons$start, trim_to)
    exons$start <- exons$start - pad
    exons$end <- exons$end + pad
    
    exons$gene_id <- gene_id
    exons$transcript_id <- transcript_id
    if (extension <= 0)
        extensions <- NULL
    else
        extensions <- dplyr::tibble(
            type="extension", gene_id, transcript_id, seqnames, strand,
            start=last(exons$end)+1, end=start+extension-1)
    
    list(
        exons=exons,
        extensions=extensions)
}


gff_to_site_assigner <- function(gff_features, pad=10, extension=750, noncoding_prop=1/3) {
    # Convert GRanges to signed-position ranges, and ensure sorted by signed-start position.
    gff_features <- gff_features |>
        granges_to_sranges() |>
        dplyr::arrange(seqnames, start)
    
    # Parent is a list column. This is a quick way to de-listify it.
    children <- tidyr::unnest(gff_features, Parent)
    
    # Get exons and nest them into transcripts
    exons <- children |>
        dplyr::select(type, transcript_id=Parent, seqnames, start, end, strand) |>
        dplyr::filter(type=="exon") |>
        tidyr::nest(.by=transcript_id, .key="exons")
    
    # Transcripts are things that have exons!
    transcript_ids <- exons$transcript_id
    
    # Get CDS ranges and nest them into transcripts
    cds <- children |>
        dplyr::select(type, transcript_id=Parent, seqnames, start, end, strand) |>
        dplyr::filter(type=="CDS") |>
        tidyr::nest(.by=transcript_id, .key="cds")
    
    assertthat::assert_that(all( cds$transcript_id %in% transcript_ids ))
    
    transcripts <- children |>
        dplyr::select(transcript_id=ID, gene_id=Parent, biotype) |>
        dplyr::filter(transcript_id %in% transcript_ids) |>
        dplyr::left_join(exons, by="transcript_id") |>
        dplyr::left_join(cds, by="transcript_id")
    
    # Genes are things that have transcripts!
    genes <- gff_features |>
        dplyr::select(gene_id=ID, symbol=Name, product=description, gene_biotype=biotype) |>
        dplyr::filter(gene_id %in% transcripts$gene_id)
    
    transcripts <- dplyr::left_join(transcripts, genes, by="gene_id")
    
    # Using exons and CDS per transcript, produce ranges we will use to assign sites to transcripts.
    result <- purrr::pmap(transcripts, .progress=TRUE,
        \(gene_id,transcript_id,exons,cds,...) 
        transcript_to_ranges(gene_id, transcript_id, exons, cds, 
            pad=pad, extension=extension, noncoding_prop=noncoding_prop))
    exon_ranges <- purrr::map(result,"exons") |> dplyr::bind_rows()
    ext_ranges <- purrr::map(result,"extensions") |> dplyr::bind_rows()
    
    # Extension ranges are trimmed to stop before any overlapping exon.
    # If an exon overlaps the start of an extension, this will remove it.
    if (nrow(ext_ranges) > 0) {
        overlaps <- sranges_find_overlaps(ext_ranges, exon_ranges)
        for(i in seq_len(nrow(overlaps))) {
            i1 <- overlaps$index1[i]
            i2 <- overlaps$index2[i]
            ext_ranges$end[i1] <- min(ext_ranges$end[i1], exon_ranges$start[i2]-1)
        }
        ext_ranges <- dplyr::filter(ext_ranges, end >= start)
    }
    
    # If pad==0, we might have some empty exons. We can discard them now.
    exon_ranges <- dplyr::filter(exon_ranges, end >= start)
    
    exon_ranges$color <- "#888800"
    ext_ranges$color <- "#880088"
    result <- dplyr::bind_rows(exon_ranges, ext_ranges)
    
    # Convert back to conventional GRanges representation
    sranges_to_granges(result)
}


sites_assign <- function(sites, assigner, gff) {
    gff <- as.data.frame(gff) |>
        dplyr::mutate(strand = strand_to_int(strand))
    assigner <- as.data.frame(assigner) |>
        dplyr::mutate(strand = strand_to_int(strand))
    
    overlaps <- find_overlaps(
        sites$chr, sites$pos, sites$pos, sites$strand,
        assigner$seqnames, assigner$start, assigner$end, assigner$strand) |>
        tidyr::nest(.by=index1)
    
    sites$gene_id <- NA
    sites$gene_ids <- list(NULL)
    for(i in seq_len(nrow(overlaps))) {
        i1 <- overlaps$index1[i]
        i2 <- overlaps$data[[i]]$index2
        genes <- unique(assigner$gene_id[ i2 ])
        if (length(genes) == 1) {
            sites$gene_id[i1] <- genes
        }
        sites$gene_ids[[i1]] <- genes
    }
    
    sites
}


#' @export
sites_give_id <- function(sites) {
    ## Name unassigned sites by position
    
    sites_unassigned <- dplyr::filter(sites, is.na(gene_id)) |>
        dplyr::mutate(site = paste0(chr,":",pos,":",strand_to_char(strand)))
    
    ## Name assigned sites by gene symbol
    
    # Order by position
    sites_assigned <- dplyr::filter(sites, !is.na(gene_id)) |> 
        dplyr::arrange(chr, strand*pos)
    
    # Get gene ids and symbols
    genes <- sites_assigned |>
        dplyr::distinct(gene_id, name) |>
        dplyr::mutate(name = ifelse(is.na(name), gene_id, name))
    
    assertthat::assert_that(length(unique(genes$gene_id)) == nrow(genes))
    
    # Rename symbols if non-unique
    while(length(unique(genes$name)) < nrow(genes)) {
        genes <- genes |>
            dplyr::mutate(.by=name,
                modified = (dplyr::n() > 1),
                name = if (dplyr::n() > 1) paste0(name,"-",gene_id) else name)
        warning("Renaming genes to deduplicate: ", paste(genes$name[genes$modified],collapse=" "))
        genes$modified <- NULL
    }
    
    # Name assigned sites with gene symbol and position number
    sites_assigned <- sites_assigned |>
        dplyr::left_join(dplyr::select(genes,gene_id,site=name), by="gene_id") |>
        dplyr::mutate(.by=site,
            site = paste0(site,":",dplyr::row_number(),"of",dplyr::n()))
    
    ## Combine and check uniqueness
    sites <- dplyr::bind_rows(sites_assigned, sites_unassigned)
    assertthat::assert_that(length(unique(sites$site)) == nrow(sites))
    sites
}
