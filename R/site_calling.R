
# Sites are called using a triangular convolution.
#
# Sites suppress nearby sites with lower scores.
# Even suppressed sites perform suppression.
#
call_sites_inner <- function(vec, min_reads=50, site_pad=10, suppression_pad=50) {
    n <- length(vec)
    
    padder <- rep(0, site_pad)
    
    #uniform <- rep(1,site_pad*2+1)
    triangle <- c(seq_len(site_pad), site_pad+1, rev(seq_len(site_pad))) / (site_pad+1)
    
    # Convolution. This is where most of the time in this function is spent. Lots of pointlessly convolving zeros!
    spread <- stats::filter(c(padder,as.numeric(vec),padder), triangle)[(site_pad+1):(site_pad+n)]
    
    # Consider positions exceeding depth threshold, in order of depth.
    ordering <- which(spread >= min_reads)
    ordering <- ordering[order(spread[ordering], decreasing=TRUE)]
    
    is_suppressed <- rep(FALSE, n)
    is_site <- rep(FALSE, n)
    
    for(i in ordering) {
        is_site[i] <- !is_suppressed[i]
        is_suppressed[max(1,i-suppression_pad) : min(n,i+suppression_pad)] <- TRUE
    }
    
    dplyr::tibble(
        pos = which(is_site),
        depth = spread[pos] )
}

call_sites <- function(
        bigwig_dir, genome, 
        min_reads=10, site_pad=10, suppression_pad=50) {
    depth_fwd <- rtracklayer::import(file.path(bigwig_dir, "total.fwd.end.bigwig"), "bigwig", as="RleList")
    depth_rev <- rtracklayer::import(file.path(bigwig_dir, "total.rev.end.bigwig"), "bigwig", as="RleList")
    
    sites <- list(dplyr::tibble(chr=character(0), pos=numeric(0), strand=numeric(0)))
    
    for(name in names(depth_fwd)) {
        #message(name, " +")
        this_sites <- call_sites_inner(
            depth_fwd[[name]], min_reads=min_reads, site_pad=site_pad, suppression_pad=suppression_pad)
            
        sites[[length(sites)+1]] <- dplyr::tibble(chr=name, pos=this_sites$pos, strand=1, depth=this_sites$depth)
    }
    
    for(name in names(depth_rev)) {
        #message(name, " -")
        this_sites <- call_sites_inner(
            rev(depth_rev[[name]]), min_reads=min_reads, site_pad=site_pad, suppression_pad=suppression_pad)
        this_sites$pos <- length(depth_rev[[name]])+1 - this_sites$pos
            
        sites[[length(sites)+1]] <- dplyr::tibble(chr=name, pos=this_sites$pos, strand=-1, depth=this_sites$depth)
    }
    
    sites <- dplyr::bind_rows(sites)
    
    sites$location <- paste0(sites$chr,":",sites$pos," ",strand_to_char(sites$strand))
    
    # Assign a temporary name
    sites$site <- paste0("site", seq_len(nrow(sites)))
    
    sites
}

qc_sites <- function(sites, genome, sited_read_parquets, tail_excess_required=5, min_tail=13, downstrand_length=200, a_prop=0.6) {
    # Collect observed tail length
    tail_means <- arrow::open_dataset(sited_read_parquets) |>
        dplyr::filter(tail >= .env$min_tail) |>
        dplyr::summarize(mean_tail=mean(tail), .by=c(site)) |>
        dplyr::collect()
    
    sites$mean_tail <- tail_means$mean_tail[ match(sites$site, tail_means$site) ]
    
    # Collect genomic tail length
    ranges <- sites |>
        dplyr::transmute(seqnames=chr, start=pos, end=pos, strand=ifelse(strand>=0,"+","-")) |>
        GenomicRanges::GRanges(seqinfo=Biostrings::seqinfo(genome))
    
    downstrands <- suppressWarnings(GenomicRanges::flank(ranges, downstrand_length, start=FALSE)) |>
        GenomicRanges::trim()
    
    sites$downstrand <- BSgenome::getSeq(genome, downstrands) |> as.character()
    
    penalty <- a_prop / (1-a_prop)
    
    sites$genomic_a_length <- purrr::map_dbl(sites$downstrand, \(seq) {
        scan_from(seq,"A",1,nchar(seq),penalty)[2]
    })
    
    #sites$a_count <- stringr::str_count(sites$downstrand, "A")
    
    sites$keep <- sites$mean_tail >= sites$genomic_a_length + tail_excess_required
    
    sites
}


# For all exons overlapping the extension, set extension$end to exon$start-1
# Filter any negative-length extensions 

last <- function(vec) vec[length(vec)]

transcript_to_ranges <- function(gene_id, transcript_id, exons, cds, pad, extension, noncoding_prop) {
    # Sanity check
    assertthat::assert_that( length(unique(c(exons$seqnames, cds$seqnames))) == 1 )
    assertthat::assert_that( length(unique(c(exons$strand, cds$strand))) == 1 )
    
    # Fixup we needed for some dodgy yeast annotation
    start <- min(exons$start, cds$start)
    end <- max(exons$end, cds$end)
    cds_end <- if (is.null(cds)) NA else max(cds$end)
    
    if (exons$start[1] > start) {
        warning("Exon start extended to cover CDS start for ", transcript_id)
        exons$start[1] <- start
    }
    n <- nrow(exons)
    if (exons$end[n] < end) {
        warning("Exon end extended to cover CDS end for ", transcript_id)
        exons$end[n] <- end
    }
    
    # Choose a trim point, either after the CDS or a certain proportion if non-coding.
    if (!is.null(cds)) {
        trim_to <- last(cds$end)+1
        relation <- "3'UTR"
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
        relation <- "3'Noncoding"
    }
    
    assertthat::assert_that(last(exons$end)+1 >= trim_to) 
    
    seqnames <- exons$seqnames[1]
    strand <- exons$strand[1]
    
    end_exons <- dplyr::filter(exons, end+1 >= trim_to)
    end_exons$start <- pmax(end_exons$start, trim_to)
    end_exons$start <- end_exons$start - pad
    end_exons$end <- end_exons$end + pad
    end_exons$relation <- relation
    end_exons$priority <- 1
    
    if (extension <= 0)
        extension <- NULL
    else
        extension <- dplyr::tibble(
            type="extension", 
            relation="Downstrand",
            priority=1,
            seqnames, strand,
            start=last(end_exons$end)+1, 
            end=start+extension-1)
    
    exons$start <- exons$start - pad
    exons$end <- exons$end + pad
    exons$relation <- "Exon"
    exons$priority <- 3
    
    cds$start <- cds$start - pad
    cds$end <- cds$end + pad
    cds$relation <- "CDS"
    cds$priority <- 2
    
    transcript <- dplyr::tibble( # Any part of the transcript not covered by other things is an intron.
        type="transcript",
        relation="Intron",
        priority=4,
        seqnames, strand, 
        start=start-pad, 
        end=end+pad)
    
    with_attrs <- function(sr) {
        if (is.null(sr)) return(sr)
        sr$gene_id <- gene_id
        sr$transcript_id <- transcript_id
        sr$transcript_end <- end
        sr$cds_end <- cds_end
        sr
    }
    
    list(
        transcript=with_attrs(transcript),
        exons=with_attrs(exons),
        cds=with_attrs(cds),
        end_exons=with_attrs(end_exons),
        extension=with_attrs(extension))
}


# Choose the first matching column name from a data frame
find_a_column <- function(df, options, default) {
    result <- default
    for(name in options) {
        if (name %in% names(df)) {
            result <- df[[name]]
            break
        }
    }
    result
}

# Find the highest ranking option in values
first_match <- function(values, options) {
    i <- na.omit(match(values, options))
    if (length(i) == 0) return(NA_character_)
    options[min(i)]
}


gff_to_site_assigner <- function(gff_features, extension, pad=10, noncoding_prop=1/3) {
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
        dplyr::select(transcript_id=ID, gene_id=Parent) |>
        dplyr::filter(transcript_id %in% transcript_ids) |>
        dplyr::left_join(exons, by="transcript_id") |>
        dplyr::left_join(cds, by="transcript_id")
    
    gene_ids_with_cds <- unique(transcripts$gene_id[ !purrr::map_dbl(transcripts$cds,is.null) ])
    
    # Genes are things that have transcripts!
    # Try to find the information we want about each gene.
    genes_raw <- dplyr::filter(gff_features, ID %in% .env$transcripts$gene_id) 
    genes <- dplyr::tibble(gene_id = genes_raw$ID)
    genes$name <- find_a_column(genes_raw, c("Name"), NA_character_)
    genes$product <- find_a_column(genes_raw, c("description", "Product"), NA_character_)
    genes$biotype <- find_a_column(genes_raw, c("biotype", "Biotype"), NA_character_)
    genes$has_cds <- genes$gene_id %in% gene_ids_with_cds
    
    transcripts <- dplyr::left_join(transcripts, genes, by="gene_id")
    
    # Using exons and CDS per transcript, produce ranges we will use to assign sites to transcripts.
    result <- purrr::pmap(transcripts, .progress=TRUE,
        \(gene_id,transcript_id,exons,cds,...) 
        transcript_to_ranges(gene_id, transcript_id, exons, cds, 
            pad=pad, extension=extension, noncoding_prop=noncoding_prop))
    
    transcript_ranges <- purrr::map(result,"transcript") |> dplyr::bind_rows()
    exon_ranges <- purrr::map(result,"exons") |> dplyr::bind_rows()
    cds_ranges <- purrr::map(result,"cds") |> dplyr::bind_rows()
    end_exon_ranges <- purrr::map(result,"end_exons") |> dplyr::bind_rows()
    ext_ranges <- purrr::map(result,"extension") |> dplyr::bind_rows()
    
    # Extension ranges are trimmed to stop before any overlapping end_exon.
    # If an exon overlaps the start of an extension, this will remove it.
    if (nrow(ext_ranges) > 0) {
        overlaps <- sranges_find_overlaps(ext_ranges, end_exon_ranges)
        for(i in seq_len(nrow(overlaps))) {
            i1 <- overlaps$index1[i]
            i2 <- overlaps$index2[i]
            ext_ranges$end[i1] <- min(ext_ranges$end[i1], end_exon_ranges$start[i2]-1)
        }
        ext_ranges <- dplyr::filter(ext_ranges, end >= start)
    }
    
    # If pad==0, we might have some empty exons. We can discard them now.
    end_exon_ranges <- dplyr::filter(end_exon_ranges, end >= start)
    
    result <- dplyr::bind_rows(transcript_ranges, exon_ranges, cds_ranges, end_exon_ranges, ext_ranges)
    
    # Annotate with gene information
    result <- dplyr::left_join(result, genes, by="gene_id")
    
    # Remove prefix used in ENSEMBL GFFs
    result$gene_id <- gsub("^gene:", "", result$gene_id)
    
    
    # Some of the ranges in result may overlap, and may overlap between genes.
    # We will ignore such regions for assignment purposes.
    # We perform a disjoin operation so such regions can be identified.
    
    # Disjoin and summarize
    combine <- function(vec) paste(unique(sort(vec)), collapse = "/")

    dj <- sranges_disjoin(result)
    hits <- sranges_find_overlaps(dj, result)
    df <- dplyr::select(result[hits$index2,], !c(seqnames,start,end,strand))
    df$index1 <- hits$index1
    df <- df |>
        dplyr::filter( # Only keep highest priority ranges (i.e. prefer 3' end and extension)
            .by=index1,
            priority == min(priority)) |> 
        dplyr::summarize(
            .by=index1,
            good = length(unique(gene_id)) == 1,
            gene_id = combine(gene_id),
            name = combine(name),
            product = combine(product),
            biotype = combine(biotype),
            has_cds = any(has_cds),
            transcript_end = max(transcript_end),
            cds_end = if (all(is.na(cds_end))) NA else max(cds_end, na.rm=TRUE),
            relation = first_match(relation, c("3'UTR","3'Noncoding","Downstrand","CDS","Exon","Intron")))
    df <- dplyr::bind_cols(dj[df$index1,], dplyr::select(df, !index1))
    
    df$type <- "region"
    df$color <- c(
        "Downstrand"="#880088", 
        "3'UTR"="#888800",
        "3'Noncoding"="#880000",
        "CDS"="#008800",
        "Exon"="#000088",
        "Intron"="#888888"
        )[df$relation]
    
    assigner <- df |>
        dplyr::filter(good) |>
        dplyr::select(!good)
    
    ambiguous <- df |>
        dplyr::filter(!good) |>
        dplyr::select(!good)
    
    # Convert back to conventional GRanges representation
    list(
        assigner = sranges_to_granges(assigner, sort=TRUE),
        ambiguous = sranges_to_granges(ambiguous, sort=TRUE))
}


sites_assign <- function(sites, assigner) {
    assigner <- as.data.frame(assigner) |>
        dplyr::mutate(strand = strand_to_int(strand))
    
    if (is.null(assigner$product)) 
        assigner$product <- NA_character_
    if (is.null(assigner$name)) 
        assigner$name <- NA_character_
    if (is.null(assigner$biotype)) 
        assigner$biotype <- NA_character_
    if (is.null(assigner$relation)) 
        assigner$relation <- NA_character_
    if (is.null(assigner$has_cds)) 
        assigner$has_cds <- NA
    if (is.null(assigner$transcript_end))
        assigner$transcript_end <- NA_real_
    if (is.null(assigner$cds_end))
        assigner$cds_end <- NA_real_
    
    # Guaranteed 1:1 or 1:0
    overlaps <- find_overlaps(
        sites$chr, sites$pos, sites$pos, sites$strand,
        assigner$seqnames, assigner$start, assigner$end, assigner$strand)
    
    i1 <- overlaps$index1
    i2 <- overlaps$index2
    sites$gene_id <- NA_character_
    sites$gene_id[i1] <- assigner$gene_id[i2]
    sites$name <- NA_character_
    sites$name[i1] <- assigner$name[i2]
    sites$product <- NA_character_
    sites$product[i1] <- assigner$product[i2]
    sites$biotype <- NA_character_
    sites$biotype[i1] <- assigner$biotype[i2]
    sites$has_cds <- NA
    sites$has_cds[i1] <- as.logical(assigner$has_cds[i2])
    sites$relation <- NA_character_
    sites$relation[i1] <- assigner$relation[i2]
    sites$pos_vs_transcript_end <- NA_real_
    sites$pos_vs_transcript_end[i1] <- sites$pos[i1]*sites$strand[i1] - as.numeric(assigner$transcript_end[i2])
    sites$pos_vs_cds_end <- NA_real_
    sites$pos_vs_cds_end[i1] <- sites$pos[i1]*sites$strand[i1] - as.numeric(assigner$cds_end[i2])
    
    sites
}


#' @export
sites_give_id <- function(sites) {
    # Remove any existing site name
    sites$site <- NULL
    
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
        dplyr::mutate(name = ifelse(is.na(name) | name=="", gene_id, name))
    
    assertthat::assert_that(length(unique(genes$gene_id)) == nrow(genes))
    
    # Rename symbols if non-unique
    if (length(unique(genes$name)) < nrow(genes)) {
        genes <- genes |>
            dplyr::mutate(.by=name,
                modified = (dplyr::n() > 1),
                name = if (dplyr::n() > 1) paste0(name,"-",gene_id) else name)
        warning("Renaming genes to deduplicate: ", paste(sort(genes$name[genes$modified]),collapse=" "))
        genes$modified <- NULL
    }
    assertthat::assert_that(length(unique(genes$name)) == nrow(genes), msg="Failed to deduplicate gene names.")
    
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
