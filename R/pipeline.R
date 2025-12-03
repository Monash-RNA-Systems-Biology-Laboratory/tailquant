
# 
# Main tailquant pipeline
#

tq_create_final_files <- function(tq) {
    sites <- tq@sites |>
        dplyr::select(site, location, relation, gene_id, name, biotype) |>
        dplyr::collect()
    
    genes <- tq_genes(tq) |>
        dplyr::select(gene_id, name, biotype) 
    
    ensure_dir(tq@dir, "csv")
    
    dplyr::bind_cols(sites, round(tq_counts(tq))) |>
        readr::write_csv(file.path(tq@dir, "csv", "site_counts.csv"), na="")
    
    gene_counts <- tq_counts(tq) |> counts_genesums(tq) |> round()
    dplyr::bind_cols(genes[match(rownames(gene_counts), genes$gene_id),], gene_counts) |>
        readr::write_csv(file.path(tq@dir, "csv", "gene_counts.csv"), na="")
}


#' Run tailquant pipeline
#'
#' @param tail_excess_required A site will only be used if the mean tail length exceeds the downstrand genomic A length by this much. This provides a filter for mispriming.
#'
#' @param a_prop When finding the genomic A length downstrand of a A, require this level of A purity.
#'
#' @export
run_tq <- function(
        out_dir, 
        in_dir, 
        sample_names,
        reference,
        extension,
        noncoding_prop=1/3,
        clip_min_untemplated=12,
        clip_min_length=20,
        max_t_read_1=20,
        min_t=0,
        min_tail=13,
        tail_source="read2",
        length_trim=0,
        site_upstrand=300,
        site_pad=10,
        site_suppression_pad=50,
        site_min_reads=50,
        tail_excess_required=5,
        a_prop=0.6,
        keep_multimappers=TRUE,
        steps=1:12) {
    
    args <- as.list(environment())
    
    reference <- ensure_reference(reference)
    
    ensure_dir(out_dir)
    
    with_log(file.path(out_dir,"log.txt"), \(log) {
        
        cat(file=log, paste0(
            "Working directory\n\n", getwd(),"\n\n",
            "Arguments\n\n"))
        for(name in names(args)) {
            cat(file=log, paste0(name, "=", paste(deparse(args[[name]]), collapse="\n    "), "\n"))
        }
        cat(file=log, "\n")
        
        # We can use more tail lengths if getting them from read 2 (the default)
        must_be_close_to_site <- tail_source != "read2"
        
        dplyr::tibble(sample=sample_names) |>
            arrow::write_parquet(file.path(out_dir,"samples.parquet"))
        
        if (1 %in% steps) {
            message("Step 1: demultiplex")
            demux_reads(
                out_dir=file.path(out_dir, "demux"),
                in_dir=in_dir, sample_names=sample_names,
                clip=TRUE, clip_min_untemplated=12, clip_min_length=20, max_t_read_1=20, min_t=0)
        }
        
        if (2 %in% steps) {
            message("Step 2: align")
            run_aligner(
                out_dir=file.path(out_dir, "bam"),
                in_dir=file.path(out_dir, "demux"),
                star_dir=reference$star,
                sample_names=sample_names)
        }
        
        if (3 %in% steps) {
            message("Step 3: generate bigwigs")
            make_bigwigs(
                out_dir=file.path(out_dir, "bigwig"),
                in_dir_alignments=file.path(out_dir, "bam"),
                in_dir_parquets=file.path(out_dir, "demux"),
                sample_names=sample_names,
                min_tail=min_tail,
                keep_multimappers=keep_multimappers) # Only affects "ends" bigwig
        }
        
        if (4 %in% steps) {
            message("Step 4: reads")
            ensure_dir(out_dir, "reads")
            parallel_walk(sample_names, \(sample) {
                bam_file <- file.path(out_dir,"bam",paste0(sample,".bam"))
                load_bam_into(
                    file.path(out_dir,"reads",paste0(sample,".reads.parquet")),
                    bam_file,
                    tail_source=tail_source, 
                    read_pairs_dir=in_dir,
                    keep_secondary=FALSE, # Working with raw STAR output, so need to filter these
                    keep_multimappers=keep_multimappers)
            })
        }
        
        if (5 %in% steps) {
            message("Step 5: call sites")
            sites <- call_sites(
                bigwig_dir=file.path(out_dir, "bigwig"),
                site_pad=site_pad,
                suppression_pad=site_suppression_pad,
                min_reads=site_min_reads)
            arrow::write_parquet(sites, file.path(out_dir, "sites_draft.parquet"))
        }
        
        if (6 %in% steps) {
            message("Step 6: site reads for site QC")
            ensure_dir(out_dir, "sited_reads_qc")
            parallel_walk(sample_names, \(sample) {
                sites <- arrow::open_dataset(file.path(out_dir, "sites_draft.parquet")) |> 
                    dplyr::collect()
                site_reads_into(
                    file.path(out_dir,"sited_reads_qc",paste0(sample,".sited_reads_qc.parquet")),
                    file.path(out_dir,"reads",paste0(sample,".reads.parquet")),
                    sites, site_pad=site_pad, site_upstrand=site_pad) # Only use reads immediately near site!
            })
        }
        
        if (7 %in% steps) {
            message("Step 7: make site assigner")
            gff <- ensure_granges(reference$annotation)
            result <- gff_to_site_assigner(gff,
                extension=extension, pad=site_pad, noncoding_prop=noncoding_prop)
            rtracklayer::export(result$assigner, file.path(out_dir, "assigner.gff3"), format="gff3", index=TRUE)
            rtracklayer::export(result$ambiguous, file.path(out_dir, "ambiguous.gff3"), format="gff3", index=TRUE)
            # Note: index=TRUE appends .bgz to the filename
        }
        
        if (8 %in% steps) {
            message("Step 8: QC sites")
            
            assigner <- rtracklayer::import(file.path(out_dir, "assigner.gff3.bgz"))
            
            sites <- arrow::open_dataset(file.path(out_dir, "sites_draft.parquet")) |> 
                dplyr::collect()
            sites <- sites_assign(sites, assigner)
            
            genome <- ensure_dna(reference$genome)
            names(genome) <- gsub("\\s.*","", names(genome))
            
            sites <- qc_sites(sites, 
                genome=genome,
                sited_read_parquets=file.path(out_dir,"sited_reads_qc",paste0(sample_names,".sited_reads_qc.parquet")),
                tail_excess_required=tail_excess_required,
                min_tail=min_tail,
                a_prop=a_prop)
            
            rm(genome)
            
            arrow::write_parquet(sites, file.path(out_dir, "sites_draft.parquet"))
            
            # Filter and rename sites to obtain final site list.
            sites_good <- sites |>
                dplyr::filter(keep) |>
                dplyr::select(!keep) |>
                sites_give_id()
            arrow::write_parquet(sites_good, file.path(out_dir, "sites.parquet"))
            
            sites |>
                dplyr::transmute(
                    seqnames=chr,
                    start=pos,
                    end=pos,
                    strand=strand_to_char(strand),
                    ID=site, depth, mean_tail, genomic_a_length, keep,
                    gene_id, name, product, biotype, has_cds,
                    relation, pos_vs_transcript_end, pos_vs_cds_end) |>
                dplyr::arrange(seqnames, start) |>
                GenomicRanges::GRanges() |>
                rtracklayer::export(file.path(out_dir, "sites_draft.gff3"), format="gff3", index=TRUE)
            
            sites_good |>
                dplyr::transmute(
                    seqnames=chr,
                    start=pos,
                    end=pos,
                    strand=strand_to_char(strand),
                    ID=site, depth, mean_tail, genomic_a_length,
                    gene_id, name, product, biotype, has_cds,
                    relation, pos_vs_transcript_end, pos_vs_cds_end) |>
                dplyr::arrange(seqnames, start) |>
                GenomicRanges::GRanges() |>
                rtracklayer::export(file.path(out_dir, "sites.gff3"), format="gff3", index=TRUE)
        }
        
        if (9 %in% steps) {
            message("Step 9: site reads")
            ensure_dir(out_dir, "sited_reads")
            parallel_walk(sample_names, \(sample) {
                sites <- arrow::open_dataset(file.path(out_dir,"sites.parquet"))
                site_reads_into(
                    file.path(out_dir,"sited_reads",paste0(sample,".sited_reads.parquet")),
                    file.path(out_dir,"reads",paste0(sample,".reads.parquet")),
                    sites, site_pad=site_pad, site_upstrand=site_upstrand)
                NULL
            })
        }
        
        if (10 %in% steps) {
            message("Step 10: tail counts")
            ensure_dir(out_dir,"tail_counts")
            parallel_walk(sample_names, \(sample) {
                load_parquet(out_dir,"sited_reads",sample,".sited_reads.parquet") |>
                    count_tails(min_tail=min_tail, length_trim=length_trim, must_be_close_to_site=must_be_close_to_site) |>
                    arrow::write_parquet(file.path(out_dir,"tail_counts",paste0(sample,".tail_counts.parquet")))
            })
        }
        
        if (11 %in% steps) {
            message("Step 11: counts")
            ensure_dir(out_dir,"counts")
            parallel_walk(sample_names, \(sample) {
                load_parquet(out_dir,"sited_reads",sample,".sited_reads.parquet") |>
                    count_umis() |>
                    arrow::write_parquet(file.path(out_dir,"counts",paste0(sample,".counts.parquet")))
            })
        }
        
        # Delete any cached files, as they may be out of date
        clean_up_files(file.path(out_dir, "cache"), "\\.qs2$")
        
        if (12 %in% steps) {
            message("Step 12: create final files and warm up cache")
            tq <- load_tq(out_dir)
            tq_create_final_files(tq)
            tq_warmup(tq)
        }
        
        message("Finished")
    })
}

