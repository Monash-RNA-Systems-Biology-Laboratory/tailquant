
# QC and exploration shiny interface

rle_frame <- function(mat) {
    mat |> 
        reshape2::melt(varnames=c("feature","sample")) |> 
        dplyr::as_tibble() |>
        dplyr::mutate(median=median(value), .by=feature)
}

group_mean_columns <- function(mat, grouping) {
    if (is.null(grouping))
        return(mat)
    
    mat <- mat[, grouping$sample, drop=FALSE]
    groups <- factor(grouping$group)
    t(rowsum(t(mat), groups) / as.vector(table(groups)))
}


# This needs to go somewhere else...
#! @export
tq_tail_bins <- function(tq, breaks, genesums=FALSE, grouping=NULL) {
    n <- length(breaks)
    i <- seq_len(n-1)
    bins <- dplyr::tibble(
        start = breaks[i],
        end = breaks[i+1]-(i<n),
        bases = end-start+1)
        
    bins$counts <- purrr::pmap(bins,\(start,end,...) tq_counts_tail_ended_in_range(tq, start, end))
    
    if (genesums) {
        bins$counts <- purrr::map(bins$counts, counts_genesums, tq)
    }
    
    lib_sizes <- Reduce(`+`, purrr::map(bins$counts, colSums))
    
    # Counts per million UMIs, and apply grouping
    bins$cpms <- purrr::pmap(bins,\(start,end,counts,...) { 
        result <- sweep(counts,2,1e6/lib_sizes,"*") |>
            group_mean_columns(grouping)
        colnames(result) <- paste0(colnames(result), " ", start, "-", end)
        result
    })
    
    bins$cpmbs <- purrr::pmap(bins,\(bases,cpms,...) cpms/bases)
    
    bins
}

#! @export
tail_bins_plot <- function(bins, n=50, moderation=10, diversity=100, samples_outer=TRUE, naming=NULL) {
    row_cpms <- Reduce(`+`, purrr::map(bins$cpms, rowSums)) / nrow(bins)
    
    patterns <- do.call(cbind, bins$cpmbs)
    
    #keep <- rank(-row_cpms) <= top_n
    #patterns <- patterns[keep,,drop=FALSE]
    
    #foo <- scale(patterns)
    #ranges <- apply(foo,1,max) - apply(foo,1,min)
    #keep <- rank(-ranges) <= 50
    
    #patterns <- t(scale(t(patterns)))
    
    patterns <- sweep(patterns,1,rowMeans(patterns),"-")
    patterns <- sweep(patterns,1,sqrt(rowSums(patterns**2)+moderation**2),"/")
    
    decomp <- svd(sweep(patterns,2,colMeans(patterns),"-"))
#    u <- svd(patterns)$u
    #nd <- round(pd*(ncol(patterns)-1))
    #u <- decomp$u[, seq_len(nd), drop=FALSE]
    
    #print(decomp$d)
    d <- decomp$d
    u <- sweep(decomp$u,2,sqrt(d**2/((d*diversity)**2+d[1]**2)),"*")
    #u <- sweep(patterns,2,colMeans(patterns),"-")
    
    keep <- rank(-rowSums(u*u), ties.method="first") <= n
    
    patterns <- patterns[keep,,drop=FALSE]
    
    if (samples_outer) {
        remap <- matrix(seq_len(ncol(patterns)), ncol=nrow(bins)) |> t() |> as.vector()
        patterns <- patterns[, remap, drop=FALSE]
    }
    
    if (!is.null(naming)) {
        rownames(patterns) <- naming[ rownames(patterns) ]
    }
    
    varistran::plot_heatmap(patterns, baseline=0, show_baseline=FALSE, show_tree=FALSE, scale_label="z")
}


explore_ui <- function(tq) {
    
    configure <- shiny::tagList(
        shiny::p(),
        shiny::checkboxInput("explore_genes", "Aggregate to genes", value=FALSE),
        shiny::numericInput("explore_tail_min", "Use reads with tail length at least", value=NA, min=13, step=1),
        shiny::numericInput("explore_const", "Moderation log2(CPM+nnn)", value=1, min=0, step=0.5),
        shiny::textInput("explore_select", "Sample selection (regular expression)", value=""),
        shiny::textInput("explore_group", "Sample grouping (regular expression)", value="(.*)"),
        shiny::tableOutput("explore_grouping"))
    
    heatmap <- shiny::tagList(
        shiny::numericInput("explore_n", "Features to show", value=50, min=10, max=1000, step=10),
        plot_ui("explore_heatmap", width=1000, height=1000, margin_controls=FALSE))
    
    tail_bin <- shiny::tagList(
        shiny::textInput("tailbin_breaks", "Tail length breaks", value="13 20 30 40 50 60 70 80"),
        shiny::numericInput("tailbin_moderation", "Moderation", value=10, min=0, step=1),
        shiny::numericInput("tailbin_diversity", "Diversity", value=100, min=0, step=1),
        shiny::numericInput("tailbin_n", "Features to show", value=50, min=10, max=1000, step=10),
        shiny::checkboxInput("tailbin_samples_outer", "Tail bins within samples", value=TRUE),
        plot_ui("explore_tailbin_heatmap", width=1000, height=1000, margin_controls=FALSE))
    
    bslib::navset_underline(
        header=shiny::p(),
        bslib::nav_panel("Configure", configure),
        bslib::nav_panel("Biplot", plot_ui("explore_biplot", width=1000, height=600, margin_controls=FALSE)),
        bslib::nav_panel("RLE plot", plot_ui("explore_rle", width=1000, height=600, margin_controls=FALSE)),
        bslib::nav_panel("Sample vs median plot", plot_ui("explore_sample_median", width=1000, height=600, margin_controls=FALSE)),
        bslib::nav_panel("Heatmap", heatmap),
        bslib::nav_panel("Tail-bin heatmap", tail_bin))
}

explore_server <- function(input, output, session, tq) {
    title <- shiny::reactive({
        if (!is.na(input$explore_tail_min)) {
            paste0("Site expression from reads with tail ≥ ", input$explore_tail_min)
        } else {
            "Site expression from all reads"
        }
    })
    
    grouping <- shiny::reactive({
        samples <- tq@samples$sample
        if (input$explore_select != "") {
            keep <- stringr::str_detect(samples, input$explore_select)
            samples <- samples[keep]
        }
        
        match <- stringr::str_match(samples, input$explore_group)
        keep <- !is.na(match[,1])
        samples <- samples[keep]
        match <- match[keep,,drop=FALSE]
        
        if (ncol(match) == 1)
            group <- match[,1]
        else {
            group <- match[,2]
            for(i in seq_len(ncol(match)-2)) {
                group <- paste(group, match[,i+2])
            }
        } 
        
        dplyr::tibble(sample=samples, group=forcats::fct_inorder(group))
    })
    
    output$explore_grouping <- shiny::renderTable({ grouping() })
    
    naming <- shiny::reactive({
        if (input$explore_genes) {
            genes <- tq_genes(tq)
            result <- ifelse(is.na(genes$name), genes$gene_id, genes$name)
            names(result) <- genes$gene_id
        } else {
            df <- tq@sites |>
                dplyr::select(site, name) |>
                dplyr::collect()
            result <- ifelse(is.na(df$name), df$site, paste0(df$name, " (",df$site,")"))
            names(result) <- df$site
        }
        result
    })
    
    units <- shiny::reactive({ 
        paste0("log2(CPM+",input$explore_const,")") 
    })
    
    lcpm <- shiny::reactive({
        if (!is.na(input$explore_tail_min)) {
            counts <- tq_counts_tail_at_least(tq, input$explore_tail_min)
        } else {
            counts <- tq_counts(tq)
        }
        
        if (input$explore_genes) {
            counts <- counts_genesums(counts, tq)
        }
        
        result <- log2(edgeR::cpm(counts) + input$explore_const)
        rownames(result) <- naming()[ rownames(result) ]
        
        result <- group_mean_columns(result, grouping())
        
        result
    })
    
    breaks <- shiny::reactive({
        result <- stringr::str_split_1(input$tailbin_breaks, "[ ,]") |>
            as.integer() |>
            sort() |> 
            unique()
        assertthat::assert_that(length(result) >= 2, msg="Insufficient tail length breaks.")
        result
    })
    
    bins <- shiny::reactive({
        tq_tail_bins(tq, breaks(), genesums=input$explore_genes, grouping=grouping())
    })
    
    
    plot_server("explore_biplot", \() {
        varistran::plot_biplot(lcpm())
    })
    
    plot_server("explore_heatmap", \() {
        varistran::plot_heatmap(lcpm(), n=input$explore_n, show_tree=FALSE, scale_label="log2 CPM")
    })
    
    plot_server("explore_tailbin_heatmap", \() {
        tail_bins_plot(bins(), 
            n=input$tailbin_n,
            moderation=input$tailbin_moderation,
            diversity=input$tailbin_diversity,
            samples_outer=input$tailbin_samples_outer,
            naming=naming())
    })
    
    plot_server("explore_rle", \() {
        lcpm() |>
        rle_frame() |>
        ggplot2::ggplot() + 
            ggplot2::aes(y=value-median,x=forcats::fct_rev(sample)) + 
            ggplot2::geom_hline(yintercept=0) +
            ggplot2::geom_jitter(width=1/3, height=0, color="#000088", size=0.25, stroke=0) +
            ggplot2::geom_boxplot(coef=0, outliers=FALSE, width=1/3, color="#cc0000") +
            ggplot2::labs(x="", y=paste0(units(), " - Median"), title=title()) +
            ggplot2::coord_flip() +
            ggplot2::theme_bw()
    })
    
    plot_server("explore_sample_median", \() {
        lcpm() |>
        rle_frame() |>
        ggplot2::ggplot() + 
            ggplot2::aes(x=median,y=value) + 
            ggplot2::facet_wrap(~sample) + 
            ggplot2::geom_abline(color="#cc0000") + 
            ggplot2::geom_point(color="#000088", size=0.25, stroke=0) + 
            ggplot2::labs(x="Median", y=units(), title=title()) +
            ggplot2::coord_fixed() +
            ggplot2::theme_bw()
    })
}