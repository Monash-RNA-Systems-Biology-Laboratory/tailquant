


site_table_search <- function(df, query) {
    parts <- stringr::str_split_1(query, "\\s+")
    parts <- parts[nchar(parts) > 0]
    
    keep <- rep(TRUE, nrow(df))
    for(part in parts) {
        this_query <- stringr::fixed(part, ignore_case=TRUE)
        cols <- colnames(df)
        found <- rep(FALSE, nrow(df))
        search <- TRUE
        match <- stringr::str_match(part,"(.*?)=(.*)")
        if (!is.na(match[2])) {
            if (match[2] == "top") {
                n <- as.integer(match[3])
                if (!is.na(n))
                    found <- seq_len(nrow(df)) <= n
                search <- FALSE
            } else if (match[2] %in% cols && match[3] != "") {
                cols <- match[2]
                this_query <- stringr::regex(match[3], ignore_case=TRUE)
            }
        }
        
        if (search) {
            for(col in cols) {
                hits <- stringr::str_detect(df[[col]], this_query)
                hits[is.na(hits)] <- FALSE
                found <- found | hits
            }
        }
        keep <- keep & found
    }
    df <- df[keep,]
    df
}



site_ui <- function(tq, max_tail=NA) {
    max_tail_upper <- tq_tail_range(tq)[2]
    if (is.na(max_tail)) {
        max_tail <- max_tail_upper
    }
    
    options_panel <- bslib::nav_panel("Options",
        shiny::p(),
        shiny::div(
            style="width: 50%",
            #shiny::numericInput("top_n", "List this many top sites by n_tail_ended (0=all)", value=0, min=0, step=1),
            shiny::numericInput("heatmap_rows", "Show this many sites in heatmap", value=50, min=1, max=2000, step=1),
            shiny::numericInput("max_tail", "Maximum tail length in plots", value=max_tail, min=1, max=max_tail_upper, step=1),
            shiny::numericInput("step", "Density plot bin size", value=1, min=1, step=1),
            shiny::checkboxInput("cpm", "Plots show expression level.", value=TRUE),
            shiny::checkboxInput("show_samples", "Show individual samples in plots.", value=TRUE),
            shiny::checkboxInput("assume_all_died", "Use old method, treating tail-to-end-of-read as actual tail length."),
        )
    )
    
    samples_panel <- bslib::nav_panel("Sample selection",
        shiny::p(),
        shiny::checkboxGroupInput("which_samples", "Samples to show", 
            choices=tq@samples$sample,
            selected=tq@samples$sample),
        shiny::p("Note: This sample selection does not affect aggregated tail length statistics, just plots showing individual samples.")
    )
    
    table_panel <- bslib::nav_panel("Site selection",
        shiny::textInput("search", label="", value="", placeholder="Search, see explanation, examples: top=1000 name=^PAP", width="50%"),
        DT::DTOutput("table", fill=FALSE)
    )
    
    explanation_panel <- bslib::nav_panel("Explanation",
        shiny::p(),
        shiny::p("Searching:"),
        shiny::p(shiny::strong("top=NNN"), " for the top NNN genes by n_tail_ended."),
        shiny::p(shiny::strong("column_name=regex "), " to search particular columns."),
        shiny::p(),
        shiny::p("n_reads = Number of reads."),
        shiny::p("n_tail_reads = Number of reads ending near the actual site rather than upstrand, and with a poly(A) tail."),
        shiny::p("n = Number of UMIs."),
        shiny::p("n_tail = Number of UMIs with a poly(A) tail."),
        shiny::p("n_tail_ended = Number of UMIs with a poly(A) tail that ended before the end of the read."),
        shiny::p("tail90, tail50, tail10 = 90%/50%/10% of tails are longer than this. tail50 is the median tail length."),
        #shiny::p("tightness = tail90/tail10. Sort by this column to see sites with tight or wide tail length distributions. Best to limit the list to some number of top sites by n_tail_ended first."),
        shiny::p("cpm = Counts Per Million, calculated from n.")
    )
    
    ui <- shiny::tagList(
        shiny::wellPanel(
            bslib::navset_underline(
                table_panel,
                options_panel,
                samples_panel,
                explanation_panel
            )
        ),
        shiny::p(),
        shiny::uiOutput("site_info"),
        bslib::navset_underline(
            header=shiny::p(),
            bslib::nav_panel("Site table", DT::DTOutput("read_counts", fill=FALSE)),
            bslib::nav_panel("Site reverse cumulative distribution", plot_ui("survival_plot", width=1000, height=600)),
            bslib::nav_panel("Site density", plot_ui("density_plot", width=1000, height=600)),
            bslib::nav_panel("Site ridgeline", plot_ui("ridgeline_plot", width=1000, height=600)),
            bslib::nav_panel("Site heatmap", plot_ui("density_heatmap", width=1000, height=600)),
            bslib::nav_panel("Site read details", plot_ui("detail_plot", width=1000, height=800)),
            bslib::nav_panel("Mult-site heatmap", plot_ui("heatmap", width=1000, height=800))
        ),
        
        shiny::div(style="height: 800px;")
    )
    
    ui
}

site_server <- function(input, output, session, tq) { #, get_sites_wanted=function() NULL) {
    sites <- tq@sites
    samples <- tq@samples
    
    min_tail <- tq_tail_range(tq)[1]
    
    # Possibly this should be done in tq_load
    lib_sizes <- tq_lib_sizes(tq)
    samples <- dplyr::left_join(samples, lib_sizes, by="sample")
    
    search <- debounce(reactive(input$search), 100)
    
    df <- shiny::reactive({
        result <- sites
        # Backwards compatability
        if (!"all_n_read" %in% names(result))
            result <- dplyr::mutate(result, all_n_read=NA)
        if (!"all_n_read_multimapper" %in% names(result))
            result <- dplyr::mutate(result, all_n_read_multimapper=NA)
        if (!"all_n" %in% names(result))
            result <- dplyr::mutate(result, all_n=NA)
        
        result <- result |>
            dplyr::transmute(
                site, name, location, n=all_n, tail_n, tail_n_ended=tail_n_died, 
                tail90, tail50, tail10, 
                multimapping=all_n_read_multimapper / all_n_read,
                #tightness=tail90/tail10, 
                relation, biotype, gene_id, product) |>
            dplyr::arrange(-tail_n_ended)
        
        #want <- get_sites_wanted()
        #if (length(want) > 0) {
        #    result <- dplyr::filter(result, site %in% want)
        #} else if (input$top_n > 0) {
        #    result <- dplyr::slice_head(result, n=input$top_n)
        #}
        
        result <- dplyr::collect(result)
        
        #if (length(want) > 0) {
        #    rows <- match(result$site, want) |> na.omit()
        #    result <- result[rows,]
        #}
        
        result <- site_table_search(result, search())
        
        #if (input$top_n > 0) {
        #    result <- dplyr::slice_head(result, n=input$top_n)
        #}
        
        result
    })
    
    selected <- reactive({
        req(input$table_rows_current)
        req(nrow(df()) > 0)
        
        row <- input$table_rows_current[1]
        if (length(input$table_rows_selected) == 1 &&
            input$table_rows_selected %in% input$table_rows_current)
            row <- input$table_rows_selected
        site <- df()$site[row]
        result <- sites |>
            dplyr::filter(site == .env$site) |>
            dplyr::collect()
        
        req(nrow(result) == 1)
        result
    })
    
    selected_samples <- reactive({
        result <- dplyr::filter(samples, sample %in% .env$input$which_samples)
        
        if (!"color" %in% colnames(result)) {
            result$color <- scales::hue_pal()(nrow(result))
        }
        
        result
    })
    
    selected_kms <- reactive({
        site <- selected()$site
        
        if (input$show_samples) {
            names <- selected_samples()$sample
            tail_counts <- selected_samples()$tail_counts |>
                purrr::map(dplyr::filter, site == .env$site)
        } else {
            names <- "Combined"
            tail_counts <- dplyr::filter(sites, site == .env$site) |> 
                dplyr::collect()
            tail_counts <- tail_counts$tail_counts
        }
        
        dplyr::tibble(
            name = names,
            tail_counts = tail_counts,
            km = purrr::map(tail_counts, calc_km, assume_all_died=input$assume_all_died))
    })
    
    selected_reads <- reactive({
        site <- selected()$site
        purrr::map(samples$sited_reads, \(item) {
            dplyr::filter(item, site == .env$site)
        })
    })
    
    output$table <- DT::renderDT(server=TRUE, {
        df_show <- dplyr::select(df(), !product)
        DT::datatable(
                df_show,
                selection='single',
                rownames=FALSE, 
                #width="100%", 
                #class='compact cell-border hover',
                #extensions='Buttons'
                options=list(searching=FALSE)) |>
            DT::formatStyle(names(df_show), lineHeight="100%", paddingTop="5px", paddingBottom="5px") |>
            #DT::formatRound(c("tightness"),2) |>
            DT::formatRound(c("n","tail_n","tail_n_ended"),0) |>
            DT::formatRound(c("multimapping"),2) |>
            DT::formatStyle(names(df_show), "white-space"="nowrap")
    })
    
    output$site_info <- shiny::renderUI({
        req(selected())
        has_gene <- !is.na(selected()$gene_id)
        shiny::div(
            shiny::strong("Selected site ", selected()$site, " at ", selected()$location),
            shiny::br(),
            if (has_gene) paste0(selected()$name, " (", selected()$gene_id, ")"),
            if (has_gene) selected()$product,
            shiny::br(),
            shiny::br()
        )
    })
    
    plot_server("survival_plot", \() {
        p <- plot_km_survival(selected_kms()$km, selected_kms()$name, 
            min_tail=min_tail,
            max_tail=input$max_tail,
            cpm=input$cpm && input$show_samples,
            lib_sizes=selected_samples()$lib_size)
        if (input$show_samples) 
            p <- p + ggplot2::scale_color_discrete(type=selected_samples()$color)
        else
            p <- p + ggplot2::guides(color="none") + ggplot2::scale_color_discrete(type="black")
            
        print(p)
    })
    
    plot_server("density_plot", \() {
        p <- plot_km_density(selected_kms()$km, selected_kms()$name, 
            min_tail=min_tail,
            max_tail=input$max_tail,
            step=input$step,
            cpm=input$cpm && input$show_samples,
            lib_sizes=selected_samples()$lib_size)
        if (input$show_samples) 
            p <- p + ggplot2::scale_color_discrete(type=selected_samples()$color)
        else
            p <- p + ggplot2::guides(color="none") + ggplot2::scale_color_discrete(type="black")
        
        print(p)
    })
    
    plot_server("ridgeline_plot", \() {
        p <- plot_km_density_ridgeline(selected_kms()$km, selected_kms()$name, 
            min_tail=min_tail,
            max_tail=input$max_tail,
            step=input$step,
            cpm=input$cpm && input$show_samples,
            lib_sizes=selected_samples()$lib_size)
        if (input$show_samples) 
            p <- p + ggplot2::scale_fill_discrete(type=selected_samples()$color)
        else
            p <- p + ggplot2::scale_fill_discrete(type="black")
        
        print(p)
    })
    
    plot_server("density_heatmap", \() {
        p <- plot_km_density_heatmap(selected_kms()$km, selected_kms()$name, 
            min_tail=min_tail,
            max_tail=input$max_tail,
            step=input$step,
            normalize_max=FALSE,
            cpm=input$cpm && input$show_samples,
            lib_sizes=selected_samples()$lib_size)
        print(p)
    })
    
    plot_server("detail_plot", \() {
        ## Could do individual plots, but it's quite slow, and also hard to make out.
        #if (input$show_samples) {
        #    plots <- purrr::map2(selected_reads(), samples$sample, \(reads,name) 
        #        plot_tail_detail(reads) + ggplot2::labs(title=name))
        #    patchwork::wrap_plots(plots, ncol=ceiling(sqrt(length(plots)))) |> print()
        #} else {
            plot_tail_detail(selected_reads()) |> print()
        #}
    })
    
    plot_server("heatmap", \() { 
        req(input$table_rows_all)
        sites_wanted <- df()$site[ head(input$table_rows_all, input$heatmap_rows) ]
        this_sites <- sites |>
            dplyr::filter(site %in% .env$sites_wanted) |> 
            dplyr::collect()
        this_sites <- this_sites[match(sites_wanted, this_sites$site),]
        this_sites <- this_sites[!purrr::map_lgl(this_sites$tail_counts,is.null),]
        req(nrow(this_sites) > 0)
        
        kms <- purrr::map(this_sites$tail_counts, calc_km, assume_all_died=input$assume_all_died)
        names <- paste(tidyr::replace_na(this_sites$name,""), this_sites$site)
        
        p <- plot_km_density_heatmap(kms, names,
            min_tail=min_tail,
            max_tail=input$max_tail,
            step=input$step)
        print(p)
    })
    
    output$read_counts <- DT::renderDT({
        df <- selected_kms()
        df$tail_counts <- purrr::map(df$tail_counts, dplyr::collect)
        
        # Awkward
        if (input$show_samples) {
            counts <- purrr::map(selected_samples()$counts,\(item)
                dplyr::filter(item, site == .env$selected()$site) |> dplyr::collect())
            n <- purrr::map_dbl(counts, \(item) sum(item$n))
            n_reads <- purrr::map_dbl(counts, \(item) sum(item$n_read))
            n_multimapper_reads <- purrr::map_dbl(counts, \(item) sum(item$n_read_multimapper %||% NA))
            
            cpm <- n * 1e6 / selected_samples()$lib_size
        } else {
            n <- NULL
            n_reads <- NULL
            n_multimapper_reads <- NULL
            cpm <- NULL
        }
        
        result <- dplyr::tibble(
            name=df$name,
            tail90=purrr::map_dbl(df$km,km_quantile,0.9),
            tail50=purrr::map_dbl(df$km,km_quantile,0.5),
            tail10=purrr::map_dbl(df$km,km_quantile,0.1),
            cpm=cpm,
            n=n,
            n_tail=purrr::map_dbl(df$tail_counts,\(item) sum(item$n_event)),
            n_tail_ended=purrr::map_dbl(df$tail_counts,\(item) sum(item$n_died)),
            n_reads=n_reads,
            n_tail_reads=purrr::map_dbl(df$tail_counts,\(item) 
                if ("n_read_event" %in% names(item)) sum(item$n_read_event) else NA),
            n_multimapper_reads,
            multimapping=n_multimapper_reads / n_reads)
        
        DT::datatable(
            result,
            rownames=FALSE,
            options=list(pageLength=100)) |>
            DT::formatStyle(names(result), lineHeight="100%", paddingTop="5px", paddingBottom="5px") |>
            DT::formatRound(c("cpm","n","n_tail","n_tail_ended","tail90","tail50","tail10"),1) |>
            DT::formatRound(c("n_reads","n_tail_reads","n_multimapper_reads"),0) |>
            DT::formatRound(c("multimapping"),3) |>
            DT::formatStyle(
                'cpm',
                background = DT::styleColorBar(c(0,result$cpm)*1.1, '#aaaaff')) |>
            DT::formatStyle(
                c("tail90","tail50","tail10"),
                background = DT::styleColorBar(c(0,result$tail90,result$tail50,result$tail10)*1.1, "#aaffaa")) |>
            DT::formatStyle(
                c("n","n_tail","n_tail_ended"),
                background = DT::styleColorBar(
                    c(0,result$n,result$n_tail,result$n_tail_ended)*1.1, "#cccccc")) |>
            DT::formatStyle(
                c("multimapping"),
                background = DT::styleColorBar(c(0,1.1), "#ffaaaa"))
    })
}

