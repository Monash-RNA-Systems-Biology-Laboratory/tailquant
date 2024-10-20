
site_examiner_ui <- function(title) {
    ui <- shiny::fluidPage(
        shiny::titlePanel(title),
        shiny::sidebarLayout(
            shiny::sidebarPanel(
                shiny::numericInput("top_n", "List this many top sites by n_tail_ended (0=all)", value=0, min=0, step=1),
                shiny::numericInput("heatmap_rows", "Show this many sites in heatmap", value=50, min=1, max=2000, step=1),
                shiny::numericInput("step", "Density plot bin size", value=1, min=1, step=1),
                shiny::checkboxInput("show_samples", "Show individual samples in plots.", value=TRUE),
                shiny::checkboxInput("assume_all_died", "Use old method, treating tail-to-end-of-read as actual tail length.")
            ), 
            shiny::mainPanel(
                DT::DTOutput("table")
            )
        ),
        shiny::tabsetPanel(
            shiny::tabPanel("Mult-site heatmap", plot_ui("heatmap", width=1000, height=800)),
            shiny::tabPanel("Site reverse cumulative distribution", plot_ui("survival_plot", width=1000, height=600)),
            shiny::tabPanel("Site density", plot_ui("density_plot", width=1000, height=600)),
            shiny::tabPanel("Site heatmap", plot_ui("density_heatmap", width=1000, height=600)),
            shiny::tabPanel("Site read details", plot_ui("detail_plot", width=1000, height=800)),
            shiny::tabPanel("Site table", shiny::tableOutput("read_counts"))
        ),
        shiny::h2("Explanation"),
        shiny::p("n_tail = Number of reads with a poly(A) tail."),
        shiny::p("n_tail_ended = Number of reads with a poly(A) tail that ended before the end of the read."),
        shiny::p("tail90, tail50, tail10 = 90%/50%/10% of tails are longer than this. tail50 is the median tail length."),
        shiny::p("tightness = tail90/tail10. Sort by this column to see sites with tight or wide tail length distributions. Best to limit the list to some number of top sites by n_tail_ended first.")
    )
    
    ui
}

site_examiner_server <- function(tq, input,output,session) {
    sites <- tq$sites
    samples <- tq$samples
    
    if (!"color" %in% colnames(samples)) {
        samples$color <- scales::hue_pal()(nrow(samples))
    }
    
    
    df <- shiny::reactive({
        result <- sites |>
            dplyr::transmute(
                site, name, location, n_tail=n, n_tail_ended=n_died, 
                tail90, tail50, tail10, tightness=tail90/tail10, gene_id, biotype) |>
            dplyr::arrange(-n_tail_ended)
        if (input$top_n > 0)
            result <- dplyr::slice_head(result, n=input$top_n)
        dplyr::collect(result)
    })
    
    selected <- reactive({
        req(length(input$table_rows_selected) == 1)
        site <- df()$site[ input$table_rows_selected ]
        sites |>
            dplyr::filter(site == .env$site) |>
            dplyr::collect()
    })
    
    selected_kms <- reactive({
        site <- selected()$site
        
        if (input$show_samples) {
            names <- samples$sample
            tail_counts <- samples$tail_counts |>
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
        DT::datatable(
            df(),
            selection='single',
            rownames=FALSE, #width="100%", 
            class='compact cell-border hover',
            extensions='Buttons'
        ) |>
          DT::formatRound(c("tightness"),2)
    })
    
    plot_server("survival_plot", \() {
        p <- plot_km_survival(selected_kms()$km, selected_kms()$name, 
            min_tail=get_attr(sites, "min_tail", 0),
            max_tail=get_attr(sites, "max_tail"))
        if (input$show_samples) 
            p <- p + ggplot2::scale_color_discrete(type=samples$color)
        else
            p <- p + ggplot2::guides(color="none") + ggplot2::scale_color_discrete(type="black")
            
        print(p)
    })
    
    plot_server("density_plot", \() {
        p <- plot_km_density(selected_kms()$km, selected_kms()$name, 
            min_tail=get_attr(sites, "min_tail", 0),
            max_tail=get_attr(sites, "max_tail"),
            step=input$step)
        if (input$show_samples) 
            p <- p + ggplot2::scale_color_discrete(type=samples$color)
        else
            p <- p + ggplot2::guides(color="none") + ggplot2::scale_color_discrete(type="black")
        
        print(p)
    })
    
    plot_server("density_heatmap", \() {
        p <- plot_km_density_heatmap(selected_kms()$km, selected_kms()$name, 
            min_tail=get_attr(sites, "min_tail", 0),
            max_tail=get_attr(sites, "max_tail"),
            step=input$step,
            normalize_max=FALSE)
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
            min_tail=get_attr(sites, "min_tail", 0),
            max_tail=get_attr(sites, "max_tail"),
            step=input$step)
        print(p)
    })
    
    output$read_counts <- shiny::renderTable(digits=0, {
        df <- selected_kms()
        df$tail_counts <- purrr::map(df$tail_counts, dplyr::collect)
        dplyr::tibble(
            name=df$name,
            n_tail=purrr::map_dbl(df$tail_counts,\(item) sum(item$n_event)),
            n_tail_ended=purrr::map_dbl(df$tail_counts,\(item) sum(item$n_died)),
            tail90=purrr::map_dbl(df$km,km_quantile,0.9),
            tail50=purrr::map_dbl(df$km,km_quantile,0.5),
            tail10=purrr::map_dbl(df$km,km_quantile,0.1)
        )
    })
}


#' @export
shiny_site_examiner <- function(tq, title="Tail distribution examiner") {
    shiny::shinyApp(
        site_examiner_ui(title)
        ,\(...) site_examiner_server(tq, ...)
    )
}

