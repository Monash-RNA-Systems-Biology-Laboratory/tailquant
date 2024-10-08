
site_examiner_ui <- function(title) {
    ui <- shiny::fluidPage(
        shiny::titlePanel(title),
        #,shiny::sidebarLayout(shiny::div(), shiny::div())
        shiny::numericInput("step", "Density plot bin size", value=1, min=1, step=1),
        shiny::numericInput("min_n", "Minimum n_tail", value=1, min=1, step=1),
        shiny::numericInput("heatmap_rows", "Heatmap sites to show", value=50, min=1, step=1),
        shiny::checkboxInput("show_samples", "Show individual samples in plots.", value=TRUE),
        shiny::checkboxInput("assume_all_died", "Use old method, treating tail-to-end-of-read as actual tail length."),
        DT::DTOutput("table"),
        shiny::tabsetPanel(
            shiny::tabPanel("Mult-site heatmap", plot_ui("heatmap", width=1000, height=800)),
            shiny::tabPanel("Site reverse cumulative distribution", plot_ui("survival_plot", width=1000, height=600)),
            shiny::tabPanel("Site density", plot_ui("density_plot", width=1000, height=600)),
            shiny::tabPanel("Site heatmap", plot_ui("density_heatmap", width=1000, height=600)),
            shiny::tabPanel("Site read details", plot_ui("detail_plot", width=1000, height=800))
        ),
        shiny::div(style="height: 1000px")
    )
    
    ui
}

site_examiner_server <- function(tq, input,output,session) {
    sites <- dplyr::collect(tq$sites)
    samples <- tq$samples
    
    
    df <- shiny::reactive({
        sites |>
            dplyr::select(
                site, name, location, n_tail=n, n_tail_ended=n_died, 
                tail90, tail50, tail10, gene_id, biotype, product) |>
            dplyr::filter(n_tail >= .env$input$min_n)
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
                purrr::map(filter, site == .env$site)
        } else {
            names <- "Combined"
            tail_counts <- dplyr::filter(sites, site == .env$site)$tail_counts
        }
        
        tibble(
            name = names,
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
        ) #|>
          #  formatRound(c("cpm","median","iqr"),0) |>
          #  formatRound(c("bimodality_ratio","madmed","mean","sd"),3)
    })
    
    plot_server("survival_plot", \() {
        p <- plot_km_survival(selected_kms()$km, selected_kms()$name, 
            min_tail=get_attr(sites, "min_tail", 0),
            max_tail=get_attr(sites, "max_tail"))
        if (!input$show_samples) p <- p + ggplot2::guides(color="none")
        print(p)
    })
    
    plot_server("density_plot", \() {
        p <- plot_km_density(selected_kms()$km, selected_kms()$name, 
            min_tail=get_attr(sites, "min_tail", 0),
            max_tail=get_attr(sites, "max_tail"),
            step=input$step)
        if (!input$show_samples) p <- p + ggplot2::guides(color="none")
        print(p)
    })
    
    plot_server("density_heatmap", \() {
        p <- plot_km_density_heatmap(selected_kms()$km, selected_kms()$name, 
            min_tail=get_attr(sites, "min_tail", 0),
            max_tail=get_attr(sites, "max_tail"),
            step=input$step)
        print(p)
    })
    
    plot_server("detail_plot", \() {
        plot_tail_detail(selected_reads()) |> print()
    })
    
    plot_server("heatmap", \() { 
        req(input$table_rows_all)
        sites_wanted <- df()$site[ head(input$table_rows_all, input$heatmap_rows) ]
        this_sites <- sites[match(sites_wanted, sites$site),]
        
        kms <- purrr::map(this_sites$tail_counts, calc_km, assume_all_died=input$assume_all_died)
        names <- paste(replace_na(this_sites$name,""), this_sites$site)
        
        p <- plot_km_density_heatmap(kms, names,
            min_tail=get_attr(sites, "min_tail", 0),
            max_tail=get_attr(sites, "max_tail"),
            step=input$step)
        print(p)
    })
}


#' @export
shiny_site_examiner <- function(tq, title="Tail distribution examiner") {    
    shiny::shinyApp(
        site_examiner_ui(title)
        ,\(...) site_examiner_server(tq, ...)
    )
}

