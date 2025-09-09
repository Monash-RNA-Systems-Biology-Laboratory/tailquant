
# QC and exploration shiny interface

rle_frame <- function(mat) {
    mat |> 
        reshape2::melt(varnames=c("feature","sample")) |> 
        dplyr::as_tibble() |>
        dplyr::mutate(median=median(value), .by=feature)
}


explore_ui <- function(tq) {
    
    configure <- shiny::tagList(
        shiny::numericInput("explore_tail_min", "Use reads with tail length at least", value=NA, min=13, step=1),
        shiny::numericInput("explore_const", "Constant in log2(CPM+constant)", value=1, min=0, step=0.5),
        shiny::numericInput("explore_n", "Sites to show in heatmap", value=50, min=10, max=1000, step=10))
    
    bslib::navset_underline(
        header=shiny::p(),
        bslib::nav_panel("Configure", configure),
        bslib::nav_panel("Heatmap", plot_ui("explore_heatmap", width=1000, height=800, margin_controls=FALSE)),
        bslib::nav_panel("RLE plot", plot_ui("explore_rle", width=1000, height=600, margin_controls=FALSE)),
        bslib::nav_panel("Sample vs median plot", plot_ui("explore_sample_median", width=1000, height=600, margin_controls=FALSE)))
}

explore_server <- function(input, output, session, tq) {
    title <- reactive({
        if (!is.na(input$explore_tail_min)) {
            paste0("Site expression from reads with tail ≥ ", input$explore_tail_min)
        } else {
            "Site expression from all reads"
        }
    })
    
    units <- reactive({ 
        paste0("log2(CPM+",input$explore_const,")") 
    })
    
    lcpm <- reactive({
        if (!is.na(input$explore_tail_min)) {
            counts <- tq_counts_tail_at_least(tq, input$explore_tail_min)
        } else {
            counts <- tq_counts(tq)
        }
        log2(edgeR::cpm(counts) + input$explore_const)
    })
    
    plot_server("explore_heatmap", \() {
        varistran::plot_heatmap(lcpm(), n=input$explore_n)
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