

summary_ui <- function(tq) {
    df <- tq_sample_stats(tq)
    
    site_counts <- tq@sites |>
        dplyr::count(relation) |>
        dplyr::collect()
    sites_total <- sum(site_counts$n)
    sites_gene <- sum(site_counts$n[ !is.na(site_counts$relation) ])
    
    shiny::div(
        shiny::h2("Samples"),
        DT::datatable(
                df,
                rownames=FALSE,
                options=list(pageLength=100)) |>
            DT::formatStyle(names(df), lineHeight="100%", paddingTop="5px", paddingBottom="5px") |>
            DT::formatRound(c("n","n_tail","n_tail_ended","n_read","n_read_multimapper"), 0) |>
            DT::formatRound(c("mean_tail"), 1) |>
            DT::formatRound(c("reads_per_umi", "multimapping"), 3) |>
            DT::formatStyle(
                c("n","n_tail","n_tail_ended"),
                background = DT::styleColorBar(c(0,df$n)*1.1, "#cccccc")) |>
            DT::formatStyle(
                c("mean_tail"),
                background = DT::styleColorBar(c(0,df$mean_tail)*1.1, "#aaffaa")),
        shiny::p("n = Total UMIs"),
        shiny::h2("Sites"),
        shiny::p(paste0(sites_total, " sites.")),
        shiny::p(paste0(sites_gene, " sites associated with a gene.")))
}


summary_server <- function(input, output, session, tq) {
}