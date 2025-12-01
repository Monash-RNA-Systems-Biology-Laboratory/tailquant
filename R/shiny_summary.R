

summary_ui <- function(tq) {
    df <- tq_sample_stats(tq)
    
    levels <- c("3'UTR","3'Noncoding","Downstrand","CDS","Exon","Intron","Unassigned")
    site_df <- tq@sites |>
        dplyr::select(relation) |>
        dplyr::collect() |>
        dplyr::mutate(n=rowSums(tq_counts(tq))) |>
        dplyr::summarize(.by=relation, sites=dplyr::n(), n=sum(n)) |>
        dplyr::mutate(
            n_proportion=n/sum(n),
            relation=ifelse(is.na(relation),"Unassigned",relation)) |>
        dplyr::arrange(factor(relation, levels))
    
    sites_total <- sum(site_df$sites)
    
    shiny::div(
        shiny::fluidRow(
            shiny::column(width=6,
                shiny::h2(paste0(nrow(df), " samples"))),
            shiny::column(width=1, offset=5, style="text-align: right",
                shiny::downloadButton("samples_download", "CSV"))),
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
        shiny::p(shiny::br()),
        shiny::h2(paste0(sites_total, " sites")),
        DT::datatable(
                site_df,
                width="50%",
                rownames=FALSE,
                options=list(dom='t')) |>
            DT::formatStyle(names(site_df), lineHeight="100%", paddingTop="5px", paddingBottom="5px") |>
            DT::formatRound(c("sites","n"), 0) |>
            DT::formatPercentage("n_proportion", 2) |>
            DT::formatStyle(
                c("sites"),
                background = DT::styleColorBar(c(0,site_df$sites)*1.1, "#cccccc")) |>
            DT::formatStyle(
                c("n"),
                background = DT::styleColorBar(c(0,site_df$n)*1.1, "#cccccc")),
        shiny::p(shiny::br()),
        shiny::markdown("**Note:** Sites are classified from specific to broad categoires. Each site is only counted once in this table."))
}


summary_server <- function(input, output, session, tq) {
    output$samples_download <- shiny::downloadHandler(
        filename="samples.csv",
        content=\(filename) {
            readr::write_csv(tq_sample_stats(tq), filename)
        })
}