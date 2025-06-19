
test_ui <- function(tq, tests) {
    choices <- names(tests)
    names(choices) <- purrr::imap_chr(tests, \(item, name) item[["title"]] %||% name)
    
    test_type_choices <- names(test_types)
    names(test_type_choices) <- purrr::map_chr(test_types, "title")
    
    configure <- shiny::tagList(
        shiny::selectizeInput("test_name", "Test", choices=choices, selected=NULL, width="75%", options=list()),
        shiny::selectizeInput("test_type", "Test type", choices=test_type_choices, selected=NULL, width="75%"),
        shiny::numericInput("test_min_count", "Minimum count", value=10, min=0, step=1),
        shiny::numericInput("test_min_count_in", "Minimum samples", value=2, min=1, step=1),
        shiny::numericInput("test_fdr", "False Discovery Rate", value=0.05, min=0.0, max=1.0, step=0.01),
        shiny::p("Site filtering: Sites are used if, after equalizing library sizes between samples, there are at least \"minimum samples\" samples having a count of at least \"minimum count\"."),
        shiny::p("FDR: The FDR specified here is used for confect calculation. The output also includes an FDR adjusted p-value for the null hypothesise that the effect is zero, which is unaffected by this."))
        
    
    summary <- shiny::tagList(
        shiny::uiOutput("test_summary"))
    
    results <- shiny::tagList(
        DT::DTOutput("test_table", fill=FALSE))
    
    diagnostics <- shiny::tagList(
        shiny::uiOutput("test_diagnostics"))
    
    shiny::tagList(
        # Not sure this is the best place for this...
        shiny::tags$style(type='text/css', ".selectize-dropdown-content {max-height: 600px; }"), 
        bslib::navset_underline(
            header=shiny::p(),
            bslib::nav_panel("Configure", configure),
            bslib::nav_panel("Result summary", summary),
            bslib::nav_panel("Result table", results),
            bslib::nav_panel("Result diagnostics", diagnostics)))
}

test_server <- function(input, output, session, tq, tests) {
    results <- reactive(withProgress(message="Testing", value=NA, {
        spec <- tests[[ input$test_name ]]
        spec$min_count <- input$test_min_count
        spec$min_count_in <- input$test_min_count_in
        spec$fdr <- input$test_fdr
        cache_key <- paste0(
            input$test_name,
            "_min",spec$min_count,
            "_in",spec$min_count_in,
            "_fdr",spec$fdr)
        tq_test(tq, input$test_type, cache_key=cache_key, spec=spec)
    }))
    
    output$test_summary <- shiny::renderUI({
        table <- results()$table
        shiny::div(
            shiny::h3(results()$title),
            shiny::p(nrow(table), " ", results()$what),
            shiny::p(sum(!is.na(table$confect)), "significant ", results()$what),
            shiny::p(format(results()$df_prior, digits=3), "prior degrees of freedom"))
    })
    
    output$test_table <- DT::renderDT(server=TRUE, {
        spec <- attr(results(),"version")
        contrasts <- colnames(spec$contrasts)
        
        df_show <- results()$table
        cols_keep <- setdiff(colnames(df_show), c(
                "index", "F", "dispersion_seen", "dispersion_used", "df", "df_seen", "df_used",
                paste0(contrasts,"_se")))
        df_show <- df_show[, cols_keep]
        
        digits <- max(0, -floor(log10(results()$step)))
        
        priority <- c("rank","name","confect","effect","se","fdr_zero","row_mean",contrasts)
        reordering <- order(match(colnames(df_show), priority))
        df_show <- df_show[,reordering,drop=FALSE]
        
        round_cols <- c("confect","effect","se","row_mean",contrasts)
        signif_cols <- colnames(df_show)[purrr::map_lgl(df_show, is.numeric)] |>
            setdiff(c("n_present","rank",round_cols))
        
        DT::datatable(
            df_show,
            selection='single',
            rownames=FALSE, #width="100%", 
            #class='compact cell-border hover',
            extensions='Buttons') |>
            DT::formatRound(round_cols, digits=digits) |>
            DT::formatSignif(signif_cols, digits=3)
    })
    
    output$test_diagnostics <- shiny::renderUI({
        spec <- attr(results(), "version")
        
        text <- capture.output({
            cat("design=")
            print(knitr::kable(spec$design))
            cat("\ncontrasts=")
            rownames(spec$contrasts) <- colnames(spec$design)
            print(knitr::kable(spec$contrasts))
        }) |> paste(collapse="\n")
        
        shiny::div(
            shiny::pre(text))
    })
    
    sites_wanted <- reactive({
        results()$table$name[ input$test_table_rows_selected ]
    })
    
    list(sites_wanted=sites_wanted)
}
