

tq_shiny_ui <- function(tq, title, tests, max_tail=NA) {
    max_tail_upper <- tq_tail_range(tq)[2]
    if (is.na(max_tail)) {
        max_tail <- max_tail_upper
    }
    
    tests_panel <- NULL
    if (!is.null(tests)) {
        tests_panel <- bslib::nav_panel(title="Tests", test_ui(tq, tests))
    }
    
    bslib::page_navbar(
        title=title,
        theme = bslib::bs_theme(bootswatch="cosmo"),
        navbar_options = bslib::navbar_options(bg="#ddffee"),
        selected = "Summary",
        bslib::nav_item("/"),
        bslib::nav_panel(title="Summary", summary_ui(tq)),
        bslib::nav_panel(title="QC & Exploration", explore_ui(tq)),
        tests_panel,
        bslib::nav_panel(title="Sites", site_ui(tq, max_tail)))
}


tq_shiny_server <- function(input, output, session, tq, tests) {
    
    summary_server(input, output, session, tq=tq)
    
    explore_server(input, output, session, tq=tq)
    
    if (!is.null(tests)) {
        test_get <- test_server(input, output, session, tq=tq, tests=tests)
    } else {
        test_get <- list(sites_wanted=function() NULL)
    }
    
    site_server(input, output, session, tq=tq) #, get_sites_wanted=test_get$sites_wanted)
    
    shiny::observe({
        want <- test_get$sites_wanted()
        if (length(want$name) == 1) {
            if (want$what == "genes")
                query <- paste0("gene_id=^",stringr::str_escape(want$name),"$")
            else
                query <- paste0("site=^",stringr::str_escape(want$name),"$")
            shiny::updateTextInput(session, "search", value=query)
        }
    })
}

#' Shiny tailquant interface
#'
#' @export
tq_shiny <- function(tq, title="Tail distribution examiner", tests=NULL, max_tail=NA) {
    shiny::shinyApp(
        tq_shiny_ui(tq, title=title, tests=tests, max_tail=max_tail),
        \(input, output, session) tq_shiny_server(input, output, session, tq=tq, tests=tests))
}

#' @export
shiny_site_examiner <- tq_shiny