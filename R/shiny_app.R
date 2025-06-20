

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
        title = title,
        theme = bslib::bs_theme(bootswatch="minty"),
        navbar_options = bslib::navbar_options(bg="#ddffee"),
        selected = "Sites",
        tests_panel,
        bslib::nav_panel(title="Sites", site_ui(tq, max_tail)))
}

tq_shiny_server <- function(input, output, session, tq, tests) {
    if (!is.null(tests)) {
        test_get <- test_server(input, output, session, tq=tq, tests=tests)
    } else {
        test_get <- list(sites_wanted=function() NULL)
    }
    
    site_server(input, output, session, tq=tq, get_sites_wanted=test_get$sites_wanted)
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