

# Some functions will either take a single item or a (simple) list of items.
# Turn single items into a list of length 1.
ensure_list <- function(x) {
    if (!"list" %in% class(x))
        x <- list(x)
    x
}

map_progress <- function(vec, message, func, ...) {
    env <- environment()
    cli::cli_progress_bar(message, total=length(vec))
    result <- purrr::map(vec, \(item) {
        result <- func(item, ...)
        cli::cli_progress_update(.envir=env)
        result
    })
    cli::cli_progress_done()
    result
}

#map_progress(1:5,"test",Sys.sleep)

