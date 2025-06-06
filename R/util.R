

# Some functions will either take a single item or a (simple) list of items.
# Turn single items into a list of length 1.
ensure_list <- function(x) {
    if (!"list" %in% class(x))
        x <- list(x)
    x
}

