

# Some functions will either take a single item or a (simple) list of items.
# Turn single items into a list of length 1.
ensure_list <- function(x) {
    if (!"list" %in% class(x))
        x <- list(x)
    x
}


# Save disk space and loading time for diagnostic plots
compact_plot <- function(p) {
    patchwork::wrap_elements(cowplot::as_grob(p))
}

