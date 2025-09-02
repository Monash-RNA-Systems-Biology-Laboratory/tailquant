

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


# Allocate a fixed proportion of margin on the left or right
# This is essentially a more manual version of what cowplot does.
plot_ensure_margin <- function(p, left=NA, right=NA) {
    if (is.na(left) && is.na(right))
        return(p)
    
    p <- cowplot::as_gtable(p)
    
    panel <- p$layout[ p$layout$name == "panel", ]
    assertthat::assert_that(nrow(panel) == 1)
    
    if (!is.na(left)) {
        j <- seq_len(panel$l-2) + 1
        p$widths[1] <- grid::unit(left,"npc") - sum(p$widths[j])
    }
    
    if (!is.na(right)) {
        n <- length(p$widths)
        j <- seq_len(n-panel$r-1) + panel$r
        p$widths[n] <- grid::unit(right,"npc") - sum(p$widths[j])
    }
    
    p
}