
#' TailQuant class
#'
methods::setClass("TailQuant", 
    slots = c(
        dir = "character",
        sites = "ANY", 
        samples = "ANY"
    )
)

methods::setMethod("show", "TailQuant", \(object) {
    cat(paste0(
        "TailQuant object\n",
        "@dir     ", object@dir, "\n",
        "@sites   ", nrow(object@sites), "\n",
        "@samples ", nrow(object@samples), "\n"))
})
