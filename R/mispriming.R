

get_mispriming_str <- function(str, base, prop_a, switch_cost) {
    input <- as.integer(charToRaw(str) == charToRaw(base))
    
    non_a_cost <- prop_a / (1-prop_a)
    
    symbol_cost <- rbind(
        c( 0,           0),
        c( non_a_cost, -1))
    
    transition_cost <- rbind(
        c( 0,           switch_cost),
        c( switch_cost, 0))
    
    initial_cost <- c(0,0)
    final_cost <- c(0,0)
    
    states <- viterbi(input, symbol_cost, transition_cost, initial_cost, final_cost)
    states <- S4Vectors::Rle(states)
    ones <- S4Vectors::runValue(states) == 1
    dplyr::tibble(
        start=S4Vectors::start(states)[ones],
        end=S4Vectors::end(states)[ones])
}


#' @export
get_mispriming_regions <- function(genome, prop_a=0.6, switch_cost=6) {
    genome <- Biostrings::DNAStringSet(genome)
    
    if (is.null(names(genome)))
        names(genome) <- paste0("Seq", seq_len(length(genome)))
    
    result <- purrr::map_dfr(seq_len(length(genome)), \(i) {
        str <- as.character(genome[i])
        str <- stringr::str_to_upper(str)
        fwd <- get_mispriming_str(str, "A", prop_a=prop_a, switch_cost=switch_cost)
        fwd$seqnames <- names(genome)[i]
        fwd$strand <- "+"
        rev <- get_mispriming_str(str, "T", prop_a=prop_a, switch_cost=switch_cost)
        rev$seqnames <- names(genome)[i]
        rev$strand <- "-"
        dplyr::bind_rows(fwd, rev)
    })
    
    GenomicRanges::GRanges(result)
}
