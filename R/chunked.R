
# Utilities to process files in chunks, including parquet

#' @export
local_write_parquet <- function(filename, callback) {
    sink <- NULL
    writer <- NULL
    
    yield <- function(df) {
        table <- arrow::Table$create(df)
        
        if (is.null(sink)) {
            schema <- table$schema
            sink <<- arrow::FileOutputStream$create(filename)
            writer <<- arrow::ParquetFileWriter$create(
                schema=schema,
                sink=sink,
                properties=arrow::ParquetWriterProperties$create(names(schema), version="2.4", compression="zstd"))
        }
        
        if (nrow(table) > 0)
            writer$WriteTable(table, nrow(table)) # Might want to tune chunk size...
    }
    
    withr::defer_parent({
        if (!is.null(sink)) {
            writer$Close()
            sink$close()
        }
    })
    
    yield
}


#' Scan a parquet file
#'
#' Note it may be necessary for your callback to invoke gc() to avoid excess memory usage.
#'
#' @export
scan_parquet <- function(filename, callback, columns=NULL, message="Processing") {
    #This doesn't seem to be needed:
    #source <- arrow::mmap_open(filename, "read")
    #withr::defer(source$close())
    #reader <- arrow::ParquetFileReader$create(source)
    
    # Using memory mapping seems to take up more and more memory, so disable
    reader <- arrow::ParquetFileReader$create(filename, mmap=FALSE)
    
    col_names <- names(reader$GetSchema())
    if (is.null(columns)) columns <- col_names
    column_indices <- match(columns, col_names)-1 # 0 based!
    assertthat::assert_that(!any(is.na(column_indices)))
    
    cli::cli_progress_bar(message, total=reader$num_row_groups)
    
    for(i in seq_len(reader$num_row_groups)-1L) {
        # Note: important not to store loaded data in a variable, so it can be garbage collected ASAP
        callback(dplyr::collect(reader$ReadRowGroup(i, column_indices=column_indices)))
        cli::cli_progress_update(1)
    }
}


#' @export
scan_query <- function(query, callback) {
    scanner <- arrow::Scanner$create(query)
    reader <- scanner$ToRecordBatchReader()
    repeat {
        item <- reader$read_next_batch()
        if (is.null(item)) break
        callback(dplyr::collect(item))
    }
}


#' Scan GAlignments chunks
#'
#' Note: If limit is larger than chunk size, it is rounded up to a multiple of chunk size.
#'
#' @export
scan_bam_chunks <- function(filename, param, callback, limit=NA, chunk=1e6) {
    if (is.na(limit))
        limit <- Inf
    chunk <- min(chunk, limit)
    
    result <- list()
    total <- 0
    
    f <- Rsamtools::BamFile(filename, yieldSize=chunk)
    open(f)
    withr::defer(close(f))
    
    while(total < limit) {
        #TODO: Allow callback to garbage collect alignments
        alignments <- GenomicAlignments::readGAlignments(f, param=param)
        callback(alignments)
        total <- total + length(alignments)
        if (length(alignments) == 0) break
        rm(alignments)
    }
}

#' @export
map_bam_chunks <- function(filename, param, callback, limit=NA, chunk=1e6) {
    if (is.na(limit))
        limit <- Inf
    chunk <- min(chunk, limit)
    
    result <- list()
    
    callback2 <- function(alignments) {
        result[[ length(result)+1 ]] <<- callback(alignments)
    }
    
    scan_bam_chunks(filename, param, callback2, limit, chunk)
    
    result
}


#' @export
local_queue <- function(workers = future::nbrOfWorkers()) {
    queue <- list()
    
    drain_one <- function() {
        do.call(queue[[1]][[1]], queue[[1]][[2]])
        queue <<- queue[-1]
    }
    
    withr::defer_parent({
        while(length(queue) > 0) 
            drain_one()
    })
    
    function(func, ...) {
        while(length(queue) >= workers) drain_one()
        queue[[ length(queue)+1 ]] <<- list(func, list(...))
    }
}

#' @export
parallel_walk <- function(vec, func, workers = future::nbrOfWorkers()) {
    queue <- local_queue()
    for(item in vec) {
        future_result <- future::future(seed=NULL, func(item))
        queue(future::value, future_result)
    }
}

# local({
#     queue <- local_queue(2)
#     queue(print, 1)
#     queue(print, 2)
#     queue(print, 3)
#     print("Ok")
# })
