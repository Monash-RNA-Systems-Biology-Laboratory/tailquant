
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

# local({
#     queue <- local_queue(2)
#     queue(print, 1)
#     queue(print, 2)
#     queue(print, 3)
#     print("Ok")
# })
