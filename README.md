# tailquant

**Under development, anything may change.**

* Default settings may change in future.

Improved quantification of poly(A) tail length from PAT-Seq data. Designed to be a replacement for 
[Tail Tools](https://github.com/Monash-RNA-Systems-Biology-Laboratory/tail-tools).

* Kaplan-Meier survival curve analysis (for older version of PAT-Seq).
* Faster to access data representation using parquet files.
* Calling poly(A) tail length from poly(T) sequence observed in read 2.


## Installation

tailquant is an R package intended to be used with Linux. It may also work with other Unix operating systems such as Mac OS.

First install:

* R (version 4.4 or higher)
* STAR (version 2.7.11b or higher)
* pigz

Then within R install the `BiocManager` and `remotes` packages, and use `BiocManager` to install this Git repository:

```r
install.packages(c("BiocManager", "remotes"))

BiocManager::install("Monash-RNA-Systems-Biology-Laboratory/tailquant")
```

## Usage after Tail Tools

tailquant can ingest the output of our earlier Tail Tools pipeline.

```r
library(tailquant)

# Enable multiprocessing.
# - Optional! Fragile!
# - Do not use from within RStudio, except by invoking R in the terminal.
# - future::multisession may also work, if multicore isn't available, or it may also fail.
future::plan(future::multicore, workers=8)

# Process an existing Tail Tools pipeline.
# Creates a new directory called my_output_dir
ingest_tt(
    out_dir="my_output_dir", 
    in_dir="my_tail_tools_pipeline_dir",
    tail_source="tt", # Where to get tail lengths from. "tt" means from Tail Tools output.
    site_pad=10,      # (default) Reads alignments ending +/-10 bases of a site are examined.
    min_tail=13,      # (default) Minimum "A"s to consider as having a poly(A) tail.
    length_trim=10    # (default) poly(A) tail reaching within 10 bases of the end 
                      #           we treat as not having ended.
)
```


## Demultiplexing and poly(T) length calling

tailquant can be used to demultiplex reads before using Tail Tools. It can then make use of poly(T) length information from read 2 if this is available.

We expect read pairs with:

* In read 1, sequence transcribed from the genome, possibly followed by a poly(A) tail sequence, and then possibly followed by UMI and barcode sequence.
* In read 2, a barcode, a UMI, and poly(T) sequence in read 2.

Due to the primers used, the poly(T) sequence is expected to be of length at least 12. Seeing a poly(T) of 13 bases or more in read 2 indicates we are seeing real poly(A) tail and not just primer sequence.

We also sometimes observe read pairs where there is a long span of "T"s in read 1. We think these are not real, and should be filtered. 

```r
library(tailquant)

# Enable multiprocessing.
# - Optional! Fragile!
# - Do not use from within RStudio, except by invoking R in the terminal.
# - future::multisession may also work, if multicore isn't available, or it may also fail.
future::plan(future::multicore, workers=6)

# ... 
# Create a data frame of samples with columns:
#     sample  - a sample name
#     barcode - barcode sequence (as seen in read 2)
# ...

# Have a look at your FASTQ file to check an appropriate quality cutoff.
# Ends of reads will be clipped at the point 1 in 5 non-G bases start falling below this quality.

# Step 1 loads all of the read pairs into a parquet file, 
# and calls samples, quality clipping, poly(A), and poly(T) lengths.
ingest_read_pairs(
     out_prefix="reads"   # <- A directory called reads will be created
    ,reads1="reads1.fastq.gz"
    ,reads2="reads2.fastq.gz"
    ,clip_quality_char="I" # <- Check this is appropriate
    ,samples=samples
    )

# You can now examine reads.peek.txt and reads.report.html to assess quality.

# Step 2 produces demultiplexed fastq files from the parquet file, suitable for Tail Tools.
demux_reads(
     out_dir="my_output_dir" 
    ,in_dir="reads"
    ,max_t_read_1=20   # (default) Discard any reads with a run of Ts longer than this in read 1
    #,min_t=13         # Optionally discard any read pairs without poly(T) in read 2 
    )

```

We've observed in read 1 there are often reads with long spans of "T"s. These are filtered by default by `demux_reads`. We would also always expect to see a poly(T) in read 2 corresponding the the poly(A) tail. `demux_reads` can also optionally filter reads that lack this.

`reads.parquet` contains the input reads, sample assignment, UMIs, and information about poly(A) and poly(T) lengths. This is a big file, so be sure to delete it once you are done with it.

`reads.peek.txt` shows a random selection of reads, and tailquant's interpretation of where the poly(A) and poly(T) tails are. `reads.report.html` is an HTML report giving a summary of tail identification, sample assignment, and UMI usage.

`demux_reads` produces fastq files suitable for Tail Tools. Tail Tools should be run with the clipping setting `adaptor="rc_umi_rc_barcode"`.

The poly(T) lengths can be used when ingesting Tail Tools output:

```r
# Process an existing Tail Tools pipeline, but use poly(T) lengths.
# Creates a new directory called my_output_dir
ingest_tt(
     out_dir="my_output_dir"
    ,in_dir="my_tail_tools_pipeline_dir"
    ,tail_source="read2" # Where to get tail lengths from.
    ,read_pairs_dir="reads"
    ,min_tail=13      # (default) Minimum "T"s to consider as having a poly(A) tail.
    ,length_trim=0     # poly(T) should be fine reaching to the end of read2.
    )
```


## Shiny app

```r
# Open Shiny app
tq <- load_tq("my_output_dir")
app <- tq_shiny(tq, title="Tail length examiner for this dataset")
app
```

The tailquant Shiny app supports a wide range of linear model based differential tests. 

Tests are specified in a named list:

* The name of each item will be used in a cache filenames.
* Each item is itself a named list containing:
    * "design" - matrix with row names corresponding to sample names - a design matrix as in limma etc.
    * "contrasts" - matrix with column names - a contrast matrix as in limma etc. Should have as many rows as there are columns in the design matrix. A single column specifies a t-test, and multiple columns specifies an F-test.
    * "title" - string - a title for the test.

For common experimental designs tailquant provides helper functions `make_tests_oneway` and `make_tests_twoway`.

For example, suppose your sample names were named like `group_replicate`, you could use:

```r
library(tidyverse)

# Use a regular expression to extract parts from the sample names
# and convert them to factors.
sample_names <- tq@samples$sample
parts <- str_match(sample_names, "(.*)_(.*)")
group <- fct_inorder(parts[,2])
rep <- fct_inorder(parts[,3])

# Construct tests for an experiment with one experimental factor.
# The batches argument is optional, only include it if you believe there is a batch effect.
tests <- make_tests_oneway(sample_names, groups=group, batches=rep)

app <- tq_shiny(tq, tests=tests, title="App with differential tests")
app
```

### Using tests with the old Tail Tools differential test app

Tests for tailquant can be converted to to work with the older Tail Tools software using the function `tests_to_tt()`.

```r
library(tailtools)

tt_app <- shiny_tests( 
    tests_to_tt(tests, pipeline_dir="...Tail Tools output directory..."),
    title="Tail Tools differential tests")

tt_app
```

