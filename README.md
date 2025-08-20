# tailquant

**Under development, may change.**

Improved quantification of poly(A) tail length from PAT-Seq data using Kaplan-Meier survival curves.

## Usage

```r
# Package needs to be installed for multiprocessing to work
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

# Open Shiny app
tq <- load_tq("my_output_dir")
tq_shiny(tq, title="Tail length examiner for this dataset")
```


## Demultiplexing and poly(A)/poly(T) length calling

```r
# Package needs to be installed for multiprocessing to work
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
# and calls samples, quality clipping poly(A) and poly(T) lengths
ingest_read_pairs(
    out_file="reads.parquet",   # <- This file will be created
    reads1="reads1.fastq.gz",
    reads2="reads2.fastq.gz",
    clip_quality_char="I", # <- Check this is appropriate
    samples=samples)

# Step 2 produces demultiplexed fastq files from the parquet file.
demux_reads(
    out_dir="my_output_dir", 
    in_file="reads.parquet")

# Optionally in step 2 you can filter by poly(T) length
demux_reads(
    out_dir="my_output_dir", 
    in_file="reads.parquet",
    min_t=12)

# You can examine tailquant's interpretation of a random selection of reads with:
reads_peek(
    in_file="reads.parquet",
    out_file="peek.txt")

```

`reads.parquet` contains the input reads, sample assignment, UMIs, and information about poly(A) and poly(T) lengths. This is a big file, so be sure to delete it once you are done with it.

`demux_reads` produces fastq files suitable for Tail Tools. Tail Tools should be run with the clipping setting `adaptor="rc_umi_rc_barcode"`.
