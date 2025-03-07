# tailquant

**Under development, may change.**

Improved quantification of poly(A) tail length from PAT-Seq data using Kaplan-Meier survival curves.

## Usage

```r
# Install and load tailquant package
library(tailquant)
# or
devtools::load_all("tailquant", export_all=FALSE)

# Process an existing Tail Tools pipeline.
# Creates a new directory called my_output_dir
ingest_tt(
    out_dir="my_output_dir", 
    in_dir="my_tail_tools_pipeline_dir",
    site_pad=10,    #(default) Reads alignments ending +/-10 bases of a site are examined
    min_tail=19,    #(default) Minimum "A"s to consider as having a poly(A) tail
    length_trim=10  #(default) poly(A) tail reaching within 10 bases of the end 
                    #          we treat as not having ended.
)

# Open Shiny app
tq <- load_tq("my_output_dir")
shiny_site_examiner(tq, title="Tail length examiner for this dataset")
```


## Demultiplexing and poly(A)/poly(T) length calling

```r
# Package needs to be installed for multiprocessing to work
library(tailquant)

# Optional! Fragile! Enable multiprocessing.
future::plan(future::multisession, workers=8)

# ... 
# Create a data frame of samples with columns:
#     sample  - a sample name
#     barcode - barcode sequence (as seen in read 2)
# ...

ingest_read_pairs(
    "reads.parquet",   # <- This file will be created
    "reads1.fastq.gz",
    "reads2.fastq.gz",
    samples=samples)

demux_reads(
    "my_output_dir", 
    "reads.parquet")
```

`reads.parquet` contains the input reads, sample assignment, UMIs, and information about poly(A) and poly(T) lengths. This is a big file, so be sure to delete it once you are done with it.

`demux_reads` produces fastq files suitable for Tail Tools. Tail Tools should be run with the clipping setting `adaptor="rc_umi_rc_barcode"`.
