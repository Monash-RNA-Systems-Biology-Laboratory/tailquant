# tailquant

**Under development, may change.**

Improved quantification of poly(A) tail length from PAT-Seq data using Kaplan-Meier survival curves.

TODO: UMI handling.

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

