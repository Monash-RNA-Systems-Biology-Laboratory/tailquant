# tailquant

**Under development, may change.**

Improved quantification of poly(A) tail length from PAT-Seq data using Kaplan-Meier survival curves.

TODO: UMI handling.

## Usage

```r
# Install an load tailquant package
library(tailquant)
devtools::load_all("tailquant", export_all=FALSE)

# Process an existing Tail Tools pipeline.
# Creates a new directory called my_output_dir
ingest_tt("my_output_dir", "my_tail_tools_pipeline_dir")

# Open Shiny app
tq <- load_tq("my_output_dir")
shiny_site_examiner(tq, title="Tail length examiner for this dataset")
```

