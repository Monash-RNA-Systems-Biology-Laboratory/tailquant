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
* samtools
* pigz

Then within R install the `BiocManager` and `remotes` packages, and use `BiocManager` to install this Git repository:

```r
install.packages(c("BiocManager", "remotes"))

BiocManager::install("Monash-RNA-Systems-Biology-Laboratory/tailquant")
```
