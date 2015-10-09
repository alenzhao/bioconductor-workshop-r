<!-- R Markdown Documentation, DO NOT EDIT THE PLAIN MARKDOWN VERSION OF THIS FILE -->

<!-- Copyright 2015 Google Inc. All rights reserved. -->

<!-- Licensed under the Apache License, Version 2.0 (the "License"); -->
<!-- you may not use this file except in compliance with the License. -->
<!-- You may obtain a copy of the License at -->

<!--     http://www.apache.org/licenses/LICENSE-2.0 -->

<!-- Unless required by applicable law or agreed to in writing, software -->
<!-- distributed under the License is distributed on an "AS IS" BASIS, -->
<!-- WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. -->
<!-- See the License for the specific language governing permissions and -->
<!-- limitations under the License. -->

# Analyzing Variants with BigQuery

Google Genomics can import variant calls from VCF files or Complete Genomics masterVar files so that you can query them with a simple API as we saw earlier in vignette [Working with Variants](http://bioconductor.org/packages/devel/bioc/vignettes/GoogleGenomics/inst/doc/AnnotatingVariants.html).  You can also export Variants to BigQuery for interactive analysis of these large datasets.  For more detail, see https://cloud.google.com/genomics/v1/managing-variants

In this example we will work with the [Illumina Platinum Genomes](http://googlegenomics.readthedocs.org/en/latest/use_cases/discover_public_data/platinum_genomes.html) dataset.

### Helpful BigQuery links

For this example, we'll also be working with [Google BigQuery](https://cloud.google.com/bigquery/). It's often helpful to have a [link to the docs](https://cloud.google.com/bigquery/what-is-bigquery) handy, and especially the [query reference](https://cloud.google.com/bigquery/query-reference).

## Run a query from R

The [bigrquery](https://github.com/hadley/bigrquery) package written by Hadley Wickham implements an R interface to [Google BigQuery](https://cloud.google.com/bigquery/).


```r
library(dplyr)
library(bigrquery)
```


```r
######################[ TIP ]########################################
## Set the Google Cloud Platform project id under which these queries will run.
## 
## If you are using the Google Bioconductor workshop docker image, this is already
## set for you in your .Rprofile and you can skip this step.

# project <- "YOUR-PROJECT-ID"
#####################################################################
```


```r
# Change the table here if you wish to run these queries against a different Variants table.
theTable <- "genomics-public-data:platinum_genomes.variants"
# theTable <- "genomics-public-data:1000_genomes.variants"
# theTable <- "genomics-public-data:1000_genomes_phase_3.variants"
```

Let's start by just counting the number of records in the table:

```r
querySql <- paste("SELECT COUNT(1) FROM [", theTable, "]", sep="")
querySql
```

```
## [1] "SELECT COUNT(1) FROM [genomics-public-data:platinum_genomes.variants]"
```

And send the query to the cloud for execution:

```r
result <- query_exec(querySql, project=project)
```

[bigrquery](https://github.com/hadley/bigrquery) uses the package [httr](https://github.com/hadley/httr) to perform OAuth.


```r
######################[ TIP ]########################################
## If you have any trouble with OAuth and need to redo/reset OAuth,
## run the following code.

# if(FALSE != getOption("httr_oauth_cache")) {
#  file.remove(getOption("httr_oauth_cache"))
#}
#message("Restart R to redo/reset OAuth.")
#####################################################################
```

And we see that the table has 688167235 rows - wow!

```r
result
```

```
##         f0_
## 1 688167235
```

## Run a query using the BigQuery Web User Interface

So what is actually in this table?  Click on [this link](https://bigquery.cloud.google.com/table/genomics-public-data:platinum_genomes.variants) to view the schema in the BigQuery web user interface.

We can also run the exact same query using the BigQuery web user interface.  In the BigQuery web user interface:

 (1) click on the *"Compose Query"* button
 (2) paste in the SQL for the query we just ran via R
 (3) click on *"Run Query"*.

## Run a query stored in a file from R

Instead of typing SQL directly into our R code, we can use a convenience function to read SQL from a file.

```r
library(GoogleGenomicsBioc2015Workshop)
DisplayAndDispatchQuery
```

```
function(queryUri, project, replacements=list()) {
  if (missing(queryUri)) {
    stop("Pass the file path or url to the file containing the query.")
  }
  if(missing(project)) {
    stop("Pass the project id of your Google Cloud Platform project.")
  }

  if (grepl("^https.*", queryUri)) {
    # Read the query from a remote location.
    querySql <- RCurl::getURL(queryUri, ssl.verifypeer=FALSE)
  } else {
    # Read the query from the local filesystem.
    querySql <- readChar(queryUri, nchars=1e6)
  }

  # If applicable, substitute values in the query template.
  for(replacement in names(replacements)) {
    querySql <- gsub(replacement, replacements[[replacement]], querySql, fixed=TRUE)
  }

  # Display the query to the terminal.
  cat(querySql)

  # Dispatch the query to BigQuery.
  bigrquery::query_exec(querySql, project)
}
<environment: namespace:GoogleGenomicsBioc2015Workshop>
```

This allows queries to be more easily shared among analyses and also reused for different datasets.  For example, in the following file we have a query that will retrieve data from any table exported from a Google Genomics Variant Set.

```r
file.show(file.path(system.file(package = "GoogleGenomicsBioc2015Workshop"),
                    "sql",
                    "variant-level-data-for-brca1.sql"))
```

Now let's run the query to retrieve variant data for BRCA1:

```r
result <- DisplayAndDispatchQuery(file.path(system.file(package = "GoogleGenomicsBioc2015Workshop"),
                                            "sql",
                                            "variant-level-data-for-brca1.sql"),
                                  project=project,
                                  replacements=list("_THE_TABLE_"=theTable))
```

```
# Retrieve variant-level information for BRCA1 variants.
SELECT
  reference_name,
  start,
  end,
  reference_bases,
  GROUP_CONCAT(alternate_bases) WITHIN RECORD AS alternate_bases,
  quality,
  GROUP_CONCAT(filter) WITHIN RECORD AS filter,
  GROUP_CONCAT(names) WITHIN RECORD AS names,
  COUNT(call.call_set_name) WITHIN RECORD AS num_samples,
FROM
  [genomics-public-data:platinum_genomes.variants]
WHERE
  reference_name CONTAINS '17' # To match both 'chr17' and '17'
  AND start BETWEEN 41196311 AND 41277499
# Skip non-variant segments if the source data was gVCF or CGI data
OMIT RECORD IF EVERY(alternate_bases IS NULL) OR EVERY(alternate_bases = '<NON_REF>')
ORDER BY
  start,
  alternate_bases
```
Number of rows returned by this query: 335.

Results from [bigrquery](https://github.com/hadley/bigrquery) are dataframes:

```r
mode(result)
```

```
## [1] "list"
```

```r
class(result)
```

```
## [1] "data.frame"
```

```r
summary(result)
```

```
##  reference_name         start               end          
##  Length:335         Min.   :41196407   Min.   :41196408  
##  Class :character   1st Qu.:41218836   1st Qu.:41218837  
##  Mode  :character   Median :41239914   Median :41239916  
##                     Mean   :41238083   Mean   :41238085  
##                     3rd Qu.:41258167   3rd Qu.:41258168  
##                     Max.   :41277186   Max.   :41277187  
##  reference_bases    alternate_bases       quality       
##  Length:335         Length:335         Min.   :   0.00  
##  Class :character   Class :character   1st Qu.:  48.56  
##  Mode  :character   Mode  :character   Median : 317.76  
##                                        Mean   : 445.50  
##                                        3rd Qu.: 737.80  
##                                        Max.   :6097.85  
##     filter             names            num_samples    
##  Length:335         Length:335         Min.   : 1.000  
##  Class :character   Class :character   1st Qu.: 2.000  
##  Mode  :character   Mode  :character   Median : 7.000  
##                                        Mean   : 5.304  
##                                        3rd Qu.: 7.000  
##                                        Max.   :17.000
```

```r
head(result)
```

```
##   reference_name    start      end reference_bases alternate_bases quality
## 1          chr17 41196407 41196408               G               A  733.47
## 2          chr17 41196820 41196822              CT               C   63.74
## 3          chr17 41196820 41196823             CTT            C,CT  314.59
## 4          chr17 41196840 41196841               G               T   85.68
## 5          chr17 41197273 41197274               C               A 1011.08
## 6          chr17 41197938 41197939               A              AT   86.95
##                                       filter names num_samples
## 1                                       PASS  <NA>           7
## 2                                      LowQD  <NA>           1
## 3                                       PASS  <NA>           3
## 4 TruthSensitivityTranche99.90to100.00,LowQD  <NA>           2
## 5                                       PASS  <NA>           7
## 6                                      LowQD  <NA>           3
```

## Visualize Query Results

The prior query was basically a data retrieval similar to what we performed earlier in this workshop when we used the GoogleGenomics Bioconductor package to retrieve data from the Google Genomics Variants API.  

But BigQuery really shines when it is used to perform an actual *analysis* - do the heavy-lifting on the big data resident in the cloud, and bring back the result of the analysis to R for further downstream analysis and visualization.

Let's do that now with a query that computes the Transition Transversion ratio for the variants within genomic region windows.

```r
result <- DisplayAndDispatchQuery(file.path(system.file(package = "GoogleGenomicsBioc2015Workshop"),
                                            "sql",
                                            "ti-tv-ratio.sql"),
                                  project=project,
                                  replacements=list("_THE_TABLE_"=theTable,
                                                    "_WINDOW_SIZE_"=100000))
```

```
# Compute the Ti/Tv ratio for variants within genomic region windows.
SELECT
  reference_name,
  window * 1e+05 AS window_start,
  transitions,
  transversions,
  transitions/transversions AS titv,
  num_variants_in_window,
FROM (
  SELECT
    reference_name,
    window,
    SUM(mutation IN ('A->G', 'G->A', 'C->T', 'T->C')) AS transitions,
    SUM(mutation IN ('A->C', 'C->A', 'G->T', 'T->G',
                     'A->T', 'T->A', 'C->G', 'G->C')) AS transversions,
    COUNT(mutation) AS num_variants_in_window
  FROM (
    SELECT
      reference_name,
      reference_bases,
      alternate_bases,
      INTEGER(FLOOR(start / 1e+05)) AS window,
      CONCAT(reference_bases, CONCAT(STRING('->'), alternate_bases)) AS mutation,
      COUNT(alternate_bases) WITHIN RECORD AS num_alts,
    FROM
      [genomics-public-data:platinum_genomes.variants]
    # Optionally add clause here to limit the query to a particular
    # region of the genome.
    #_WHERE_
    HAVING
      # Skip 1/2 genotypes _and non-SNP variants
      num_alts = 1
      AND reference_bases IN ('A','C','G','T')
      AND alternate_bases IN ('A','C','G','T'))
  GROUP BY
    reference_name,
    window)
ORDER BY
  reference_name,
  window_start
Retrieving data:  2.9s
```
Number of rows returned by this query: 28734.


```r
summary(result)
```

```
##  reference_name      window_start        transitions     transversions   
##  Length:28734       Min.   :        0   Min.   :   0.0   Min.   :   0.0  
##  Class :character   1st Qu.: 34600000   1st Qu.: 145.0   1st Qu.: 109.0  
##  Mode  :character   Median : 71000000   Median : 194.0   Median : 142.0  
##                     Mean   : 79795469   Mean   : 214.8   Mean   : 158.2  
##                     3rd Qu.:115400000   3rd Qu.: 252.0   3rd Qu.: 182.0  
##                     Max.   :249200000   Max.   :6313.0   Max.   :6703.0  
##                                                                          
##       titv       num_variants_in_window
##  Min.   :0.000   Min.   :    1.0       
##  1st Qu.:1.173   1st Qu.:  258.0       
##  Median :1.362   Median :  338.0       
##  Mean   :1.375   Mean   :  373.1       
##  3rd Qu.:1.557   3rd Qu.:  432.0       
##  Max.   :5.000   Max.   :12001.0       
##  NA's   :7
```

```r
head(result)
```

```
##   reference_name window_start transitions transversions      titv
## 1           chr1        0e+00         293           198 1.4797980
## 2           chr1        1e+05         147            76 1.9342105
## 3           chr1        2e+05          64            61 1.0491803
## 4           chr1        3e+05           2            11 0.1818182
## 5           chr1        4e+05          30            14 2.1428571
## 6           chr1        5e+05         207            76 2.7236842
##   num_variants_in_window
## 1                    491
## 2                    223
## 3                    125
## 4                     13
## 5                     44
## 6                    283
```

Since [bigrquery](https://github.com/hadley/bigrquery) results are dataframes, we can make use of all sorts of other great R packages to do our downstream work.  Here we use a few more packages from the Hadleyverse: dplyr for data filtering and ggplot2 for visualization.


```r
# Change this filter if you want to visualize the result of this analysis for a different chromosome.
chromosomeOneResults <- filter(result, reference_name == "chr1" | reference_name == "1")
```


```r
library(scales)
library(ggplot2)
```


```r
ggplot(chromosomeOneResults, aes(x=window_start, y=titv)) +
  geom_point() +
  stat_smooth() +
  scale_x_continuous(labels=comma) +
  xlab("Genomic Position") +
  ylab("Ti/Tv") +
  ggtitle("Ti/Tv by 100,000 base pair windows on Chromosome 1")
```

<img src="figure/titv-1.png" title="plot of chunk titv" alt="plot of chunk titv" style="display: block; margin: auto;" />

## Provenance

```r
sessionInfo()
```

```
R version 3.2.0 (2015-04-16)
Platform: x86_64-apple-darwin13.4.0 (64-bit)
Running under: OS X 10.10.5 (Yosemite)

locale:
[1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

attached base packages:
[1] stats4    parallel  stats     graphics  grDevices utils     datasets 
[8] methods   base     

other attached packages:
 [1] ggbio_1.16.0                           
 [2] TxDb.Hsapiens.UCSC.hg19.knownGene_3.1.2
 [3] GenomicFeatures_1.20.1                 
 [4] AnnotationDbi_1.30.1                   
 [5] Biobase_2.28.0                         
 [6] BSgenome.Hsapiens.UCSC.hg19_1.4.0      
 [7] BSgenome_1.36.0                        
 [8] rtracklayer_1.28.4                     
 [9] mgcv_1.8-6                             
[10] nlme_3.1-120                           
[11] ggplot2_1.0.1                          
[12] scales_0.2.5                           
[13] bigrquery_0.1.0                        
[14] dplyr_0.4.2                            
[15] GoogleGenomicsBioc2015Workshop_0.1     
[16] GoogleGenomics_1.1.3                   
[17] VariantAnnotation_1.14.1               
[18] GenomicAlignments_1.4.1                
[19] Rsamtools_1.20.4                       
[20] Biostrings_2.36.1                      
[21] XVector_0.8.0                          
[22] GenomicRanges_1.20.4                   
[23] GenomeInfoDb_1.4.0                     
[24] IRanges_2.2.2                          
[25] S4Vectors_0.6.0                        
[26] BiocGenerics_0.14.0                    

loaded via a namespace (and not attached):
 [1] httr_1.0.0           jsonlite_0.9.17      splines_3.2.0       
 [4] Formula_1.2-1        assertthat_0.1       latticeExtra_0.6-26 
 [7] RBGL_1.44.0          biovizBase_1.16.0    RSQLite_1.0.0       
[10] lattice_0.20-31      digest_0.6.8         RColorBrewer_1.1-2  
[13] colorspace_1.2-6     htmltools_0.2.6      httpuv_1.3.3        
[16] Matrix_1.2-0         plyr_1.8.3           OrganismDbi_1.10.0  
[19] XML_3.98-1.2         biomaRt_2.24.0       zlibbioc_1.14.0     
[22] BiocParallel_1.2.2   nnet_7.3-9           lazyeval_0.1.10     
[25] proto_0.3-10         survival_2.38-1      magrittr_1.5        
[28] mime_0.4             evaluate_0.7.2       GGally_0.5.0        
[31] MASS_7.3-40          foreign_0.8-63       graph_1.46.0        
[34] tools_3.2.0          formatR_1.2          stringr_1.0.0       
[37] munsell_0.4.2        cluster_2.0.1        lambda.r_1.1.7      
[40] futile.logger_1.4.1  grid_3.2.0           RCurl_1.95-4.6      
[43] dichromat_2.0-0      rstudioapi_0.3.1     rjson_0.2.15        
[46] rmarkdown_0.7        bitops_1.0-6         labeling_0.3        
[49] gtable_0.1.2         DBI_0.3.1            reshape_0.8.5       
[52] curl_0.9.3           markdown_0.7.7       reshape2_1.4.1      
[55] R6_2.1.1             gridExtra_0.9.1      knitr_1.10.5        
[58] Hmisc_3.16-0         futile.options_1.0.0 stringi_0.5-5       
[61] Rcpp_0.12.0          rpart_4.1-9          acepack_1.3-3.3     
```
