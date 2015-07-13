## ---- echo=FALSE, results="hide"-----------------------------------------
# Ensure that any errors cause the Vignette build to fail.
library(knitr)
opts_chunk$set(error=FALSE)

## ----message=FALSE-------------------------------------------------------
require(bigrquery)

## ----eval=FALSE----------------------------------------------------------
## ######################[ TIP ]########################################
## ## Set the Google Cloud Platform project id under which these queries will run.
## ##
## ## If you are using the Google Bioconductor workshop docker image, this is already
## ## set for you in your .Rprofile and you can skip this step.
## 
## # project <- "YOUR-PROJECT-ID"
## #####################################################################

## ------------------------------------------------------------------------
# Change the table here if you wish to run these queries against a different Variants table.
theTable <- "genomics-public-data:platinum_genomes.variants"
# theTable <- "genomics-public-data:1000_genomes.variants"
# theTable <- "genomics-public-data:1000_genomes_phase_3.variants"

## ------------------------------------------------------------------------
querySql <- paste("SELECT COUNT(*) FROM [", theTable, "]", sep="")
querySql

## ------------------------------------------------------------------------
result <- query_exec(querySql, project=project)

## ----eval=FALSE----------------------------------------------------------
## ######################[ TIP ]########################################
## ## If you have any trouble with OAuth and need to redo/reset OAuth,
## ## run the following code.
## 
## # if(FALSE != getOption("httr_oauth_cache")) {
## #  file.remove(getOption("httr_oauth_cache"))
## #}
## #message("Restart R to redo/reset OAuth.")
## #####################################################################

## ------------------------------------------------------------------------
result

## ------------------------------------------------------------------------
DisplayAndDispatchQuery <- function(queryUri, project, replacements=list()) {
  if (missing(queryUri)) {
    stop("Pass the file path or url to the file containing the query.")
  }
  if(missing(project)) {
    stop("Pass the project id of your Google Cloud Platform project.")
  }

  if (grepl("^https.*", queryUri)) {
    # Read the query from a remote location.
    querySql <- getURL(queryUri, ssl.verifypeer=FALSE)
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
  query_exec(querySql, project)
}

## ------------------------------------------------------------------------
file.show(file.path(system.file(package = "Bioc2015Workshop"),
                    "sql",
                    "variant-level-data-for-brca1.sql"))

## ----comment=NA----------------------------------------------------------
result <- DisplayAndDispatchQuery(file.path(system.file(package = "Bioc2015Workshop"),
                                            "sql",
                                            "variant-level-data-for-brca1.sql"),
                                  project=project,
                                  replacements=list("_THE_TABLE_"=theTable))

## ------------------------------------------------------------------------
mode(result)
class(result)
summary(result)
head(result)

## ----comment=NA----------------------------------------------------------
result <- DisplayAndDispatchQuery(file.path(system.file(package = "Bioc2015Workshop"),
                                            "sql",
                                            "ti-tv-ratio.sql"),
                                  project=project,
                                  replacements=list("_THE_TABLE_"=theTable,
                                                    "_WINDOW_SIZE_"=100000))

## ------------------------------------------------------------------------
summary(result)
head(result)

## ----message=FALSE-------------------------------------------------------
require(dplyr)

## ------------------------------------------------------------------------
# Change this filter if you want to visualize the result of this analysis for a different chromosome.
chromosomeOneResults <- filter(result, reference_name == "chr1" | reference_name == "1")

## ----message=FALSE-------------------------------------------------------
require(ggplot2)
require(scales)

## ----titv----------------------------------------------------------------
ggplot(chromosomeOneResults, aes(x=window_start, y=titv)) +
  geom_point() +
  stat_smooth() +
  scale_x_continuous(labels=comma) +
  xlab("Genomic Position") +
  ylab("Ti/Tv") +
  ggtitle("Ti/Tv by 100,000 base pair windows on Chromosome 1")

## ----provenance, comment=NA----------------------------------------------
sessionInfo()

