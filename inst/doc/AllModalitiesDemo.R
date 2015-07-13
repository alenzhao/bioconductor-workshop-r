## ---- echo=FALSE, results="hide"-----------------------------------------
# Ensure that any errors cause the Vignette build to fail.
library(knitr)
opts_chunk$set(error=FALSE)

## ----message=FALSE, comment=NA-------------------------------------------
pca_1kg <- read.table(file.path(system.file(package = "GoogleGenomicsBioc2015Workshop"), "extdata", "1kg-pca.tsv"), col.names=c("Sample", "PC1", "PC2"))

## ----pca, fig.align="center", fig.width=10, message=FALSE, comment=NA----
require(ggplot2)
ggplot(pca_1kg) +
  geom_point(aes(x=PC1, y=PC2)) +
  xlab("principal component 1") +
  ylab("principal component 2") +
  ggtitle("Principal Coordinate Analysis upon 1,000 Genomes")

## ----message=FALSE, warning=FALSE, comment=NA----------------------------
sample_info <- read.csv("http://storage.googleapis.com/genomics-public-data/1000-genomes/other/sample_info/sample_info.csv")
require(dplyr)
pca_1kg <- inner_join(pca_1kg, sample_info)

## ----pca-with-ethnicity, fig.align="center", fig.width=10, message=FALSE, comment=NA----
ggplot(pca_1kg) +
  geom_point(aes(x=PC1, y=PC2, color=Super_Population)) +
  xlab("principal component 1") +
  ylab("principal component 2") +
  ggtitle("Principal Coordinate Analysis upon 1,000 Genomes")

## ----eval=FALSE----------------------------------------------------------
## ######################[ TIP ]##################################
## ## Set the Google Cloud Platform project id under which these queries will run.
## ## If you are using the Bioconductor workshop docker image, this is already
## ## set for you in your .Rprofile.
## 
## # project <- "YOUR-PROJECT-ID"
## #####################################################################

## ----comment=NA, message=FALSE, warning=FALSE----------------------------
# Setup for BigQuery access
require(bigrquery)
DisplayAndDispatchQuery <- function(queryUri, replacements=list()) {
  querySql <- readChar(queryUri, nchars=1e6)
  cat(querySql)
  for(replacement in names(replacements)) {
    querySql <- sub(replacement, replacements[[replacement]], querySql, fixed=TRUE)
  }
  query_exec(querySql, project)
}

## ----message=FALSE, comment=NA-------------------------------------------
sample_alt_counts <- DisplayAndDispatchQuery(file.path(system.file(package = "GoogleGenomicsBioc2015Workshop"), "sql", "sample-alt-counts.sql"))

## ----alt-counts, fig.align="center", fig.width=10, message=FALSE, warning=FALSE, comment=NA----
sample_alt_counts <- inner_join(sample_alt_counts, sample_info)
require(scales) # for scientific_format()
ggplot(sample_alt_counts) +
  geom_point(aes(x=single, y=double, color=Super_Population)) +
  scale_x_continuous(label=scientific_format()) +
  scale_y_continuous(label=scientific_format()) +
  xlab("Variants with a single non-reference allele") +
  ylab("Variants with two non-reference alleles") +
  ggtitle("Heterozygosity Counts within 1,000 Genomes")

## ----message=FALSE, comment=NA-------------------------------------------
pca_1kg_brca1 <- read.table(file.path(system.file(package = "GoogleGenomicsBioc2015Workshop"), "extdata", "1kg-brca1-pca.tsv"), col.names=c("Sample", "PC1", "PC2"))

## ----brca1-pca, fig.align="center", fig.width=10, message=FALSE, comment=NA----
ggplot(pca_1kg_brca1) +
  geom_point(aes(x=PC1, y=PC2)) +
  xlab("principal component 1") +
  ylab("principal component 2") +
  ggtitle("Principal Coordinate Analysis upon BRCA1 within 1,000 Genomes")

## ----brca1-pca-with-gender, fig.align="center", fig.width=10, message=FALSE, warning=FALSE, comment=NA----
pca_1kg_brca1 <- inner_join(pca_1kg_brca1, sample_info)
ggplot(pca_1kg_brca1) +
  geom_point(aes(x=PC1, y=PC2, color=Gender)) +
  xlab("principal component 1") +
  ylab("principal component 2") +
  ggtitle("Principal Coordinate Analysis upon BRCA1 within 1,000 Genomes")

## ----brca1-pca-with-ethnicity, fig.align="center", fig.width=10, message=FALSE, comment=NA----
ggplot(pca_1kg_brca1) +
  geom_point(aes(x=PC1, y=PC2, color=Super_Population)) +
  xlab("principal component 1") +
  ylab("principal component 2") +
  ggtitle("Principal Coordinate Analysis upon BRCA1 within 1,000 Genomes")

## ----message=FALSE, comment=NA-------------------------------------------
pca_1kg_brca1 <- mutate(pca_1kg_brca1, 
                        case = 0 > PC1)

## ----brca1-pca-case-control, fig.align="center", fig.width=10, message=FALSE, comment=NA----
ggplot(pca_1kg_brca1) +
  geom_point(aes(x=PC1, y=PC2, color=Super_Population, shape=case), size=3) +
  xlab("principal component 1") +
  ylab("principal component 2") +
  ggtitle("Principal Coordinate Analysis upon BRCA1 within 1,000 Genomes")

## ----message=FALSE, comment=NA-------------------------------------------
case_sample_ids <- paste("'", filter(pca_1kg_brca1, case==TRUE)$Sample, "'", sep="", collapse=",")
result <- DisplayAndDispatchQuery(file.path(system.file(package = "GoogleGenomicsBioc2015Workshop"), "sql", "gwas-brca1-pattern.sql"),
                                  list(CASE_SAMPLE_IDS__=case_sample_ids))

## ----message=FALSE, comment=NA-------------------------------------------
head(result)

## ----message=FALSE, comment=NA-------------------------------------------
require(GoogleGenomics)

## ----eval=FALSE----------------------------------------------------------
## ######################[ TIP ]##################################
## ## You may be prompted to authenticate if you have not already done so.
## ## If you are using the Bioconductor workshop docker image, a package
## ## load hook will do this automagically.
## #
## # authenticate(file="YOUR/PATH/TO/client_secrets.json", invokeBrowser=FALSE)
## ###############################################################

## ----message=FALSE, comment=NA-------------------------------------------
top_results_sorted_by_start <- arrange(head(result, 20), start)
variants <- Reduce(c, apply(top_results_sorted_by_start,
                           1,
                           function(var) {
                             getVariants(datasetId="10473108253681171589",
                                         chromosome=as.character(var["reference_name"]),
                                         start=as.integer(var["start"]),
                                         end=as.integer(var["end"]))
                             }))
length(variants)

## ----message=FALSE, comment=NA-------------------------------------------
granges <- variantsToGRanges(variants)
granges

## ----comment=NA, message=FALSE-------------------------------------------
require(VariantAnnotation)
require(BSgenome.Hsapiens.UCSC.hg19)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
codingVariants <- locateVariants(granges, txdb, CodingVariants())
codingVariants

coding <- predictCoding(rep(granges, elementLengths(granges$ALT)),
                        txdb,
                        seqSource=Hsapiens,
                        varAllele=unlist(granges$ALT, use.names=FALSE))
coding

## ----message=FALSE, comment=NA-------------------------------------------
galignments <- getReads(readGroupSetId="CMvnhpKTFhDnk4_9zcKO3_YB", chromosome="17",
                     start=41218200, end=41218500, converter=readsToGAlignments)
galignments

## ----alignments, fig.align="center", fig.width=10, message=FALSE, warning=FALSE, comment=NA----
require(ggbio)
strand_plot <- autoplot(galignments, aes(color=strand, fill=strand))
coverage_plot <- ggplot(as(galignments, "GRanges")) + stat_coverage(color="gray40",
                                                      fill="skyblue")
tracks(strand_plot, coverage_plot, xlab="chr17")

## ----message=FALSE, comment=NA-------------------------------------------
sessionInfo()

