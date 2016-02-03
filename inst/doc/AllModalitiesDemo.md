<!-- R Markdown Documentation, DO NOT EDIT THE PLAIN MARKDOWN VERSION OF THIS FILE -->

<!-- Copyright 2014 Google Inc. All rights reserved. -->

<!-- Licensed under the Apache License, Version 2.0 (the "License"); -->
<!-- you may not use this file except in compliance with the License. -->
<!-- You may obtain a copy of the License at -->

<!--     http://www.apache.org/licenses/LICENSE-2.0 -->

<!-- Unless required by applicable law or agreed to in writing, software -->
<!-- distributed under the License is distributed on an "AS IS" BASIS, -->
<!-- WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. -->
<!-- See the License for the specific language governing permissions and -->
<!-- limitations under the License. -->

Data Analysis using Google Genomics
===================================

The following example makes use of the [Phase 1 variants](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20110521/README.phase1_integrated_release_version3_20120430) from the [1,000 Genomes Project](http://www.1000genomes.org/).  For more detail about how this data was loaded into the Google Genomics API, please see [Google Genomics Public Data](http://googlegenomics.readthedocs.org/en/latest/use_cases/discover_public_data/1000_genomes.html).

The VCFs comprising this dataset are **1.12 TB** when uncompressed and provide information about **39,706,715** variants for **1,092** individuals.

* [Working at Scale](#working-at-scale)
  * [Cluster Computing](#cluster-computing)
  * [Querying](#querying)
* [Zooming-In](#zooming-in)
  * [Simplistic GWAS](#simplistic-gwas)
  * [Annotate Variants with BioConductor](#annotate-variants-with-bioconductor)
* [Zooming-In Even Further](#zooming-in-even-further)
  * [Visualize Reads with BioConductor](#visualize-reads-with-bioconductor)
* [Provenance](#provenance)

Working at Scale
-------------------

### Cluster Computing

Suppose we have a new dataset.  One of the first things we might do is a basic visualization.  Let's start by projecting the relevant data into 2-dimensional space by performing a [Principal Coordinate Analysis](http://occamstypewriter.org/boboh/2012/01/17/pca_and_pcoa_explained/) based on the number of variants shared by each pair of individuals.

In this example we are reading in previously computed results, but its easy to spin up an [Apache Spark](http://spark.apache.org/) cluster on [Google Compute Engine](https://cloud.google.com/hadoop/what-is-hadoop) and run this analysis.


```r
library(GoogleGenomicsBioc2015Workshop)

# The directory in which the raw data files reside.
#extDataDir = "/PATH/TO/GIT/CLONE/OF/bioconductor-workshop-r/inst/extdata"
extdataDir = system.file("extdata", package = "GoogleGenomicsBioc2015Workshop")
```


```r
pca_1kg <- read.table(file.path(extdataDir, "1kg-pca.tsv"), col.names=c("Sample", "PC1", "PC2"))
```
This analysis performed an `O(N^2)` computation upon the relevant fields within the *terabyte* of data by running an [Apache Spark](http://spark.apache.org/) job which used the [Google Genomics Variants API](https://cloud.google.com/genomics/reference/rest/v1/variants) for its input.  See the Google Genomics [PCA cookbook entry](http://googlegenomics.readthedocs.org/en/latest/use_cases/compute_principal_coordinate_analysis/index.html) for implementation details and instructions as to how to run this job.

Visualizing the results, we see quite distinct clusters:

```r
library(ggplot2)
ggplot(pca_1kg) +
  geom_point(aes(x=PC1, y=PC2)) +
  xlab("principal component 1") +
  ylab("principal component 2") +
  ggtitle("Principal Coordinate Analysis upon 1,000 Genomes")
```

<img src="figure/pca-1.png" title="plot of chunk pca" alt="plot of chunk pca" style="display: block; margin: auto;" />

Let's pull in the [supplementary information](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20130606_sample_info/README_20130606_sample_info) we do have on these samples from [Google Cloud Storage](https://cloud.google.com/storage/):

```r
sample_info <- read.csv("http://storage.googleapis.com/genomics-public-data/1000-genomes/other/sample_info/sample_info.csv")
library(dplyr)
pca_1kg <- inner_join(pca_1kg, sample_info)
```

Applying sample ethnicity to the plot:

```r
ggplot(pca_1kg) +
  geom_point(aes(x=PC1, y=PC2, color=Super_Population)) +
  xlab("principal component 1") +
  ylab("principal component 2") +
  ggtitle("Principal Coordinate Analysis upon 1,000 Genomes")
```

<img src="figure/pca-with-ethnicity-1.png" title="plot of chunk pca-with-ethnicity" alt="plot of chunk pca-with-ethnicity" style="display: block; margin: auto;" />

we see that ethnicity appears to strongly correlate with the clusters.

### Querying

Let's also visualize a different aspect of this data by examining the variant counts for each individual of heterozygous reference variants (where one of the alleles is equal to the reference) and heterozygous alternate variants (where neither of the alleles is equal to the reference).


```r
# The directory in which the files containing SQL reside.
#sqlDir = "/PATH/TO/GIT/CLONE/OF/bioconductor-workshop-r/inst/sql"
sqlDir = system.file("sql", package = "GoogleGenomicsBioc2015Workshop")

######################[ TIP ]##################################
## Set the Google Cloud Platform project id under which these queries will run.
## If you are using the Bioconductor workshop docker image, this is already
## set for you in your .Rprofile.

# project <- "YOUR-PROJECT-ID"
#####################################################################
```


```r
sample_alt_counts <- DisplayAndDispatchQuery(file.path(sqlDir, "sample-alt-counts.sql"),
                                             project=project)
```

```
# Count alternate alleles for each sample.
SELECT
  Sample,
  SUM(single) AS single,
  SUM(double) AS double,
FROM (
  SELECT
    call.call_set_name AS Sample,
    SOME(call.genotype > 0) AND NOT EVERY(call.genotype > 0) WITHIN call AS single,
    EVERY(call.genotype > 0) WITHIN call AS double,
  FROM
    [genomics-public-data:1000_genomes.variants]
  OMIT RECORD IF
    reference_name IN ("X", "Y", "MT"))
GROUP BY
  Sample
ORDER BY
  Sample
```
Number of rows returned by this query: 1092.

This analysis performed an `O(N)` computation via [Google BigQuery](https://cloud.google.com/bigquery/).  Since BigQuery is a columnar data store, it scans only the columns referenced by the query.  In this case, 1 TB of data was scanned, typically within 10 seconds.

Visualizing the results, we again see quite distinct clusters:

```r
sample_alt_counts <- inner_join(sample_alt_counts, sample_info)
library(scales) # for scientific_format()
ggplot(sample_alt_counts) +
  geom_point(aes(x=single, y=double, color=Super_Population)) +
  scale_x_continuous(label=scientific_format()) +
  scale_y_continuous(label=scientific_format()) +
  xlab("Variants with a single non-reference allele") +
  ylab("Variants with two non-reference alleles") +
  ggtitle("Heterozygosity Counts within 1,000 Genomes")
```

<img src="figure/alt-counts-1.png" title="plot of chunk alt-counts" alt="plot of chunk alt-counts" style="display: block; margin: auto;" />

Zooming-In
------------------------

Suppose we are interested in examining variants within the BRCA1 gene.  We might run our PCoA a second time, zooming-in specifically to this region within the genome.

Again in this example we read in previously computed results, but since the amount of data over which we are computing is much less, it is feasible to run this Spark job on a local machine in just a few minutes.

```r
pca_1kg_brca1 <- read.table(file.path(extdataDir, "1kg-brca1-pca.tsv"), col.names=c("Sample", "PC1", "PC2"))
```

Examining this data visually:

```r
ggplot(pca_1kg_brca1) +
  geom_point(aes(x=PC1, y=PC2)) +
  xlab("principal component 1") +
  ylab("principal component 2") +
  ggtitle("Principal Coordinate Analysis upon BRCA1 within 1,000 Genomes")
```

<img src="figure/brca1-pca-1.png" title="plot of chunk brca1-pca" alt="plot of chunk brca1-pca" style="display: block; margin: auto;" />

we see distinct clusters with a structure much different than our former result upon the entire dataset.

Let's apply the sample information we do have (just gender and ethnicity) to this visualization to see if any of it appears to explain the clustering.  First, we'll try gender:

```r
pca_1kg_brca1 <- inner_join(pca_1kg_brca1, sample_info)
ggplot(pca_1kg_brca1) +
  geom_point(aes(x=PC1, y=PC2, color=Gender)) +
  xlab("principal component 1") +
  ylab("principal component 2") +
  ggtitle("Principal Coordinate Analysis upon BRCA1 within 1,000 Genomes")
```

<img src="figure/brca1-pca-with-gender-1.png" title="plot of chunk brca1-pca-with-gender" alt="plot of chunk brca1-pca-with-gender" style="display: block; margin: auto;" />

which has no apparent bearing on these variants.

Next, we'll try ethnicity:

```r
ggplot(pca_1kg_brca1) +
  geom_point(aes(x=PC1, y=PC2, color=Super_Population)) +
  xlab("principal component 1") +
  ylab("principal component 2") +
  ggtitle("Principal Coordinate Analysis upon BRCA1 within 1,000 Genomes")
```

<img src="figure/brca1-pca-with-ethnicity-1.png" title="plot of chunk brca1-pca-with-ethnicity" alt="plot of chunk brca1-pca-with-ethnicity" style="display: block; margin: auto;" />

which does appear to correlate to some amount of the clustering in the second principal component axis but not in the first principal component axis.

Let's split these individuals into two groups based on their position relative to the origin of the first principal component:

```r
pca_1kg_brca1 <- mutate(pca_1kg_brca1, 
                        case = 0 > PC1)
```
And visualize them again with their grouping:

```r
ggplot(pca_1kg_brca1) +
  geom_point(aes(x=PC1, y=PC2, color=Super_Population, shape=case), size=3) +
  xlab("principal component 1") +
  ylab("principal component 2") +
  ggtitle("Principal Coordinate Analysis upon BRCA1 within 1,000 Genomes")
```

<img src="figure/brca1-pca-case-control-1.png" title="plot of chunk brca1-pca-case-control" alt="plot of chunk brca1-pca-case-control" style="display: block; margin: auto;" />

### Simplistic GWAS

Next we perform a [simplistic GWAS](http://homes.cs.washington.edu/~suinlee/genome560/lecture7.pdf) on the BRCA1 variants to retrieve a ranked list of the variants that appear to differentiate these groups.

```r
case_sample_ids <- paste("'", filter(pca_1kg_brca1, case==TRUE)$Sample, "'", sep="", collapse=",")
result <- DisplayAndDispatchQuery(file.path(sqlDir, "gwas-brca1-pattern.sql"),
                                  project=project,
                                  replacements=list(CASE_SAMPLE_IDS__=case_sample_ids))
```

```
# A template for a simplistic GWAS query upon 1,000 Genomes phase 1 variants
# within BRCA1.  The template allows customization of the list of sample ids
# in the case group.  http://homes.cs.washington.edu/~suinlee/genome560/lecture7.pdf
SELECT
  reference_name,
  start,
  end,
  reference_bases,
  alternate_bases,
  vt,
  case_count,
  control_count,
  allele_count,
  ref_count,
  alt_count,
  case_ref_count,
  case_alt_count,
  control_ref_count,
  control_alt_count,
  # https://en.wikipedia.org/wiki/Yates%27s_correction_for_continuity
  ROUND(
    POW(ABS(case_ref_count - (ref_count/allele_count)*case_count) - 0.5,
      2)/((ref_count/allele_count)*case_count) +
    POW(ABS(control_ref_count - (ref_count/allele_count)*control_count) - 0.5,
      2)/((ref_count/allele_count)*control_count) +
    POW(ABS(case_alt_count - (alt_count/allele_count)*case_count) - 0.5,
      2)/((alt_count/allele_count)*case_count) +
    POW(ABS(control_alt_count - (alt_count/allele_count)*control_count) - 0.5,
      2)/((alt_count/allele_count)*control_count),
    3) AS chi_squared_score
FROM (
  SELECT
    reference_name,
    start,
    end,
    reference_bases,
    alternate_bases,
    vt,
    SUM(ref_count + alt_count) AS allele_count,
    SUM(ref_count) AS ref_count,
    SUM(alt_count) AS alt_count,
    SUM(IF(TRUE = is_case, INTEGER(ref_count + alt_count), 0)) AS case_count,
    SUM(IF(FALSE = is_case, INTEGER(ref_count + alt_count), 0)) AS control_count,
    SUM(IF(TRUE = is_case, ref_count, 0)) AS case_ref_count,
    SUM(IF(TRUE = is_case, alt_count, 0)) AS case_alt_count,
    SUM(IF(FALSE = is_case, ref_count, 0)) AS control_ref_count,
    SUM(IF(FALSE = is_case, alt_count, 0)) AS control_alt_count,
  FROM (
    SELECT
      reference_name,
      start,
      end,
      reference_bases,
      NTH(1, alternate_bases) WITHIN RECORD AS alternate_bases,
      vt,
      # 1000 genomes data is bi-allelic so there is only ever a single alt
      SUM(0 == call.genotype) WITHIN call AS ref_count,
      SUM(1 == call.genotype) WITHIN call AS alt_count,
      call.call_set_name IN ('HG00096','HG00099','HG00100','HG00102','HG00104','HG00106','HG00108','HG00112','HG00114','HG00118','HG00119','HG00121','HG00124','HG00125','HG00127','HG00128','HG00135','HG00136','HG00137','HG00142','HG00143','HG00149','HG00150','HG00152','HG00155','HG00156','HG00158','HG00160','HG00177','HG00179','HG00180','HG00182','HG00183','HG00185','HG00187','HG00188','HG00189','HG00190','HG00231','HG00232','HG00234','HG00236','HG00237','HG00238','HG00239','HG00243','HG00244','HG00245','HG00246','HG00247','HG00250','HG00255','HG00256','HG00257','HG00258','HG00260','HG00262','HG00264','HG00267','HG00268','HG00269','HG00270','HG00271','HG00273','HG00274','HG00275','HG00277','HG00278','HG00280','HG00285','HG00309','HG00311','HG00312','HG00315','HG00318','HG00321','HG00324','HG00325','HG00326','HG00329','HG00334','HG00335','HG00336','HG00339','HG00342','HG00345','HG00349','HG00350','HG00351','HG00360','HG00362','HG00364','HG00366','HG00367','HG00373','HG00375','HG00376','HG00377','HG00382','HG00383','HG00384','HG00403','HG00404','HG00406','HG00407','HG00419','HG00422','HG00436','HG00437','HG00446','HG00448','HG00451','HG00452','HG00457','HG00458','HG00463','HG00464','HG00472','HG00473','HG00476','HG00478','HG00525','HG00531','HG00534','HG00537','HG00542','HG00543','HG00554','HG00556','HG00557','HG00565','HG00577','HG00578','HG00580','HG00581','HG00583','HG00593','HG00595','HG00596','HG00614','HG00620','HG00626','HG00628','HG00634','HG00635','HG00637','HG00638','HG00650','HG00651','HG00654','HG00656','HG00663','HG00671','HG00683','HG00684','HG00692','HG00699','HG00701','HG00702','HG00704','HG00707','HG00708','HG00731','HG00736','HG00740','HG01048','HG01052','HG01060','HG01061','HG01069','HG01070','HG01075','HG01079','HG01082','HG01085','HG01095','HG01097','HG01104','HG01105','HG01107','HG01113','HG01125','HG01133','HG01134','HG01136','HG01137','HG01140','HG01149','HG01167','HG01173','HG01183','HG01188','HG01190','HG01191','HG01198','HG01251','HG01257','HG01271','HG01272','HG01274','HG01275','HG01342','HG01344','HG01345','HG01350','HG01351','HG01354','HG01356','HG01360','HG01365','HG01366','HG01374','HG01377','HG01378','HG01384','HG01389','HG01390','HG01437','HG01441','HG01455','HG01461','HG01462','HG01465','HG01488','HG01492','HG01494','HG01495','HG01498','HG01516','HG01518','HG01519','HG01521','HG01522','HG01550','HG01551','HG01617','HG01618','HG01619','HG01625','HG01626','NA06986','NA07000','NA07051','NA07056','NA10847','NA10851','NA11829','NA11843','NA11892','NA11894','NA11919','NA11920','NA11931','NA11932','NA11933','NA11992','NA11995','NA12003','NA12043','NA12044','NA12046','NA12058','NA12154','NA12155','NA12249','NA12272','NA12273','NA12282','NA12283','NA12340','NA12341','NA12347','NA12348','NA12546','NA12717','NA12748','NA12761','NA12763','NA12778','NA12812','NA12814','NA12815','NA12827','NA12829','NA12830','NA12843','NA12872','NA12873','NA12889','NA18501','NA18504','NA18520','NA18522','NA18525','NA18527','NA18528','NA18534','NA18535','NA18541','NA18542','NA18545','NA18547','NA18548','NA18549','NA18550','NA18552','NA18553','NA18559','NA18560','NA18562','NA18566','NA18567','NA18570','NA18572','NA18573','NA18574','NA18576','NA18582','NA18592','NA18593','NA18596','NA18597','NA18602','NA18603','NA18605','NA18609','NA18611','NA18613','NA18614','NA18615','NA18616','NA18617','NA18619','NA18624','NA18627','NA18628','NA18630','NA18631','NA18633','NA18635','NA18636','NA18637','NA18639','NA18641','NA18642','NA18740','NA18748','NA18757','NA18856','NA18868','NA18874','NA18907','NA18910','NA18916','NA18924','NA18939','NA18940','NA18944','NA18947','NA18948','NA18951','NA18956','NA18959','NA18960','NA18961','NA18962','NA18963','NA18971','NA18974','NA18977','NA18980','NA18981','NA18983','NA18984','NA18987','NA18988','NA18992','NA18998','NA18999','NA19003','NA19004','NA19005','NA19007','NA19009','NA19020','NA19028','NA19035','NA19041','NA19044','NA19056','NA19057','NA19058','NA19062','NA19063','NA19064','NA19072','NA19078','NA19079','NA19080','NA19081','NA19083','NA19095','NA19096','NA19098','NA19099','NA19118','NA19147','NA19152','NA19171','NA19172','NA19197','NA19198','NA19204','NA19207','NA19213','NA19257','NA19307','NA19308','NA19309','NA19311','NA19324','NA19327','NA19328','NA19355','NA19373','NA19374','NA19376','NA19379','NA19383','NA19393','NA19395','NA19399','NA19404','NA19428','NA19429','NA19430','NA19431','NA19435','NA19439','NA19444','NA19445','NA19451','NA19453','NA19461','NA19467','NA19468','NA19469','NA19472','NA19473','NA19625','NA19648','NA19652','NA19654','NA19655','NA19657','NA19660','NA19663','NA19664','NA19672','NA19676','NA19678','NA19682','NA19685','NA19701','NA19707','NA19716','NA19720','NA19723','NA19725','NA19726','NA19731','NA19732','NA19734','NA19737','NA19738','NA19741','NA19750','NA19752','NA19753','NA19758','NA19759','NA19761','NA19762','NA19764','NA19770','NA19777','NA19786','NA19788','NA19789','NA19794','NA19795','NA19818','NA19819','NA19835','NA19900','NA19901','NA19909','NA19914','NA19916','NA19920','NA19922','NA19982','NA20281','NA20282','NA20287','NA20289','NA20294','NA20322','NA20334','NA20336','NA20341','NA20344','NA20348','NA20356','NA20357','NA20359','NA20363','NA20414','NA20503','NA20504','NA20505','NA20506','NA20507','NA20510','NA20512','NA20513','NA20515','NA20516','NA20517','NA20518','NA20519','NA20520','NA20521','NA20524','NA20525','NA20527','NA20529','NA20531','NA20532','NA20534','NA20535','NA20538','NA20539','NA20540','NA20542','NA20544','NA20581','NA20589','NA20753','NA20757','NA20758','NA20759','NA20761','NA20765','NA20769','NA20770','NA20772','NA20773','NA20774','NA20778','NA20785','NA20786','NA20792','NA20795','NA20798','NA20804','NA20805','NA20806','NA20807','NA20808','NA20810','NA20811','NA20812','NA20813','NA20814','NA20815','NA20816','NA20819','NA20826') AS is_case,
    FROM
      [genomics-public-data:1000_genomes.variants]
    WHERE
      reference_name = '17'
      AND start BETWEEN 41196311 AND 41277499
    )
  GROUP BY
    reference_name,
    start,
    end,
    reference_bases,
    alternate_bases,
    vt)
WHERE
  # For chi-squared, expected counts must be at least 5 for each group
  (ref_count/allele_count)*case_count >= 5.0
  AND (ref_count/allele_count)*control_count >= 5.0
  AND (alt_count/allele_count)*case_count >= 5.0
  AND (alt_count/allele_count)*control_count >= 5.0
HAVING
  # Chi-squared critical value for df=1, p-value=5*10^-8 is 29.71679
  chi_squared_score >= 29.71679
ORDER BY
  chi_squared_score DESC,
  start
```
Number of rows returned by this query: 180.

Note that even though this query ran over a small region of the genome, with a minor change to the SQL we could have run this same GWAS query over all variants within a much larger region, over an entire chromosome, or even the full dataset; returning the ranked list of variants that differ between the two groups.

```r
head(result)
```

```
  reference_name    start      end reference_bases alternate_bases  vt
1             17 41218332 41218333               G               A SNP
2             17 41259048 41259049               C               T SNP
3             17 41261232 41261233               C               T SNP
4             17 41265775 41265776               A               G SNP
5             17 41268205 41268206               A               C SNP
6             17 41241389 41241390               C               A SNP
  case_count control_count allele_count ref_count alt_count case_ref_count
1       1158          1026         2184      1473       711            447
2       1158          1026         2184      1473       711            447
3       1158          1026         2184      1473       711            447
4       1158          1026         2184      1473       711            447
5       1158          1026         2184      1473       711            447
6       1158          1026         2184      1471       713            446
  case_alt_count control_ref_count control_alt_count chi_squared_score
1            711              1026                 0           931.230
2            711              1026                 0           931.230
3            711              1026                 0           931.230
4            711              1026                 0           931.230
5            711              1026                 0           931.230
6            712              1025                 1           929.544
```

### Annotate Variants with BioConductor

Now let's use the [GoogleGenomics R client](https://github.com/Bioconductor/GoogleGenomics) to retrieve the full records for the variants in which we are interested.

```r
library(GoogleGenomics)
```


```r
######################[ TIP ]##################################
## You may be prompted to authenticate if you have not already done so.
## If you are using the Bioconductor workshop docker image, a package
## load hook will do this automagically.
#
# authenticate(file="YOUR/PATH/TO/client_secrets.json", invokeBrowser=FALSE)
###############################################################
```


```r
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
```

```
[1] 20
```

We can also convert this data to [BioConductor](http://www.bioconductor.org/) datatypes such as [GRanges data type](http://www.bioconductor.org/packages/release/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesIntroduction.pdf).

```r
granges <- variantsToGRanges(variants)
granges
```

```
GRanges object with 20 ranges and 4 metadata columns:
             seqnames               ranges strand   |            REF
                <Rle>            <IRanges>  <Rle>   | <DNAStringSet>
   rs4793194    chr17 [41218333, 41218333]      *   |              G
   rs8176234    chr17 [41219780, 41219780]      *   |              T
   rs8176233    chr17 [41219804, 41219804]      *   |              T
   rs3950989    chr17 [41237953, 41237953]      *   |              G
   rs8176161    chr17 [41241390, 41241390]      *   |              C
         ...      ...                  ...    ... ...            ...
  rs12936316    chr17 [41263044, 41263044]      *   |              A
   rs8176109    chr17 [41265776, 41265776]      *   |              A
   rs8176098    chr17 [41268206, 41268206]      *   |              A
   rs8176092    chr17 [41270229, 41270229]      *   |              T
   rs8176088    chr17 [41270463, 41270463]      *   |              G
                            ALT      QUAL      FILTER
             <DNAStringSetList> <numeric> <character>
   rs4793194                  A       100        PASS
   rs8176234                  C       100        PASS
   rs8176233                  C       100        PASS
   rs3950989                  A       100        PASS
   rs8176161                  A       100        PASS
         ...                ...       ...         ...
  rs12936316                  G       100        PASS
   rs8176109                  G       100        PASS
   rs8176098                  C       100        PASS
   rs8176092                  G       100        PASS
   rs8176088                  A       100        PASS
  -------
  seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

This allows us to utilize the various BioConductor variant annotation packages:


```r
library(VariantAnnotation)
library(BSgenome.Hsapiens.UCSC.hg19)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
codingVariants <- locateVariants(granges, txdb, CodingVariants())
codingVariants
```

```
GRanges object with 22 ranges and 9 metadata columns:
      seqnames               ranges strand   | LOCATION  LOCSTART
         <Rle>            <IRanges>  <Rle>   | <factor> <integer>
    1    chr17 [41244000, 41244000]      *   |   coding      3335
    2    chr17 [41244000, 41244000]      *   |   coding      3407
    3    chr17 [41244000, 41244000]      *   |   coding      3548
    4    chr17 [41244000, 41244000]      *   |   coding      3548
    5    chr17 [41244000, 41244000]      *   |   coding      3548
  ...      ...                  ...    ... ...      ...       ...
   18    chr17 [41245466, 41245466]      *   |   coding      2082
   19    chr17 [41245466, 41245466]      *   |   coding      2082
   20    chr17 [41245466, 41245466]      *   |   coding      1941
   21    chr17 [41245466, 41245466]      *   |   coding      2004
   22    chr17 [41245466, 41245466]      *   |   coding      1194
         LOCEND   QUERYID        TXID                    CDSID      GENEID
      <integer> <integer> <character>            <IntegerList> <character>
    1      3335         8       63595 186231,186230,186233,...         672
    2      3407         8       63598 186231,186230,186233,...         672
    3      3548         8       63599 186231,186230,186233,...         672
    4      3548         8       63600 186231,186230,186233,...         672
    5      3548         8       63607 186231,186230,186233,...         672
  ...       ...       ...         ...                      ...         ...
   18      2082         9       63609 186231,186230,186233,...         672
   19      2082         9       63610 186231,186230,186233,...         672
   20      1941         9       63611 186231,186230,186233,...         672
   21      2004         9       63612 186231,186230,186233,...         672
   22      1194         9       63613 186231,186230,186233,...         672
            PRECEDEID        FOLLOWID
      <CharacterList> <CharacterList>
    1                                
    2                                
    3                                
    4                                
    5                                
  ...             ...             ...
   18                                
   19                                
   20                                
   21                                
   22                                
  -------
  seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

```r
coding <- predictCoding(rep(granges, elementLengths(granges$ALT)),
                        txdb,
                        seqSource=Hsapiens,
                        varAllele=unlist(granges$ALT, use.names=FALSE))
coding
```

```
GRanges object with 22 ranges and 16 metadata columns:
            seqnames               ranges strand   |            REF
               <Rle>            <IRanges>  <Rle>   | <DNAStringSet>
    rs16942    chr17 [41244000, 41244000]      -   |              T
    rs16942    chr17 [41244000, 41244000]      -   |              T
    rs16942    chr17 [41244000, 41244000]      -   |              T
    rs16942    chr17 [41244000, 41244000]      -   |              T
    rs16942    chr17 [41244000, 41244000]      -   |              T
        ...      ...                  ...    ... ...            ...
  rs1799949    chr17 [41245466, 41245466]      -   |              G
  rs1799949    chr17 [41245466, 41245466]      -   |              G
  rs1799949    chr17 [41245466, 41245466]      -   |              G
  rs1799949    chr17 [41245466, 41245466]      -   |              G
  rs1799949    chr17 [41245466, 41245466]      -   |              G
                           ALT      QUAL      FILTER      varAllele
            <DNAStringSetList> <numeric> <character> <DNAStringSet>
    rs16942                  C       100        PASS              G
    rs16942                  C       100        PASS              G
    rs16942                  C       100        PASS              G
    rs16942                  C       100        PASS              G
    rs16942                  C       100        PASS              G
        ...                ...       ...         ...            ...
  rs1799949                  A       100        PASS              T
  rs1799949                  A       100        PASS              T
  rs1799949                  A       100        PASS              T
  rs1799949                  A       100        PASS              T
  rs1799949                  A       100        PASS              T
                  CDSLOC    PROTEINLOC   QUERYID        TXID
               <IRanges> <IntegerList> <integer> <character>
    rs16942 [3335, 3335]          1112         8       63595
    rs16942 [3407, 3407]          1136         8       63598
    rs16942 [3548, 3548]          1183         8       63599
    rs16942 [3548, 3548]          1183         8       63600
    rs16942 [3548, 3548]          1183         8       63607
        ...          ...           ...       ...         ...
  rs1799949 [2082, 2082]           694         9       63609
  rs1799949 [2082, 2082]           694         9       63610
  rs1799949 [1941, 1941]           647         9       63611
  rs1799949 [2004, 2004]           668         9       63612
  rs1799949 [1194, 1194]           398         9       63613
                               CDSID      GENEID   CONSEQUENCE
                       <IntegerList> <character>      <factor>
    rs16942 186231,186230,186233,...         672 nonsynonymous
    rs16942 186231,186230,186233,...         672 nonsynonymous
    rs16942 186231,186230,186233,...         672 nonsynonymous
    rs16942 186231,186230,186233,...         672 nonsynonymous
    rs16942 186231,186230,186233,...         672 nonsynonymous
        ...                      ...         ...           ...
  rs1799949 186231,186230,186233,...         672    synonymous
  rs1799949 186231,186230,186233,...         672    synonymous
  rs1799949 186231,186230,186233,...         672    synonymous
  rs1799949 186231,186230,186233,...         672    synonymous
  rs1799949 186231,186230,186233,...         672    synonymous
                  REFCODON       VARCODON         REFAA         VARAA
            <DNAStringSet> <DNAStringSet> <AAStringSet> <AAStringSet>
    rs16942            AAA            AGA             K             R
    rs16942            AAA            AGA             K             R
    rs16942            AAA            AGA             K             R
    rs16942            AAA            AGA             K             R
    rs16942            AAA            AGA             K             R
        ...            ...            ...           ...           ...
  rs1799949            AGC            AGT             S             S
  rs1799949            AGC            AGT             S             S
  rs1799949            AGC            AGT             S             S
  rs1799949            AGC            AGT             S             S
  rs1799949            AGC            AGT             S             S
  -------
  seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

So a question for our users who have much experience in this domain: what should we examine next to determine potential explanations for the clustering we see?  Perhaps the relevant [ancestral haplotypes](http://hapmap.ncbi.nlm.nih.gov/originhaplotype.html)?

Zooming-in Even Further
------------------------

### Visualize Reads with BioConductor

We can also retrieve the reads from the [Genomics Reads API](https://cloud.google.com/genomics/reference/rest/v1/readgroupsets) for a given sample and examine coverage:

```r
galignments <- getReads(readGroupSetId="CMvnhpKTFhDnk4_9zcKO3_YB", chromosome="17",
                     start=41218200, end=41218500, converter=readsToGAlignments)
galignments
```

```
GAlignments object with 38 alignments and 1 metadata column:
                     seqnames strand       cigar    qwidth     start
                        <Rle>  <Rle> <character> <integer> <integer>
  ERR251040.39294893    chr17      -        100M       100  41218105
  ERR251040.23636053    chr17      -        100M       100  41218118
  ERR016214.20846952    chr17      -       4S77M        81  41218128
   ERR251039.9112219    chr17      +      69M31S       100  41218140
   ERR251040.9517733    chr17      +        100M       100  41218158
                 ...      ...    ...         ...       ...       ...
  ERR251039.29590756    chr17      -        100M       100  41218429
    ERR251039.668959    chr17      +        100M       100  41218465
   ERR016214.4338110    chr17      -      34S35M        69  41218474
  ERR251039.41699004    chr17      +        100M       100  41218484
   ERR016213.5228009    chr17      -         68M        68  41218496
                           end     width     njunc   |      flag
                     <integer> <integer> <integer>   | <numeric>
  ERR251040.39294893  41218204       100         0   |        19
  ERR251040.23636053  41218217       100         0   |       147
  ERR016214.20846952  41218204        77         0   |        19
   ERR251039.9112219  41218208        69         0   |        33
   ERR251040.9517733  41218257       100         0   |       163
                 ...       ...       ...       ... ...       ...
  ERR251039.29590756  41218528       100         0   |        19
    ERR251039.668959  41218564       100         0   |       163
   ERR016214.4338110  41218508        35         0   |       147
  ERR251039.41699004  41218583       100         0   |        35
   ERR016213.5228009  41218563        68         0   |       147
  -------
  seqinfo: 1 sequence from an unspecified genome; no seqlengths
```


```r
library(ggbio)
strand_plot <- autoplot(galignments, aes(color=strand, fill=strand))
coverage_plot <- ggplot(as(galignments, "GRanges")) + stat_coverage(color="gray40",
                                                      fill="skyblue")
tracks(strand_plot, coverage_plot, xlab="chr17")
```

<img src="figure/alignments-1.png" title="plot of chunk alignments" alt="plot of chunk alignments" style="display: block; margin: auto;" />

See also [GABrowse](http://gabrowse.appspot.com/#=&location=17%3A41218331&readsetId=CMvnhpKTFhDnk4_9zcKO3_YB&backend=GOOGLE) for an interactive Reads browser.

In summary, in this demo from the R prompt we were able to exercise both large scale and small scale data analysis using cloud-based infrastructure.

Provenance
-------------------
Lastly, let us capture version information about R and loaded packages for the sake of provenance.

```r
sessionInfo()
```

```
R version 3.2.2 (2015-08-14)
Platform: x86_64-apple-darwin13.4.0 (64-bit)
Running under: OS X 10.11.3 (El Capitan)

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
 [9] GoogleGenomics_1.1.3                   
[10] VariantAnnotation_1.14.1               
[11] GenomicAlignments_1.4.1                
[12] Rsamtools_1.20.4                       
[13] Biostrings_2.36.1                      
[14] XVector_0.8.0                          
[15] GenomicRanges_1.20.4                   
[16] GenomeInfoDb_1.4.0                     
[17] IRanges_2.2.2                          
[18] S4Vectors_0.6.0                        
[19] BiocGenerics_0.14.0                    
[20] scales_0.3.0                           
[21] bigrquery_0.1.0                        
[22] dplyr_0.4.3                            
[23] ggplot2_1.0.1                          
[24] GoogleGenomicsBioc2015Workshop_0.1     
[25] knitr_1.10.5                           

loaded via a namespace (and not attached):
 [1] Rcpp_0.12.0          biovizBase_1.16.0    lattice_0.20-33     
 [4] assertthat_0.1       digest_0.6.8         R6_2.1.1            
 [7] plyr_1.8.3           acepack_1.3-3.3      futile.options_1.0.0
[10] RSQLite_1.0.0        evaluate_0.7.2       httr_1.0.0          
[13] highr_0.5            zlibbioc_1.14.0      lazyeval_0.1.10     
[16] curl_0.9.3           rstudioapi_0.3.1     rpart_4.1-10        
[19] proto_0.3-10         labeling_0.3         splines_3.2.2       
[22] BiocParallel_1.2.2   foreign_0.8-65       stringr_1.0.0       
[25] RCurl_1.95-4.6       biomaRt_2.24.0       munsell_0.4.2       
[28] nnet_7.3-10          gridExtra_0.9.1      Hmisc_3.16-0        
[31] XML_3.98-1.2         reshape_0.8.5        MASS_7.3-43         
[34] bitops_1.0-6         RBGL_1.44.0          grid_3.2.2          
[37] jsonlite_0.9.17      GGally_0.5.0         gtable_0.1.2        
[40] DBI_0.3.1            magrittr_1.5         formatR_1.2         
[43] graph_1.46.0         stringi_0.5-5        reshape2_1.4.1      
[46] latticeExtra_0.6-26  futile.logger_1.4.1  Formula_1.2-1       
[49] RColorBrewer_1.1-2   rjson_0.2.15         lambda.r_1.1.7      
[52] tools_3.2.2          dichromat_2.0-0      OrganismDbi_1.10.0  
[55] survival_2.38-3      colorspace_1.2-6     cluster_2.0.3       
```
