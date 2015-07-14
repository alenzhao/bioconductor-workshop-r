# bioconductor-workshop-r
# 
R package containing instructional materials for using GoogleGenomics Bioconductor and bigrquery packages.

Depending upon how up-to-date your installation of Bioconductor is, the setup commands below *take quite a long time*.  See http://googlegenomics.readthedocs.org/en/latest/workshops/bioc-2015.html for a faster way to get started.

To install:
```
  # Install BiocInstaller.
  source("http://bioconductor.org/biocLite.R")
  # OPTIONAL: see http://www.bioconductor.org/developers/how-to/useDevel/
  # useDevel()
  # Install devtools which is needed for the special use of biocLite() below.
  biocLite("devtools")
  # Install the workshop material.
  biocLite("googlegenomics/bioconductor-workshop-r", build_vignettes=TRUE, dependencies=TRUE)
```

View and run the vignettes.
```
  help(package="GoogleGenomicsBioc2015Workshop")
```
