# EpimodDetector

EpimodDetector is based on the differential methylation analysis of control and m6Am methyltransferase (PCIF1) knockout conditions. 
In this package, differentially methylated transcription start sites are detected by hypergeometric test and defined as potential m6Am sites.

Installation:
Before the installation, devtools and Rsubread is required to be installed, if you haven't done that yet,
```
install.packages("devtools")

if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("Rsubread")

library(devtools)
install_github("WYeast/EpimodDetector")
```
Usage (default annotation is for human hg38 genome):
```
library(EpimodDetector)

m6Am_call(IP_BAM="IP.control.sorted.bam",
          INPUT_BAM = "Input.control.sorted.bam",
          KO_IP_BAM = "IP.PCIF1.KO.sorted.bam",
          KO_INPUT_BAM = "Input.PCIF1.KO.sorted.bam",
          IS_PAIRED_END = TRUE
          )
```         
Output file: a .csv file that lists fold change, log.p, log.fdr for reads in every TSS region across the different conditions
