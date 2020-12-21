m6Am_call <- function(
  IP_BAM,INPUT_BAM,KO_IP_BAM,KO_INPUT_BAM,GENE_ANNO_SAF=NA){
  PARAMETERS=list();
  PARAMETERS$IP_BAM=IP_BAM
  PARAMETERS$INPUT_BAM=INPUT_BAM
  PARAMETERS$KO_IP_BAM=KO_IP_BAM
  PARAMETERS$KO_INPUT_BAM=KO_INPUT_BAM
  PARAMETERS$GENE_ANNO_SAF=GENE_ANNO_SAF
  # check annotation
  if (is.na(PARAMETERS$GENE_ANNO_SAF)) {
    stop("must specify the TSS annotation .saf file/TxDb object for m6Am_peak_caller to work!",
         call. = TRUE, domain = NULL)}
  # get bam files
  bam=c(PARAMETERS$IP_BAM,PARAMETERS$INPUT_BAM,PARAMETERS$KO_IP_BAM,PARAMETERS$KO_INPUT_BAM)
  no_bam_files=length(bam)
#for (ibam in 1:no_bam_files) {
#  file=bam[ibam]
#  if (! file.exists(paste(file,'.bai',sep=""))){
#    print(paste("Stage: index bam file", file))
#    indexBam(file)
#}}

#reads count
  counts=NULL
  for (ibam in 1:no_bam_files) {
    file=bam[ibam]
    fc=featureCounts(file,annot.ext=PARAMETERS$GENE_ANNO_SAF)
    counts=cbind(counts,fc$counts)
  }
  ##individual counts and total counts
  ip_count=counts[,PARAMETERS$IP_BAM]
  input_count=counts[,PARAMETERS$INPUT_BAM]
  ko_ip_count=counts[,PARAMETERS$KO_IP_BAM]
  ko_input_count=counts[,PARAMETERS$KO_INPUT_BAM]

  sum=apply(counts,2,sum)
  ip_total=as.numeric(sum[PARAMETERS$IP_BAM])
  input_total=as.numeric(sum[PARAMETERS$INPUT_BAM])
  ko_ip_total=as.numeric(sum[PARAMETERS$KO_IP_BAM])
  ko_input_total=as.numeric(sum[PARAMETERS$KO_INPUT_BAM])

  ##rhtest

  rhtest_res=rhtest(ip_count, input_count, ko_ip_count, ko_input_count,
                    ip_total, input_total, ko_ip_total,
                    ko_input_total, minimal_count_fdr = 10)
  ##write result
  write.csv(rhtest_res,"test.m6Am.call.csv")
}
