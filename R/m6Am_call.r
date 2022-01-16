#' This is the description of m6Am_call function.
#' @title m6Am_call
#'
#' @description  This function detects m6Am sites through the comparision of samples from m6A-exo-seq data
#'
#' @details This function performs differential methylation analysis for samples under control and m6Am methyltransferase knockout conditions. Differentially methylated transcription start sites are defined as m6Am sites.
#'
#' @param IP_BAM: bam files for control ip sample
#' @param INPUT_BAM: bam files for control input sample
#' @param KO_IP_BAM: bam files for PCIF1 knockout ip sample
#' @param KO_INPUT_BAM: bam files for PCIF1 knockout input sample
#' @param GENE_ANNO_SAF: .saf file for the annotation of TSS regions
#' @return a .csv file that lists fold change, log.p, log.fdr for reads in every TSS region across the different conditions
#' @example
#' m6Am_call(IP_BAM="m6Am-Cont.sorted.bam",INPUT_BAM = "Input-Cont.sorted.bam",KO_IP_BAM = "m6Am-KO.sorted.bam",KO_INPUT_BAM = "Input-KO.sorted.bam",GENE_ANNO_SAF = "Homo_sapiens.GRCh38.95.gtf.tss.+-300.new.bed"))
#' @export
m6Am_call <- function(
  IP_BAM,INPUT_BAM,KO_IP_BAM,KO_INPUT_BAM,GENE_ANNO_SAF=NA){
  PARAMETERS=list();
  PARAMETERS$IP_BAM=IP_BAM
  PARAMETERS$INPUT_BAM=INPUT_BAM
  PARAMETERS$KO_IP_BAM=KO_IP_BAM
  PARAMETERS$KO_INPUT_BAM=KO_INPUT_BAM
  PARAMETERS$GENE_ANNO_SAF=GENE_ANNO_SAF
  PARAMETERS$IS_PAIRED_END=IS_PAIRED_END
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
    fc=featureCounts(file,annot.ext=PARAMETERS$GENE_ANNO_SAF,isPairedEnd=PARAMETERS$IS_PAIRED_END)
    counts=cbind(counts,fc$counts)
  }
  ##individual counts and total counts
  ip_count=counts[,1]
  input_count=counts[,2]
  ko_ip_count=counts[,3]
  ko_input_count=counts[,4]

  sum=apply(counts,2,sum)
  ip_total=as.numeric(sum(ip_count))
  input_total=as.numeric(sum(input_count))
  ko_ip_total=as.numeric(sum(ko_ip_count))
  ko_input_total=as.numeric(sum(ko_input_count))

  ##rhtest
  rhtest_res=rhtest(ip_count, input_count, ko_ip_count, ko_input_count,
                    ip_total, input_total, ko_ip_total,
                    ko_input_total, minimal_count_fdr = 10)
  ##write result
  rhtest_res=as.matrix(as.data.frame(rhtest_res))
  rhtest_res$gene_id=rownames(counts)
  write.csv(rhtest_res,"m6Am.csv",quote= F,sep="\t",col.names=F, row.names = F)
}
