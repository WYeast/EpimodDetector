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
#' m6Am_call(IP_BAM="m6Am-Cont.sorted.bam",INPUT_BAM = "Input-Cont.sorted.bam",KO_IP_BAM = "m6Am-KO.sorted.bam",KO_INPUT_BAM = "Input-KO.sorted.bam",GENE_ANNO_SAF = GENE_ANNO_SAF, ))
#' @export
m6Am_call <- function(
  IP_BAM,INPUT_BAM,KO_IP_BAM,KO_INPUT_BAM,GENE_ANNO_SAF=NA,TSS_ANNO_SAF=NA,IS_PAIRED_END=NA,MINIMAL_COUNT_FDR=0){
  PARAMETERS=list();
  PARAMETERS$IP_BAM=IP_BAM
  PARAMETERS$INPUT_BAM=INPUT_BAM
  PARAMETERS$KO_IP_BAM=KO_IP_BAM
  PARAMETERS$KO_INPUT_BAM=KO_INPUT_BAM
  PARAMETERS$GENE_ANNO_SAF=GENE_ANNO_SAF
  PARAMETERS$TSS_ANNO_SAF=TSS_ANNO_SAF
  PARAMETERS$IS_PAIRED_END=IS_PAIRED_END
  PARAMETERS$MINIMAL_COUNT_FDR=MINIMAL_COUNT_FDR
  # check annotation, if not provided, use default hg38 annotation
  if (is.na(PARAMETERS$GENE_ANNO_SAF)) {
    PARAMETERS$GENE_ANNO_SAF = system.file("extdata", "hg38.exons.bed", package="EpimodDetector")}
  if (is.na(PARAMETERS$TSS_ANNO_SAF)) {
    PARAMETERS$TSS_ANNO_SAF = system.file("extdata", "hg38_start_100_cage_corrected.saf", package="EpimodDetector")}
  # check single/paired-end status
  if (is.na(PARAMETERS$IS_PAIRED_END)) {
    stop("must specify the single/paired-end status for m6Am_peak_caller to work!",
         call. = TRUE, domain = NULL)}
  # get bam files
  bam=c(PARAMETERS$IP_BAM,PARAMETERS$INPUT_BAM,PARAMETERS$KO_IP_BAM,PARAMETERS$KO_INPUT_BAM)
  no_bam_files=length(bam)

  #reads count
  merged_data=NULL

  for (ibam in 1:no_bam_files) {
    file=bam[ibam]
    #if sample is input, obtain gene level count
    if(ibam%%2==0){
      fc=featureCounts(file,annot.ext=PARAMETERS$GENE_ANNO_SAF,isPairedEnd=PARAMETERS$IS_PAIRED_END)
    }
    #if sample is ip, obtain TSS count
    else{
      fc=featureCounts(file,annot.ext=PARAMETERS$TSS_ANNO_SAF,isPairedEnd=PARAMETERS$IS_PAIRED_END)
    }
    fc=cbind(fc$annotation,fc$counts)
    if(ibam==1){
      merged_data=fc
    }
    else{
      merged_data=merge(x=merged_data,y=fc,by.x="GeneID",by.y="GeneID",all=FALSE)
    }
  }
  ##individual counts and total counts
  ip_count=merged_data[,7]
  input_count=merged_data[,13]
  ko_ip_count=merged_data[,19]
  ko_input_count=merged_data[,25]

  ip_total=as.numeric(sum(ip_count))
  input_total=as.numeric(sum(input_count))
  ko_ip_total=as.numeric(sum(ko_ip_count))
  ko_input_total=as.numeric(sum(ko_input_count))

  ##rhtest
  rhtest_res=rhtest(ip_count, input_count, ko_ip_count, ko_input_count,
                    ip_total, input_total, ko_ip_total,
                    ko_input_total, minimal_count_fdr = PARAMETERS$MINIMAL_COUNT_FDR)
  ##write result
  rhtest_res=as.data.frame(rhtest_res)
  rhtest_res=cbind(merged_data$GeneID,rhtest_res)
  colnames(rhtest_res)[1]="GeneID"
  write.csv(rhtest_res,"Diff.TSS.csv",row.names = FALSE)
}
