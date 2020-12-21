#' This is the description of rhtest function.
#' @title rhtest
#'
#' @description  This function is the core funtion within m6Am_call
#'
#' @details This function performs differential methylation analysis through hypergeometric test.
#'
#' @param untreated_ip: bam files for control ip sample
#' @param untreated_input: bam files for control input sample
#' @param treated_ip: bam files for PCIF1 knockout ip sample
#' @param treated_input: bam files for PCIF1 knockout input sample
#' @param untreated_ip_total: depth of control ip sample
#' @param untreated_input_total: depth of control input sample
#' @param treated_ip_total: depth of PCIF1 knockout ip sample
#' @param treated_input_total: depth of PCIF1 knockout input sample
#' @return a list object that lists fold change, log.p, log.fdr for reads in every TSS region across the different conditions
#' @example
#' rhtest(untreated_ip, untreated_input, treated_ip, treated_input,untreated_ip_total, untreated_input_total, treated_ip_total,treated_input_total, minimal_count_fdr = 10)
#' @export

rhtest<-function (untreated_ip, untreated_input, treated_ip, treated_input,
          untreated_ip_total, untreated_input_total, treated_ip_total,
          treated_input_total, minimal_count_fdr = 10)
{
  if ((untreated_ip_total * treated_input_total) > (untreated_input_total *
                                                    treated_ip_total)) {
    if (untreated_ip_total > treated_input_total) {
      temp = (untreated_input_total * treated_ip_total)/treated_input_total
      untreated_ip = round(untreated_ip * temp/untreated_ip_total)
    }
    else {
      temp = (untreated_input_total * treated_ip_total)/untreated_ip_total
      treated_input = round(treated_input * temp/treated_input_total)
    }
  }
  if ((untreated_ip_total * treated_input_total) < (untreated_input_total *
                                                    treated_ip_total)) {
    if (untreated_input_total > treated_ip_total) {
      temp = round(untreated_ip_total * treated_input_total/treated_ip_total)
      untreated_input = round(untreated_input * temp/untreated_input_total)
    }
    else {
      temp = (untreated_ip_total * treated_input_total)/untreated_input_total
      treated_ip = round(treated_ip * temp/treated_ip_total)
    }
  }
  untreated_ip = pmax(1, untreated_ip)
  untreated_input = pmax(1, untreated_input)
  treated_ip = pmax(1, treated_ip)
  treated_input = pmax(1, treated_input)
  q = treated_ip
  m = treated_ip + untreated_ip
  n = untreated_input + treated_input
  k = treated_ip + treated_input
  log.hyper.p = phyper(q - 1, m, n, k, lower.tail = FALSE,
                       log.p = TRUE)
  log.hypo.p = phyper(q, m, n, k, lower.tail = TRUE, log.p = TRUE)
  log.fc = log((treated_ip/treated_input)/(untreated_ip/untreated_input))
  log.p = pmin(log.hyper.p, log.hypo.p)
  pvalues = exp(log.p)
  log.fdr = log(p.adjust(pvalues, method = "fdr"))
  m = untreated_ip + untreated_input + treated_ip + treated_input
  ID = which(m > minimal_count_fdr)
  log.fdr_sig = log(p.adjust(pvalues[ID], method = "fdr"))
  log.fdr[ID] = log.fdr_sig
  log.fdr = pmax(log.fdr, -1000)
  log.p = pmax(log.p, -1000)
  DIFF = list(log.fdr = log.fdr, log.p = log.p, log.fc = log.fc)
  return(DIFF)
}
