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
