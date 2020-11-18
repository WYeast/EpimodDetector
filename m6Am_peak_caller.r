setwd("C:/Users/88/Desktop/yang_shi/")
source("rhtest.r")
counts <- read.csv("counts.new.csv",header=T)
head(counts)

##reads count
untreated_ip=counts$m6Am.Cont
untreated_input=counts$Input.Cont
treated_ip=counts$m6Am.KO
treated_input=counts$Input.KO

untreated_ip_total=4.9*10^6
untreated_input_total=1.8*10^6
treated_ip_total=1.6*10^6
treated_input_total=1.7*10^6
gene_id=counts$Geneid

rhtest_res=rhtest(untreated_ip, untreated_input, treated_ip, treated_input,
      untreated_ip_total, untreated_input_total, treated_ip_total,
      treated_input_total, minimal_count_fdr = 10)
head(rhtest_res)
res=cbind(gene_id,rhtest_res[,2],rhtest_res[,3],rhtest_res[,4])
write.csv(rhtest_res,"rhtest.csv")
