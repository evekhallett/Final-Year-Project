library(ggplot2)
library(relaimpo)
setwd('/Users/evehallett/Documents/dissertation/to_submit/Figure7/')
all_codons<- read.csv('/Users/evehallett/Documents/dissertation/to_submit/Figure7/files/goodman_gc3.csv')

m1 <- lm(Prot.FCC ~ c2 + c3 + c4 + c5 + c6 +
           c7 + c8 + c9 + c10 + c11, data = all_codons)
summary(m1)

bootresults<-boot.relimp(m1, b=1000, rela=TRUE) 
ci<-booteval.relimp(bootresults, norank=T)
print(ci)
pdf('/Users/evehallett/Documents/dissertation/to_submit/Figure7/relaimpo.pdf')
print(plot(ci))
dev.off()

summary(all_codons)

#initial experimentation - without bootstrapping 
#regular approach

#rel_r <- calc.relimp(m1, rela=TRUE)
#rel_infl <- rel_r$lmg[1]
#rel_infl <- rel_r$lmg[1]
#print(rel_infl)
#print(rel_r)
#print(rel_r$lmg)
#lmg_results <- data.frame(position = c(2:11),
#                          lmg_value = rel_r$lmg,
#                          stringsAsFactors = FALSE)
#print(lmg_results)

#absolute
#rel_nr <- calc.relimp(m1)
#rel_abs_infl <- rel_nr$lmg[1]
#lmg_results_abs <- data.frame(position = c(2:11),
#                          lmg_value = rel_nr$lmg,
#                          stringsAsFactors = FALSE)
#print(lmg_results_abs)


