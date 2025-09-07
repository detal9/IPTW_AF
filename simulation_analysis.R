#### Set working directory
setwd("C:\\Users\\Denis Talbot\\Dropbox\\Travail\\Recherche\\AF")

cover.averisk = NA;

load("results_scenario0_20221026.Rdata");
bias1 = c(mean(estimates.trad), mean(estimates.trad.iptw),
          mean(estimates.trad2), mean(estimates.trad2.iptw),
          mean(estimates.averisk), mean(estimates.AF), mean(estimates.iptw)) - 0.1327762;
sd1 = c(sd(estimates.trad), sd(estimates.trad.iptw),
        sd(estimates.trad2), sd(estimates.trad2.iptw),
        sd(estimates.averisk), sd(estimates.AF), sd(estimates.iptw));  
cover1 = c(mean(cover.trad), mean(cover.trad.iptw),
           mean(cover.trad2), mean(cover.trad2.iptw),
           mean(cover.averisk), mean(cover.AF), mean(cover.iptw))*100;

load("results_scenario1_20221026.Rdata");
bias2 = c(mean(estimates.trad), mean(estimates.trad.iptw),
          mean(estimates.trad2), mean(estimates.trad2.iptw),
          mean(estimates.averisk), mean(estimates.AF), mean(estimates.iptw)) - 0.2083791;
sd2 = c(sd(estimates.trad), sd(estimates.trad.iptw),
        sd(estimates.trad2), sd(estimates.trad2.iptw),
        sd(estimates.averisk), sd(estimates.AF), sd(estimates.iptw));  
cover2 = c(mean(cover.trad), mean(cover.trad.iptw),
           mean(cover.trad2), mean(cover.trad2.iptw),
           mean(cover.averisk), mean(cover.AF), mean(cover.iptw))*100;

load("results_scenario2_20221026.Rdata");
bias3 = c(mean(estimates.trad), mean(estimates.trad.iptw),
          mean(estimates.trad2), mean(estimates.trad2.iptw),
          mean(estimates.averisk), mean(estimates.AF), mean(estimates.iptw)) - 0.2083791;
sd3 = c(sd(estimates.trad), sd(estimates.trad.iptw),
        sd(estimates.trad2), sd(estimates.trad2.iptw),
        sd(estimates.averisk), sd(estimates.AF), sd(estimates.iptw));  
cover3 = c(mean(cover.trad), mean(cover.trad.iptw),
           mean(cover.trad2), mean(cover.trad2.iptw),
           mean(cover.averisk), mean(cover.AF), mean(cover.iptw))*100;

load("results_scenario3_20221026.Rdata");
bias4 = c(mean(estimates.trad), mean(estimates.trad.iptw),
          mean(estimates.trad2), mean(estimates.trad2.iptw),
          mean(estimates.averisk), mean(estimates.AF), mean(estimates.iptw)) - 0.2083791;
sd4 = c(sd(estimates.trad), sd(estimates.trad.iptw),
        sd(estimates.trad2), sd(estimates.trad2.iptw),
        sd(estimates.averisk), sd(estimates.AF), sd(estimates.iptw));  
cover4 = c(mean(cover.trad), mean(cover.trad.iptw),
           mean(cover.trad2), mean(cover.trad2.iptw),
           mean(cover.averisk), mean(cover.AF), mean(cover.iptw))*100;

load("results_scenario4_20221027.Rdata");
bias5 = c(mean(estimates.trad), mean(estimates.trad.iptw),
          mean(estimates.trad2), mean(estimates.trad2.iptw),
          mean(estimates.averisk), mean(estimates.AF), mean(estimates.iptw)) - 0.1013821;
sd5 = c(sd(estimates.trad), sd(estimates.trad.iptw),
        sd(estimates.trad2), sd(estimates.trad2.iptw),
        sd(estimates.averisk), sd(estimates.AF), sd(estimates.iptw));  
cover5 = c(mean(cover.trad), mean(cover.trad.iptw),
           mean(cover.trad2), mean(cover.trad2.iptw),
           mean(cover.averisk), mean(cover.AF), mean(cover.iptw))*100;

res = data.frame(estimators = c("pd.cov", "pd.iptw", "pe.cov", "pe.iptw", "getAF", "AFcoxph", "IPTW.KM"),
                 bias1, sd1, cover1, bias2, sd2, cover2, bias3, sd3, cover3, bias4, sd4, cover4,
                 bias5, sd5, cover5);
res[, c(2,3,5,6,8,9,11,12,14,15)] =  res[, c(2,3,5,6,8,9,11,12,14,15)]*100;
res[, 2:16] = round(res[,2:16],1);   









