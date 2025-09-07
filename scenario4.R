#### Scenario 4
# Independent censoring
# Proportional hazard


#### Set working directory
setwd("C:\\Users\\Denis Talbot\\Dropbox\\Travail\\Recherche\\AF")


#### Source AF.iptw function
source("iptwAF.R");



#### Load required packages
require(averisk);
require(AF);



#### Determine the true value of the AF using Monte Carlo simulations of counterfactual
## Note : Data are generated according to a Weibull distribution
lambdaE = 0.00003;
aE = 2.1;
lambdaC = 0.00008;
aC = 2.4;
# set.seed(17893179);
# n = 5000000;
# L1 = runif(n, min = 20, max = 50);
# L2 = rbinom(n, 1, 0.5);
# A = rbinom(n, 1, plogis(-5 + 0.1*L1 + L2));
# TtE = 1/lambdaE**(1/aE)*(-log(1 - runif(n))/exp(0.5*A + 0.1*L1 + L2))**(1/aE);
# TtE_0 = 1/lambdaE**(1/aE)*(-log(1 - runif(n))/exp(0.5*0 + 0.1*L1 + L2))**(1/aE);
# AF_5 = 1 - mean(TtE_0 < 5)/mean(TtE < 5);
trueAF_5 = 0.2083791;
# I calibrated 
# - Sample size for men
# - The proportion of exposed ~ 25 %
# - Number of events among men ~ 500
# - Exposure effect size HR ~ 1.5



#### Prepare objects for the simulation
set.seed(7919317);
nrep = 1000;
R = 200;
n = 3000;
estimates.iptw = rep(NA, nrep);
cover.iptw = rep(NA, nrep);
estimates.trad = rep(NA, nrep);
cover.trad = rep(NA, nrep);
estimates.trad2 = rep(NA, nrep);
cover.trad2 = rep(NA, nrep);
estimates.trad.iptw = cover.trad.iptw = rep(NA, nrep);
estimates.trad2.iptw = cover.trad2.iptw = rep(NA, nrep);
estimates.averisk = cover.averisk = rep(NA, nrep);
estimates.AF = cover.AF = rep(NA, nrep);



#### Simulation
for(i in 1:nrep){
  # Generating the simulated data
  L1 = runif(n, min = 20, max = 50);
  L2 = rbinom(n, 1, 0.5);
  A = rbinom(n, 1, plogis(-5 + 0.1*L1 + L2));
  TtE = 1/lambdaE**(1/aE)*(-log(1 - runif(n))/exp(0.5*A + 0.1*L1 + L2))**(1/aE);
  TtC = 1/lambdaC**(1/aC)*(-log(1 - runif(n))/exp(0.5*A + 0.1*L1 + L2))**(1/aC);
  FollowUp = pmin(TtE, TtC);
  Event = 1*(TtE < TtC);

  # Proposed weigthed Kaplan-Meier
  res = AF.iptw(fu.time = FollowUp, event = Event, A = A,
                L = cbind(L1, L2), time = 5, R = R);
  estimates.iptw[i] = res$AF;
  cover.iptw[i] = res$LL < trueAF_5 & res$UL > trueAF_5;

  # Two traditional formulas
  modY = coxph(Surv(FollowUp, Event) ~ A + L1 + L2, ties = "breslow");
  HR = exp(coef(modY)[1]);
  LL.HR = exp(coef(modY)[1] - 1.96*sqrt(modY$var[1,1]));
  UL.HR = exp(coef(modY)[1] + 1.96*sqrt(modY$var[1,1]));
  pd = mean(A[Event == 1 & FollowUp < 5]);
  estimates.trad[i] = pd*(HR - 1)/HR;
  cover.trad[i] = pd*(LL.HR - 1)/LL.HR < trueAF_5 & pd*(UL.HR - 1)/UL.HR > trueAF_5;
  pe = mean(A);
  estimates.trad2[i] = pe*(HR - 1)/(pe*(HR - 1) + 1);
  cover.trad2[i] = pe*(LL.HR - 1)/(pe*(LL.HR - 1) + 1) < trueAF_5 &
                   pe*(UL.HR - 1)/(pe*(UL.HR - 1) + 1) > trueAF_5;

  # Modified traditional formulas where HR are estimated using IPTW
  modA = glm(A ~ L1 + L2, family = "binomial");
  w = A/modA$fitted + (1 - A)/(1 - modA$fitted);
  modY2 = coxph(Surv(FollowUp, Event) ~ A, weights = w, robust = TRUE);
  HR2 = exp(coef(modY2)[1]);
  LL.HR2 = exp(coef(modY2)[1] - 1.96*sqrt(modY2$var[1,1]));
  UL.HR2 = exp(coef(modY2)[1] + 1.96*sqrt(modY2$var[1,1]));
  estimates.trad.iptw[i] = pd*(HR2 - 1)/HR2;
  cover.trad.iptw[i] = pd*(LL.HR2 - 1)/LL.HR2 < trueAF_5 & pd*(UL.HR2 - 1)/UL.HR2 > trueAF_5;
  estimates.trad2.iptw[i] = pe*(HR2 - 1)/(pe*(HR2 - 1) + 1);
  cover.trad2.iptw[i] = pe*(LL.HR2 - 1)/(pe*(LL.HR2 - 1) + 1) < trueAF_5 &
                        pe*(UL.HR2 - 1)/(pe*(UL.HR2 - 1) + 1) > trueAF_5;

  # averisk package
  Event5 = Event*(FollowUp <= 5);
  avePAF = getAF(Event5 ~ A, cat_confounders = "L2", cont_confounders = "L1",
                 the.data = data.frame(Event5, A, L1, L2),
                 ci = TRUE, nsample_var = R, conf_level = 0.95);
  estimates.averisk[i] = avePAF[1];
  cover.averisk[i] = avePAF[2] < trueAF_5 & avePAF[3] > trueAF_5;

  # AF package
  AFmod = AFcoxph(modY, data = data.frame(FollowUp, Event, A, L1, L2),
                  exposure = "A", times = 5);
  estimates.AF[i] = AFmod$AF.est;
  cover.AF[i] = AFmod$AF.est - qnorm(0.975)*sqrt(AFmod$AF.var) < trueAF_5 &
                AFmod$AF.est + qnorm(0.975)*sqrt(AFmod$AF.var) > trueAF_5;

  print(data.frame(percent = i/nrep*100, Sys.time()));
}


mean(estimates.iptw - trueAF_5);
mean(cover.iptw);

mean(estimates.trad - trueAF_5);
mean(cover.trad);

mean(estimates.trad2 - trueAF_5);
mean(cover.trad2);

mean(estimates.trad.iptw - trueAF_5);
mean(cover.trad.iptw);

mean(estimates.trad2.iptw - trueAF_5);
mean(cover.trad2.iptw);

mean(estimates.averisk - trueAF_5);
mean(cover.averisk);

mean(estimates.AF - trueAF_5);
mean(cover.AF);

file.name = paste("results_scenario3_",gsub("-","",Sys.Date()),".Rdata", sep="");
save(estimates.iptw, cover.iptw, estimates.trad, cover.trad,
     estimates.trad2, cover.trad2, estimates.trad.iptw, cover.trad.iptw,
     estimates.trad2.iptw, cover.trad2.iptw, estimates.averisk,
     estimates.AF, cover.AF, file = file.name);



# > mean(estimates.iptw - trueAF_5);
# [1] -0.05757425
# > mean(cover.iptw);
# [1] 0.843
# > 
# > mean(estimates.trad - trueAF_5);
# [1] -0.005743018
# > mean(cover.trad);
# [1] 0.86
# > 
# > mean(estimates.trad2 - trueAF_5);
# [1] -0.04481047
# > mean(cover.trad2);
# [1] 0.776
# > 
# > mean(estimates.trad.iptw - trueAF_5);
# [1] -0.06443046
# > mean(cover.trad.iptw);
# [1] 0.669
# > 
# > mean(estimates.trad2.iptw - trueAF_5);
# [1] -0.1018648
# > mean(cover.trad2.iptw);
# [1] 0.363
# > 
# > mean(estimates.averisk - trueAF_5);
# [1] -0.06601073
# > mean(cover.averisk);
# [1] 0.967
# > 
# > mean(estimates.AF - trueAF_5);
# [1] -0.002266586
# > mean(cover.AF);
# [1] 0.952



