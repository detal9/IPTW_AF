#### Loading required packages
require(mlogit);
require(boot);
require(survival);


#### Description
# This function estimates the fraction of the cases of a given outcome in a population
# occuring within a certain period of time that is attritributable to a given exposure (A)
# when the outcome is a time-to-event variable with possible right-censoring.  
#
# The function uses an inverse probability of treatment (or exposure) weighting estimator.
# This estimator is detailed in ****REF****. 
# The estimator assumes noninformative censoring.


#### Arguments
# dat         : (optionnal) A dataframe containing the variables fu.time, event, A and L.
# fu.time     : Either the name of the follow-up time variable as a character (if argument
#               dat is used) or a numeric variable containing the follow-up times.
# event       : Either the name of the event variable as a character (if argument
#               dat is used) or a numeric variable containing the event indicator.
#               event = 1 means that the follow-up ended with the event occurence,
#               event = 0 means that the time-to-event is right-censored at the
#               end of the follow-up.
# A           : Either the name of the exposure/treatment variable as a character 
#               (if argument dat is used) or a numeric variable containing the 
#               exposure/treatment variable. A should either be binary or categorical,
#               not continuous. 
# L           : Either a vector of character naming the potential confounding variables
#               (if argument dat is used) or a matrix/data frame containing the
#               potentially confounding covariates. If some variable should be treated
#               as categorical variables, they should be of type "factor".
#               It is possible to include nonlinear terms and interaction terms by
#               creating variables for these terms beforehand (****see example X****).
# times       : A numeric variable of length 1 indicating the follow-up time
#               at which to estimate the population attributable fraction
# Aref        : (optionnnal) The reference value for the exposure. If not provided, 0
#               is assumed to be the reference value.
# R           : (optional) Number of bootstrap replicates for computing confidence
#               intervals. The default is 0 (no confidence intervals are computed).
# trunc       : (optional) Percentile at which to truncate the weight distribution.
#               default = 1 (no truncation). Standard values are 0.99, 0.995 or 0.999.


#### AF.iptw function

AF.iptw = function(dat = NULL, fu.time, event, A, L, times, Aref = 0, R = 0, trunc = 1){
  #### Error checks
  if(!is.null(dat)){ # The user provided a dat argument
    if(!is.data.frame(dat)) stop("dat must be a data frame");
    if(!is.character(fu.time)) stop("When dat is provided, fu.time must be a character of length 1");
    if(length(fu.time) > 1) stop("When dat is provided, fu.time must be a character of length 1");
    if(!is.character(event)) stop("When dat is provided, event must be a character of length 1");
    if(length(event) > 1) stop("When dat is provided, event must be a character of length 1");
    if(!is.character(A)) stop("When dat is provided, A must be a character of length 1");
    if(length(A) > 1) stop("When dat is provided, A must be a character of length 1");
    if(!is.character(L)) stop("When dat is provided, L must be a character");

    # Getting variables from the dat
    fu.time = dat[, fu.time];
    event = dat[, event];
    A = dat[, A];

    # Creating the exposure model formula
    if(nlevels(as.factor(A)) == 2){ # A is binary
       Aform = paste("A ~", paste(L, collapse = " + "));
    } else{ # A is categorical
       Aform = paste("A ~ 1 |", paste(L, collapse = " + "));
    }

    L = dat[, L];

    # Creating a new dat
    dat2 = data.frame(fu.time = fu.time, event = event, A = A, L = L);
    L.names = paste0("L", 1:(ncol(dat2) - 3)); # new names of the covariates
    names(dat2)[4:ncol(dat2)] = L.names; # assign those new names


  } else { # The user did not provide a dat argument
    if(!is.numeric(fu.time)) stop("When dat is not provided, fu.time must be a numeric variable");
    if((!is.numeric(A)) & (!is.character(A))) stop("When dat is not provided, A must be a numeric or character variable");
    if(!is.numeric(event)) stop("When dat is not provided, event must be a numeric variable");
    if(!(is.matrix(L) | is.data.frame(L))) {
      stop("When dat is not provided, L must be a matrix or a data frame");
    }
    if((length(fu.time) != length(A)) | 
       (length(fu.time) != length(event)) |
       (length(fu.time) != nrow(L)) |
       (length(A) != length(event)) |
       (length(A) != nrow(L)) |
       (length(event) != nrow(L))) stop("fu.time, A, L and event should have the same number of rows");

    dat2 = data.frame(fu.time = fu.time, event = event, A = A, L = L); # Create a dataframe dat
    L.names = paste0("L", 1:(ncol(dat2) - 3)); # new names of the covariates
    names(dat2)[4:ncol(dat2)] = L.names; # assign those new names

    # Creating the exposure model formula
    if(nlevels(as.factor(A)) == 2){ # A is binary
       Aform = paste("A ~", paste(L.names, collapse = " + "));
    } else{ # A is categorical
       Aform = paste("A ~ 1 |", paste(L.names, collapse = " + "));
    }
  }
  if(min(fu.time <= 0)) stop("Only positive follow-up times are allowed...");
  if((!identical(sort(unique(event)), c(0, 1))) & (!identical(sort(unique(event)), 1))){
     stop("event should only include 0 and 1");
  }
  if(min(table(A)) == 1){ 
   stop("One level of A has a cell count of 1. Either A is continuous (not permitted) or data is too sparse...");
  }
  if(!(Aref %in% A)) stop("The reference value for A is not observed. Verify that the Aref is correctly specified.");
  if(!is.numeric(times)) stop("times must be a numeric of length 1 or more");
  if(times <= 0) stop("times must be > 0");
  if(!is.numeric(R)) stop("R must be numeric");
  if(length(time) != 1) stop("times must be of length 1");
  if(R < 0) stop("R must be >= 0");
  if(R %% 1 != 0) warning("R was not an integer and has been rounded up"); R = ceiling(R);
  if(!is.numeric(trunc)) stop("trunc should be numeric");
  if(min(trunc) <= 0 | max(trunc) >1) stop("trunc should be > 0 and <=1 and usually close to 1");


  #### Bootstrap function
  bootf = function(dat2, x){
    dat = dat2[x,]; # resampling

    # Extracting variables from dataframe
    A = dat$A; 
    event = dat$event;
    fu.time = dat$fu.time;
    L = dat[,4:ncol(dat)];
    
    #### Computing inverse probability of treatment weights
    if(nlevels(as.factor(A)) == 2){ # A is binary
      A = 1*(A != Aref); # Make A 0/1
      Aref = 0;
      modA = glm(Aform, family = "binomial", data = dat); # logistic regression
      predA = modA$fitted; # P(A = A_i|L)
      w = A/predA + (1 - A)/(1 - predA); # Weights
      w = pmin(w, quantile(w, trunc)); # Weight truncation
    } else { # A is categorical
      dat$id = 1:nrow(dat);
      ds = mlogit.data(data = dat, shape = "wide", choice = "A", varying = NULL, idvar = id); # mlogit data
      modA = mlogit(as.formula(Aform), data = ds); # Exposure model
      ps = data.frame(predict(modA, type = "probs", newdata = ds)) # propensity score vector
      ps = ps[order(dat$id),]; #reordering rows
      Indicator = matrix(NA, nrow = nrow(dat), ncol = nlevels(factor(A))); # Indicator matrix I(A = a)
      for(i in 1:nlevels(factor(A))){
        Indicator[,i] = (A == levels(factor(A))[i]);
      }
      w = rowSums(Indicator/ps); # Weights
      w = pmin(w, quantile(w, trunc)); # Weight truncation
    }

  
    #### Kaplan-Meier estimates of survival, observed and counterfactual
    kmObs = survfit(Surv(fu.time, event) ~ 1); # Observed K-M curve
    TObs = max(which(kmObs$time < times)); # Which estimate corresponds to times
    PY = 1 - kmObs$surv[TObs]; # Estimated observed event probability
    km0 = survfit(Surv(fu.time, event) ~ 1, weights = w, subset = A == Aref); # Counterfactual K-M curve
    T0 = max(which(km0$time < times)); # Which estimate corresponds to times
    PY0 = 1 - km0$surv[T0]; # Estimated counterfactual event probability
    FA = 1 - PY0/PY; # Estimated attributable fraction
    return(FA);
  }
 
  if(R == 0) {
    FA = bootf(dat2, 1:nrow(dat2));
    LL = NA; UL = NA;
  } else { # Bootstrap
    res = boot(dat2, bootf, R);
    FA = res$t0;
    LL = quantile(res$t, 0.025);
    UL = quantile(res$t, 0.975);  
  }
  

  #### Print and return results
  return(list(AF = FA, LL = LL, UL = UL));
} 








