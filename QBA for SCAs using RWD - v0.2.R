#######################
###
### SETUP
###
#######################

# load libraries

library(mice) # For Multiple Imputation with Chained Equations in QBA for missing data
library(dplyr) # A fast, consistent tool for working with data frame like objects, both in memory and out of memory
library(tableone) # Creates 'Table 1', i.e., description of baseline patient characteristics
library(WeightIt) # Generates balancing weights for causal effect estimation in observational studies
library(cobalt) # Generates balance tables and plots for covariates of groups preprocessed through matching, weighting or subclassification
library(survey) # Calculates summary statistics, two-sample tests, rank tests, generalised linear models, cumulative link models, Cox models, loglinear models, and general maximum pseudolikelihood estimation for multistage stratified, cluster-sampled, unequally weighted survey samples.
library(survminer) # Contains the function 'ggsurvplot()' for drawing easily beautiful and 'ready-to-publish' survival curves with the 'number at risk' table and 'censoring count plot'.

# Custom functions

# 1. Wrapper for mice() function to suppress "non-integer" warnings
# df: data frame with missing values
# n : integer for number of imputations
# verbose: Logical, if TRUE, prints history in console

impute <- function(df, n=1, verbose=TRUE) {
  withCallingHandlers({
    imp <- mice(df, meth=meth, pred=pred, m=n, print=verbose)
  }, 
  warning=function(w) {
    if (grepl("non-integer", w)) {
      invokeRestart("muffleWarning")
    } # non-integer success warnings are fine,
    # mice uses weighted regression so target is not integer
  })
  return(imp)
}




# 2. Helper function to create cohorts from a data frame by filtering patients receiving either the active treatment or control regimen. It also applies the eligibility criterion of ecog < 2
# x: data frame with columns `linename`, `os`, `event_os` and `ecog`

cohortize <- function(x) {
  x %>%
    rename(drug=linename, OS=os, EVENT=event_os) %>%
    filter(drug %in% c(treated, ref)) %>%
    filter(ecog %in% c(0, 1)) %>%
    droplevels()
}




# 3. Function to run IPTW Cox model
# df: data frame
# covs: a character vector of column names in df that represent baseline covariates for calculating propensity scores for IPTW
# ref: character representing the control treatment line name 
# treated: character representing the active treatment line name 

iptw <- function(df, covs, ref, treated) {
  print(paste0("Reference: ", ref))
  print(paste0("Treated: ", treated))
  
  # Create treatment indicator
  df <- df %>%
    filter(drug %in% c(ref, treated)) %>%
    mutate(drug=ifelse(drug==treated, 1, 0)) %>%
    tidyr::drop_na(all_of(covs))
  
  # Create unadjusted baseline characteristics table with standardized mean differences (SMD)
  tabs <- CreateTableOne(vars=covs, strata="drug", data=df)
  print(tabs, smd=TRUE)
  
  # Compute propensity scores (PS) using pedictions from logistic regression for treatment given the covariates
  W <- weightit(f.build("drug", covs), data=df, method="ps", estimand="ATT", stabilize=FALSE)
  df$wts <- W$ps 
  
  # trim extreme weights (i.e. >3) to the 99th percentile
  if (any(df$wts > 3)) {
    df$wts <- trim(df$wts, at=0.99)
  }
  print(summary(df$wts))
  
  # Compute re-weighted data
  svy <- svydesign(ids=~1, data=df, weights=~df$wts)
  
  # Create baseline characteristics table with SMDs after IPTW using PS
  svyt <- svyCreateTableOne(data = svy, strata = "drug", vars=covs)
  smds <- ExtractSmd(svyt)[,1]
  print(svyt, smd=T)
  if (any(smds>0.12)) {
    message("Imbalance: ", paste(names(smds[smds>0.1]), collapse=","))
    message("Imbalance: ", paste((smds[smds>0.1]) %>% round(3), collapse=","))
  }
  
  # Generate KM curves
  km.fit <- do.call(survfit, list(Surv(OS, EVENT) ~ drug, data = df, weight = df$wts))
  gg <- suppressWarnings(survminer::ggsurvplot(km.fit,
                                               palette=RColorBrewer::brewer.pal(3, "Set1"),
                                               pval=TRUE,
                                               pval.size=4.5,
                                               size=0.7,
                                               censor.size=4,
                                               legend.labs=c(ref.form, treated),
                                               xlab="Time (months)",
                                               risk.table=T,
                                               surv.median.line="hv",
                                               pval.coord = c(45, 0.03),
                                               conf.int=T))
  
  # Compute Cox model
  message("Cox regression")
  cox.fit <- coxph(as.formula(paste0("Surv(OS, EVENT) ~ drug + ", paste0(covs, collapse="+"))),
                   data=df, weights=df$wts, robust=TRUE)
  print(summary(coxph(as.formula(paste0("Surv(OS, EVENT) ~ drug + ", paste0(covs, collapse="+"))),
                      data=df, weights=NULL, robust=TRUE)))
  print(cox.fit)
  list(
    km=gg,
    surv=suppressWarnings(survminer::surv_median(km.fit)),
    n0=nrow(df[df$drug==0,]),
    n1=nrow(df[df$drug==1,]),
    HR=cox.fit$coefficients[["drug"]],
    lower=confint(cox.fit)["drug",][[1]],
    upper=confint(cox.fit)["drug",][[2]],
    seHR=summary(cox.fit)$coefficients[1,4],
    coef=summary(cox.fit),
    meds=surv_median(km.fit))
}




# 4. Function to pool estimates after imputation using Rubin's rules
# iptw.out: list of imputed data frames with IPTW weights
# niters: number of imputation iterations

pool <- function(iptw.out, niters) {
  # Pool point estimates
  ests <- sapply(iptw.out, function(i) i$HR)
  est <- mean(ests)
  
  # Pooled within variance
  se <- sapply(iptw.out, function(i) i$seHR)
  var.wit <- mean(se**2)
  
  # Pooled between variance
  var.bet <- 1/(niters-1) * sum((est - ests)**2)
  
  # Total variance
  var.tot <- var.wit + var.bet + var.bet/niters
  
  # Degrees of freedom
  lambda <- (var.bet + (var.bet/niters)) / var.tot
  df.old <- (niters-1)/(lambda^2)
  n <- iptw.out[[1]]$n1 + iptw.out[[1]]$n0
  df.obs <- ((n-1)+1)/((n-1)+3) * (n-1) * (1-lambda)
  df.adj <- (df.old * df.obs) / (df.old + df.obs)
  
  # 95% confidence intervals
  cint <- qt(c(.025, .975), df=df.adj) * (sqrt(var.tot)) + est
  
  # p-values
  tStat <- est / sqrt(var.tot)
  pval <- 2 * pt(-abs(tStat), df=df.adj)
  
  # Mean median survival time
  medsurv <- sapply(iptw.out, function(i) i$meds$median) %>% rowMeans()
  lowsurv <- sapply(iptw.out, function(i) i$meds$lower) %>% rowMeans()
  uppsurv <- sapply(iptw.out, function(i) i$meds$upper) %>% rowMeans()
  
  # Return estimates
  list(HR=est, lower=cint[1], upper=cint[2], p=pval,
       median_surv=medsurv, lower_surv=lowsurv, upper_surv=uppsurv)
}




# 5. Function for computing E-values
# pt_est: original point estimate (number)
# lower: number indicating the lower boundary of the pt_est's 95% confidence interval
# upper: number indicating the upper boundary of the pt_est's 95% confidence interval
# type: type of pt_est, either "HR", "OR" or "RR".
evalue <- function(pt_est,
                   lower=NA,
                   upper=NA,
                   type=c("HR", "OR", "RR"),
                   common=TRUE) {
  # Estimates
  ests <- c(pt_est, lower, upper)
  
  # E-value computation
  .compute_evalue <- function(x) {
    x <- ifelse(x < 1, 1/x, x)
    x + sqrt(x * (x - 1))
  }
  type=match.arg(type)
  
  # Convert effect measures to approximate risk ratios if needed
  if (type=="RR") {
    rr <- ests
  } else if (type=="OR") {
    rr <- EValue::toRR(EValue::OR(ests, rare=!common))
  } else if (type=="HR") {
    rr <- EValue::toRR(EValue::HR(ests, rare=!common))
  }
  
  # Compute E-value from risk ratios
  e <- .compute_evalue(rr)
  
  # Return results
  data.frame(
    at=c("point", "lower", "upper"),
    value=ests, # Original estimates
    RR=as.vector(rr), # Estimates as risk ratios
    e_value=e # E-values
  )
}

# load data
df <- readRDS("data.rds")

# Define a vector of delta values to calculate delta-adjusted ECOG values for imputation
deltas <- -5:5

# Specify confounders (propensity scores, i.e., probability for each patient to receive active treatment will be conditioned on these baseline covariates)
covs <- c("age65", "sex", "smokhist2", "tsd2", "stage2", "ecog", "race2")

# Define treatment and control (i.e., reference)
treated <- "Pralsetinib"
refs <- c("Carboplatin,Pembrolizumab,Pemetrexed", "Pembrolizumab")
ref <- refs[1] # Either 1 or 2

# Formatted text for graphs etc
ref.form <- ifelse(ref=="Carboplatin,Pembrolizumab,Pemetrexed", "Pembro+chemo", "Pembrolizumab")

# Choose delta value to estimate HRs (i.e. "default" HR)
delta <- 0

#######################
###
### MULTIPLE IMPUTATION
###
#######################

# The following steps impute missing values in the data `df`, including delta-adjusted imputation for the ordinal variable, ECOG.

# Initialize mice
init <- mice(df, maxit = 0)

# Extract predictor matrix
pred <- init$pred

# Set delta adjustment method for ECOG
meth <- init$meth
meth["ecog"] <- "deltaecog"

# Run imputations for a range of delta values and save these imputed files
print("Running imputations...")
for (delta in deltas) {
    message(paste0("On delta: ", delta))
    f <- paste0("imputed_data/delta_", delta, ".rds")
    imp <- impute(df, n=3)
    print(imp$loggedEvents)
    saveRDS(imp, f)
}

#######################
###
### EFFECT ESTIMATION
###
#######################

# Read imputed data
tra <- readRDS("imputed_data/delta_0.rds") # Does not change with delta
eca <- readRDS(paste0("imputed_data/delta_", delta, ".rds")) # ECOG PS imputed by delta only in ECA

# For each imputation, concatenate trial arm and ECA and run IPTW-CoxPH
imp.results <- lapply(1:tra$m, function(j) {
    x0 <- cohortize(mice::complete(tra, j)) %>% filter(drug == ref)
    x1 <- cohortize(mice::complete(eca, j)) %>% filter(drug == treated)
    x <- rbind(x0, x1)
    iptw(x, covs, ref=ref, treated=treated)
})

# Check results
imp.results

# Pool results from imputation
out <- pool(imp.results, length(imp.results))
# Print results
print(paste0("HR: ", round(exp(out$HR), 3), "; 95% CI: ", round(exp(out$lower), 3), "-", round(exp(out$upper), 3)))

#######################
###
### E-VALUES
###
#######################

# An E-value for unmeasured confounding is minimum strength of association, on the risk ratio scale, that unmeasured confounder(s) would need to have with both the treatment and the outcome to fully explain away a specific treatment-outcome association, conditional on the measured covariates:

evalue(out$HR %>% exp(), out$lower %>% exp(), out$upper %>% exp())
