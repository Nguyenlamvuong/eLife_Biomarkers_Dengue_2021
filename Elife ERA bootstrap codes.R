#----------------------------------------------------------------------------------------------------
# title: "Combination of inflammatory and vascular markers in the febrile phase of dengue is associated with more severe outcomes"
# author: "Nguyen Lam Vuong"
# date: "29-Jul-2021"
# SOME FUNCTIONS FOR USE IN THE ANALYSIS

#----------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------
# Load packages
library(tidyverse)
library(MuMIn)
source("Elife ERA functions.R") # to include my functions
options(na.action = "na.fail") # for 'dredge' function [MuMIn]

# Load data
dat <- read_csv("Dengue_Biomarkers_data_27Jul2021.csv") %>% 
  filter(Time == "Enrolment") %>%
  mutate(VCAM = ifelse(uVCAM==1, log2(0.028), log2(VCAM1)),
         SDC = log2(SDC1),
         Ang = ifelse(uAng==1, log2(17.1), log2(Ang2)),
         IL8 = ifelse(uIL8==1, log2(1.8), log2(IL8)),
         IP10 = ifelse(uIP10==1, log2(1.18), log2(IP10)),
         IL1RA = ifelse(uIL1RA==1, log2(18), log2(IL1RA)),
         CD163 = log2(CD163),
         TREM = ifelse(uTREM==1, log2(10.65), log2(TREM1)),
         Fer = log2(Fer),
         CRP = log2(CRP), 
         Vir = log10(Viremia)) %>%
  select(Code, Age, age15, sev.or.inte, VCAM, SDC, Ang, IL8, IP10, IL1RA, CD163, TREM, Fer, CRP)

dat1 <- dat %>% filter(age15 == "No") ## Data for children
dat2 <- dat %>% filter(age15 == "Yes") ## Data for adults

#----------------------------------------------------------------------------------------------------
## The following codes were modified from Heinze G. et al (https://doi.org/10.1002/bimj.201700067)

#----------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------
# FOR CHILDREN

#----------------------------------------------------------------------------------------------------
# Initial setup

pred <- c("VCAM", "SDC", "Ang", "IL8", "IP10", "IL1RA", "CD163", "TREM", "Fer", "CRP")
full_mod <- glm(sev.or.inte ~ VCAM + SDC + Ang + IL8 + IP10 + IL1RA + CD163 + TREM + Fer + CRP, 
                family=binomial, data=dat1, x=T, y=T)
full_est <- coef(full_mod)
pred_name <- names(full_est)[-1]
bootnum <- 1000

set.seed(51086)

#----------------------------------------------------------------------------------------------------
# Bootstrap for best of all combinations for children
boot_var <- matrix(0, ncol = length(pred) + 1, nrow = bootnum, dimnames = list(NULL, c(pred, "aic")))
boot_est <-  boot_se <- matrix(0, ncol = length(pred_name) + 1, nrow = bootnum,
                               dimnames = list(NULL, c("(Intercept)", pred_name)))

for (i in 1:bootnum) {
  data_id <- sample(1:dim(dat1)[1], replace = T)
  dat_id <- dat1[data_id, ]
  boot_m <- glm(sev.or.inte ~ VCAM + SDC + Ang + IL8 + IP10 + IL1RA + CD163 + TREM + Fer + CRP, 
                family=binomial, data=dat_id, x=T, y=T)
  boot_sel <- dredge(boot_m, rank="AIC")
  boot_var_tmp <- attr(get.models(boot_sel, 1)[[1]]$terms, "term.labels")
  
  for (j in 1:(ncol(boot_var)-1)) {
    boot_var[i,j] <- ifelse(names(boot_var[i,j]) %in% boot_var_tmp, 1, 0)
  }
  
  formula <- paste("sev.or.inte~", paste(names(boot_var[i,][boot_var[i,]==1]), collapse = "+"))
  boot_mod <- glm(formula, data = dat_id, family = binomial, x = T, y = T)
  boot_var[i, ncol(boot_var)] <- AIC(boot_mod)
  boot_est[i, names(coef(boot_mod))] <- coef(boot_mod)
  boot_se[i, names(coef(boot_mod))] <- coef(summary(boot_mod))[, "Std. Error"]
}

save(boot_var, boot_est, boot_se, file="5a.RData") # Save bootstrap results for later use

#----------------------------------------------------------------------------------------------------
# Bootstrap for best combination of 2 variables for children

boot_var <- matrix(0, ncol = length(pred) + 1, nrow = bootnum, dimnames = list(NULL, c(pred, "aic")))
boot_est <-  boot_se <- matrix(0, ncol = length(pred_name) + 1, nrow = bootnum,
                               dimnames = list(NULL, c("(Intercept)", pred_name)))

for (i in 1:bootnum) {
  data_id <- sample(1:dim(dat1)[1], replace = T)
  dat_id <- dat1[data_id, ]
  boot_m <- glm(sev.or.inte ~ VCAM + SDC + Ang + IL8 + IP10 + IL1RA + CD163 + TREM + Fer + CRP, 
                family=binomial, data=dat_id, x=T, y=T)
  boot_sel <- dredge(boot_m, rank="AIC", m.lim=c(2,2))
  boot_var_tmp <- attr(get.models(boot_sel, 1)[[1]]$terms, "term.labels")
  
  for (j in 1:(ncol(boot_var)-1)) {
    boot_var[i,j] <- ifelse(names(boot_var[i,j]) %in% boot_var_tmp, 1, 0)
  }
  
  formula <- paste("sev.or.inte~", paste(names(boot_var[i,][boot_var[i,]==1]), collapse = "+"))
  boot_mod <- glm(formula, data = dat_id, family = binomial, x = T, y = T)
  boot_var[i, ncol(boot_var)] <- AIC(boot_mod)
  boot_est[i, names(coef(boot_mod))] <- coef(boot_mod)
  boot_se[i, names(coef(boot_mod))] <- coef(summary(boot_mod))[, "Std. Error"]
}

save(boot_var, boot_est, boot_se, file="5b2.RData") # Save bootstrap results for later use

#----------------------------------------------------------------------------------------------------
# Bootstrap for best combination of 3 variables for children

boot_var <- matrix(0, ncol = length(pred) + 1, nrow = bootnum, dimnames = list(NULL, c(pred, "aic")))
boot_est <-  boot_se <- matrix(0, ncol = length(pred_name) + 1, nrow = bootnum,
                               dimnames = list(NULL, c("(Intercept)", pred_name)))

for (i in 1:bootnum) {
  data_id <- sample(1:dim(dat1)[1], replace = T)
  dat_id <- dat1[data_id, ]
  boot_m <- glm(sev.or.inte ~ VCAM + SDC + Ang + IL8 + IP10 + IL1RA + CD163 + TREM + Fer + CRP, 
                family=binomial, data=dat_id, x=T, y=T)
  boot_sel <- dredge(boot_m, rank="AIC", m.lim=c(3,3))
  boot_var_tmp <- attr(get.models(boot_sel, 1)[[1]]$terms, "term.labels")
  
  for (j in 1:(ncol(boot_var)-1)) {
    boot_var[i,j] <- ifelse(names(boot_var[i,j]) %in% boot_var_tmp, 1, 0)
  }
  
  formula <- paste("sev.or.inte~", paste(names(boot_var[i,][boot_var[i,]==1]), collapse = "+"))
  boot_mod <- glm(formula, data = dat_id, family = binomial, x = T, y = T)
  boot_var[i, ncol(boot_var)] <- AIC(boot_mod)
  boot_est[i, names(coef(boot_mod))] <- coef(boot_mod)
  boot_se[i, names(coef(boot_mod))] <- coef(summary(boot_mod))[, "Std. Error"]
}

save(boot_var, boot_est, boot_se, file="5b3.RData") # Save bootstrap results for later use

#----------------------------------------------------------------------------------------------------
# Bootstrap for best combination of 4 variables for children

boot_var <- matrix(0, ncol = length(pred) + 1, nrow = bootnum, dimnames = list(NULL, c(pred, "aic")))
boot_est <-  boot_se <- matrix(0, ncol = length(pred_name) + 1, nrow = bootnum,
                               dimnames = list(NULL, c("(Intercept)", pred_name)))

for (i in 1:bootnum) {
  data_id <- sample(1:dim(dat1)[1], replace = T)
  dat_id <- dat1[data_id, ]
  boot_m <- glm(sev.or.inte ~ VCAM + SDC + Ang + IL8 + IP10 + IL1RA + CD163 + TREM + Fer + CRP, 
                family=binomial, data=dat_id, x=T, y=T)
  boot_sel <- dredge(boot_m, rank="AIC", m.lim=c(4,4))
  boot_var_tmp <- attr(get.models(boot_sel, 1)[[1]]$terms, "term.labels")
  
  for (j in 1:(ncol(boot_var)-1)) {
    boot_var[i,j] <- ifelse(names(boot_var[i,j]) %in% boot_var_tmp, 1, 0)
  }
  
  formula <- paste("sev.or.inte~", paste(names(boot_var[i,][boot_var[i,]==1]), collapse = "+"))
  boot_mod <- glm(formula, data = dat_id, family = binomial, x = T, y = T)
  boot_var[i, ncol(boot_var)] <- AIC(boot_mod)
  boot_est[i, names(coef(boot_mod))] <- coef(boot_mod)
  boot_se[i, names(coef(boot_mod))] <- coef(summary(boot_mod))[, "Std. Error"]
}

save(boot_var, boot_est, boot_se, file="5b4.RData") # Save bootstrap results for later use

#----------------------------------------------------------------------------------------------------
# Bootstrap for best combination of 5 variables for children

boot_var <- matrix(0, ncol = length(pred) + 1, nrow = bootnum, dimnames = list(NULL, c(pred, "aic")))
boot_est <-  boot_se <- matrix(0, ncol = length(pred_name) + 1, nrow = bootnum,
                               dimnames = list(NULL, c("(Intercept)", pred_name)))

for (i in 1:bootnum) {
  data_id <- sample(1:dim(dat1)[1], replace = T)
  dat_id <- dat1[data_id, ]
  boot_m <- glm(sev.or.inte ~ VCAM + SDC + Ang + IL8 + IP10 + IL1RA + CD163 + TREM + Fer + CRP, 
                family=binomial, data=dat_id, x=T, y=T)
  boot_sel <- dredge(boot_m, rank="AIC", m.lim=c(5,5))
  boot_var_tmp <- attr(get.models(boot_sel, 1)[[1]]$terms, "term.labels")
  
  for (j in 1:(ncol(boot_var)-1)) {
    boot_var[i,j] <- ifelse(names(boot_var[i,j]) %in% boot_var_tmp, 1, 0)
  }
  
  formula <- paste("sev.or.inte~", paste(names(boot_var[i,][boot_var[i,]==1]), collapse = "+"))
  boot_mod <- glm(formula, data = dat_id, family = binomial, x = T, y = T)
  boot_var[i, ncol(boot_var)] <- AIC(boot_mod)
  boot_est[i, names(coef(boot_mod))] <- coef(boot_mod)
  boot_se[i, names(coef(summary(boot_mod)))] <- coef(summary(boot_mod))[, "Std. Error"]
}

save(boot_var, boot_est, boot_se, file="5b5.RData") # Save bootstrap results for later use

#----------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------
# FOR ADULTS

#----------------------------------------------------------------------------------------------------
# Initial setup

pred <- c("VCAM", "SDC", "Ang", "IL8", "ns1(IP10)", "IL1RA", "CD163", "TREM", "Fer", "CRP")
full_mod <- glm(sev.or.inte ~ VCAM + SDC + Ang + IL8 + ns1(IP10) + IL1RA + CD163 + TREM + Fer + CRP, 
                family=binomial, data=dat2, x=T, y=T)
full_est <- coef(full_mod)
pred_name <- names(full_est)[-1]
bootnum <- 1000

set.seed(51086)

#----------------------------------------------------------------------------------------------------
# Bootstrap for best of all combinations for adults

boot_var <- matrix(0, ncol = length(pred) + 1, nrow = bootnum, dimnames = list(NULL, c(pred, "aic")))
boot_est <-  boot_se <- matrix(0, ncol = length(pred_name) + 1, nrow = bootnum,
                               dimnames = list(NULL, c("(Intercept)", pred_name)))

for (i in 1:bootnum) {
  data_id <- sample(1:dim(dat1)[1], replace = T)
  dat_id <- dat1[data_id, ]
  boot_m <- glm(sev.or.inte ~ VCAM + SDC + Ang + IL8 + ns1(IP10) + IL1RA + CD163 + TREM + Fer + CRP, 
                family=binomial, data=dat_id, x=T, y=T)
  boot_sel <- dredge(boot_m, rank="AIC")
  boot_var_tmp <- attr(get.models(boot_sel, 1)[[1]]$terms, "term.labels")
  
  for (j in 1:(ncol(boot_var)-1)) {
    boot_var[i,j] <- ifelse(names(boot_var[i,j]) %in% boot_var_tmp, 1, 0)
  }
  
  formula <- paste("sev.or.inte~", paste(names(boot_var[i,][boot_var[i,]==1]), collapse = "+"))
  boot_mod <- glm(formula, data = dat_id, family = binomial, x = T, y = T)
  boot_var[i, ncol(boot_var)] <- AIC(boot_mod)
  boot_est[i, names(coef(boot_mod))] <- coef(boot_mod)
  boot_se[i, names(coef(boot_mod))] <- coef(summary(boot_mod))[, "Std. Error"]
}

save(boot_var, boot_est, boot_se, file="5c.RData") # Save bootstrap results for later use

#----------------------------------------------------------------------------------------------------
# Bootstrap for best combination of 2 variables for adults

boot_var <- matrix(0, ncol = length(pred) + 1, nrow = bootnum, dimnames = list(NULL, c(pred, "aic")))
boot_est <-  boot_se <- matrix(0, ncol = length(pred_name) + 1, nrow = bootnum,
                               dimnames = list(NULL, c("(Intercept)", pred_name)))

for (i in 1:bootnum) {
  data_id <- sample(1:dim(dat1)[1], replace = T)
  dat_id <- dat1[data_id, ]
  boot_m <- glm(sev.or.inte ~ VCAM + SDC + Ang + IL8 + ns1(IP10) + IL1RA + CD163 + TREM + Fer + CRP, 
                family=binomial, data=dat_id, x=T, y=T)
  boot_sel <- dredge(boot_m, rank="AIC", m.lim=c(2,2))
  boot_var_tmp <- attr(get.models(boot_sel, 1)[[1]]$terms, "term.labels")
  
  for (j in 1:(ncol(boot_var)-1)) {
    boot_var[i,j] <- ifelse(names(boot_var[i,j]) %in% boot_var_tmp, 1, 0)
  }
  
  formula <- paste("sev.or.inte~", paste(names(boot_var[i,][boot_var[i,]==1]), collapse = "+"))
  boot_mod <- glm(formula, data = dat_id, family = binomial, x = T, y = T)
  boot_var[i, ncol(boot_var)] <- AIC(boot_mod)
  boot_est[i, names(coef(boot_mod))] <- coef(boot_mod)
  boot_se[i, names(coef(boot_mod))] <- coef(summary(boot_mod))[, "Std. Error"]
}

save(boot_var, boot_est, boot_se, file="5d2.RData") # Save bootstrap results for later use

#----------------------------------------------------------------------------------------------------
# Bootstrap for best combination of 3 variables for adults

boot_var <- matrix(0, ncol = length(pred) + 1, nrow = bootnum, dimnames = list(NULL, c(pred, "aic")))
boot_est <-  boot_se <- matrix(0, ncol = length(pred_name) + 1, nrow = bootnum,
                               dimnames = list(NULL, c("(Intercept)", pred_name)))

for (i in 1:bootnum) {
  data_id <- sample(1:dim(dat1)[1], replace = T)
  dat_id <- dat1[data_id, ]
  boot_m <- glm(sev.or.inte ~ VCAM + SDC + Ang + IL8 + ns1(IP10) + IL1RA + CD163 + TREM + Fer + CRP, 
                family=binomial, data=dat_id, x=T, y=T)
  boot_sel <- dredge(boot_m, rank="AIC", m.lim=c(3,3))
  boot_var_tmp <- attr(get.models(boot_sel, 1)[[1]]$terms, "term.labels")
  
  for (j in 1:(ncol(boot_var)-1)) {
    boot_var[i,j] <- ifelse(names(boot_var[i,j]) %in% boot_var_tmp, 1, 0)
  }
  
  formula <- paste("sev.or.inte~", paste(names(boot_var[i,][boot_var[i,]==1]), collapse = "+"))
  boot_mod <- glm(formula, data = dat_id, family = binomial, x = T, y = T)
  boot_var[i, ncol(boot_var)] <- AIC(boot_mod)
  boot_est[i, names(coef(boot_mod))] <- coef(boot_mod)
  boot_se[i, names(coef(boot_mod))] <- coef(summary(boot_mod))[, "Std. Error"]
}

save(boot_var, boot_est, boot_se, file="5d3.RData") # Save bootstrap results for later use

#----------------------------------------------------------------------------------------------------
# Bootstrap for best combination of 4 variables for adults

boot_var <- matrix(0, ncol = length(pred) + 1, nrow = bootnum, dimnames = list(NULL, c(pred, "aic")))
boot_est <-  boot_se <- matrix(0, ncol = length(pred_name) + 1, nrow = bootnum,
                               dimnames = list(NULL, c("(Intercept)", pred_name)))

for (i in 1:bootnum) {
  data_id <- sample(1:dim(dat1)[1], replace = T)
  dat_id <- dat1[data_id, ]
  boot_m <- glm(sev.or.inte ~ VCAM + SDC + Ang + IL8 + ns1(IP10) + IL1RA + CD163 + TREM + Fer + CRP, 
                family=binomial, data=dat_id, x=T, y=T)
  boot_sel <- dredge(boot_m, rank="AIC", m.lim=c(4,4))
  boot_var_tmp <- attr(get.models(boot_sel, 1)[[1]]$terms, "term.labels")
  
  for (j in 1:(ncol(boot_var)-1)) {
    boot_var[i,j] <- ifelse(names(boot_var[i,j]) %in% boot_var_tmp, 1, 0)
  }
  
  formula <- paste("sev.or.inte~", paste(names(boot_var[i,][boot_var[i,]==1]), collapse = "+"))
  boot_mod <- glm(formula, data = dat_id, family = binomial, x = T, y = T)
  boot_var[i, ncol(boot_var)] <- AIC(boot_mod)
  boot_est[i, names(coef(boot_mod))] <- coef(boot_mod)
  boot_se[i, names(coef(summary(boot_mod)))] <- coef(summary(boot_mod))[, "Std. Error"]
}

save(boot_var, boot_est, boot_se, file="5d4.RData") # Save bootstrap results for later use

#----------------------------------------------------------------------------------------------------
# Bootstrap for best combination of 5 variables for adults

boot_var <- matrix(0, ncol = length(pred) + 1, nrow = bootnum, dimnames = list(NULL, c(pred, "aic")))
boot_est <-  boot_se <- matrix(0, ncol = length(pred_name) + 1, nrow = bootnum,
                               dimnames = list(NULL, c("(Intercept)", pred_name)))

for (i in 1:bootnum) {
  data_id <- sample(1:dim(dat1)[1], replace = T)
  dat_id <- dat1[data_id, ]
  boot_m <- glm(sev.or.inte ~ VCAM + SDC + Ang + IL8 + ns1(IP10) + IL1RA + CD163 + TREM + Fer + CRP, 
                family=binomial, data=dat_id, x=T, y=T)
  boot_sel <- dredge(boot_m, rank="AIC", m.lim=c(5,5))
  boot_var_tmp <- attr(get.models(boot_sel, 1)[[1]]$terms, "term.labels")
  
  for (j in 1:(ncol(boot_var)-1)) {
    boot_var[i,j] <- ifelse(names(boot_var[i,j]) %in% boot_var_tmp, 1, 0)
  }
  
  formula <- paste("sev.or.inte~", paste(names(boot_var[i,][boot_var[i,]==1]), collapse = "+"))
  boot_mod <- glm(formula, data = dat_id, family = binomial, x = T, y = T)
  boot_var[i, ncol(boot_var)] <- AIC(boot_mod)
  boot_est[i, names(coef(boot_mod))] <- coef(boot_mod)
  boot_se[i, names(coef(summary(boot_mod)))] <- coef(summary(boot_mod))[, "Std. Error"]
}

save(boot_var, boot_est, boot_se, file="5d5.RData") # Save bootstrap results for later use

#----------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------
# SUMMARISE BOOTSTRAP RESULTS FOR CHILDREN

#----------------------------------------------------------------------------------------------------
# Best combination of 2 variables for children

pred <- c("VCAM", "SDC", "Ang", "IL8", "IP10", "IL1RA", "CD163", "TREM", "Fer", "CRP")

## Estimate full model
full_mod <- glm(sev.or.inte ~ VCAM + SDC + Ang + IL8 + IP10 + IL1RA + CD163 + TREM + Fer + CRP, 
                family=binomial, data=dat1, x=T, y=T)
full_est <- coef(full_mod)
full_se <- coef(summary(full_mod))[, "Std. Error"]
pred_name <- names(full_est)[-1]

## Selected model
sel_m <- dredge(full_mod, rank="AIC", m.lim=c(2,2))
sel_var <- matrix(0, ncol = length(pred), nrow = 1, dimnames = list(NULL, pred))
sel_var_tmp <- attr(get.models(sel_m, 1)[[1]]$terms, "term.labels")

for (j in 1:ncol(sel_var)) {
  sel_var[1,j] <- ifelse(names(sel_var[1,j]) %in% sel_var_tmp, 1, 0)
}

formula <- paste("sev.or.inte~", paste(names(sel_var[1,][sel_var[1,]==1]), collapse = "+"))
sel_mod <- glm(formula, data = dat1, family = binomial, x = T, y = T)
sel_est <- coef(sel_mod[1])[c("(Intercept)", pred_name)]
sel_est[is.na(sel_est)] <- 0
names(sel_est) <- c("(Intercept)", pred_name)
sel_se <- coef(summary(sel_mod))[, "Std. Error"][c("(Intercept)", pred_name)]
sel_se[is.na(sel_se)] <- 0

## Bootstrap
bootnum <- 1000
load("5b2.RData") # Load bootstrap results
boot_01 <- (boot_est != 0) * 1
boot_inclusion <- apply(boot_01, 2, function(x) sum(x) / length(x) * 100)

## Overview of estimates and measures
sqe <- (t(boot_est) - full_est) ^ 2
rmsd <- apply(sqe, 1, function(x) sqrt(mean(x)))
rmsdratio <- rmsd / full_se
boot_mean <- apply(boot_est, 2, mean)
boot_meanratio <- boot_mean / full_est
boot_relbias <- (boot_meanratio / (boot_inclusion / 100) - 1) * 100
boot_median <- apply(boot_est, 2, median)
boot_025per <- apply(boot_est, 2, function(x) quantile(x, 0.025))
boot_975per <- apply(boot_est, 2, function(x) quantile(x, 0.975))

overview <- round(cbind(full_est, full_se, boot_inclusion, sel_est, sel_se,  
                        rmsdratio, boot_relbias, boot_median, boot_025per, 
                        boot_975per), 4)
overview <- overview[order(overview[,"boot_inclusion"], decreasing=T),]
overview

## Model frequency
group_cols <- c(colnames(data.frame(boot_var))[1:length(pred)])
boot_modfreq <- data.frame(boot_var) %>%
  mutate(time = 1) %>%
  group_by_(.dots = group_cols) %>%
  summarise(aic_median = round(median(aic),1),
            aic_025 = round(quantile(aic, 0.025),1),
            aic_975 = round(quantile(aic, 0.975),1),
            count = sum(time)) %>%
  ungroup() %>%
  mutate(percent = count/bootnum * 100) %>%
  arrange(desc(count)) %>%
  slice(1:20) %>%
  as.data.frame(.)

out <- boot_modfreq
for (i in 1:(ncol(out)-5)) {
  out[,i] <- ifelse(out[,i]==0, NA, "+")
}

out

#----------------------------------------------------------------------------------------------------
# Best combination of 3 variables for children

pred <- c("VCAM", "SDC", "Ang", "IL8", "IP10", "IL1RA", "CD163", "TREM", "Fer", "CRP")

## Estimate full model
full_mod <- glm(sev.or.inte ~ VCAM + SDC + Ang + IL8 + IP10 + IL1RA + CD163 + TREM + Fer + CRP, 
                family=binomial, data=dat1, x=T, y=T)
full_est <- coef(full_mod)
full_se <- coef(summary(full_mod))[, "Std. Error"]
pred_name <- names(full_est)[-1]

## Selected model
sel_m <- dredge(full_mod, rank="AIC", m.lim=c(3,3))
sel_var <- matrix(0, ncol = length(pred), nrow = 1, dimnames = list(NULL, pred))
sel_var_tmp <- attr(get.models(sel_m, 1)[[1]]$terms, "term.labels")

for (j in 1:ncol(sel_var)) {
  sel_var[1,j] <- ifelse(names(sel_var[1,j]) %in% sel_var_tmp, 1, 0)
}

formula <- paste("sev.or.inte~", paste(names(sel_var[1,][sel_var[1,]==1]), collapse = "+"))
sel_mod <- glm(formula, data = dat1, family = binomial, x = T, y = T)
sel_est <- coef(sel_mod[1])[c("(Intercept)", pred_name)]
sel_est[is.na(sel_est)] <- 0
names(sel_est) <- c("(Intercept)", pred_name)
sel_se <- coef(summary(sel_mod))[, "Std. Error"][c("(Intercept)", pred_name)]
sel_se[is.na(sel_se)] <- 0

## Bootstrap
bootnum <- 1000
load("5b3.RData") # Load bootstrap results
boot_01 <- (boot_est != 0) * 1
boot_inclusion <- apply(boot_01, 2, function(x) sum(x) / length(x) * 100)

## Overview of estimates and measures
sqe <- (t(boot_est) - full_est) ^ 2
rmsd <- apply(sqe, 1, function(x) sqrt(mean(x)))
rmsdratio <- rmsd / full_se
boot_mean <- apply(boot_est, 2, mean)
boot_meanratio <- boot_mean / full_est
boot_relbias <- (boot_meanratio / (boot_inclusion / 100) - 1) * 100
boot_median <- apply(boot_est, 2, median)
boot_025per <- apply(boot_est, 2, function(x) quantile(x, 0.025))
boot_975per <- apply(boot_est, 2, function(x) quantile(x, 0.975))

overview <- round(cbind(full_est, full_se, boot_inclusion, sel_est, sel_se,  
                        rmsdratio, boot_relbias, boot_median, boot_025per, 
                        boot_975per), 4)
overview <- overview[order(overview[,"boot_inclusion"], decreasing=T),]
overview

## Model frequency
group_cols <- c(colnames(data.frame(boot_var))[1:length(pred)])
boot_modfreq <- data.frame(boot_var) %>%
  mutate(time = 1) %>%
  group_by_(.dots = group_cols) %>%
  summarise(aic_median = round(median(aic),1),
            aic_025 = round(quantile(aic, 0.025),1),
            aic_975 = round(quantile(aic, 0.975),1),
            count = sum(time)) %>%
  ungroup() %>%
  mutate(percent = count/bootnum * 100) %>%
  arrange(desc(count)) %>%
  slice(1:20) %>%
  as.data.frame(.)

out <- boot_modfreq
for (i in 1:(ncol(out)-5)) {
  out[,i] <- ifelse(out[,i]==0, NA, "+")
}

out

#----------------------------------------------------------------------------------------------------
# Best combination of 4 variables for children

pred <- c("VCAM", "SDC", "Ang", "IL8", "IP10", "IL1RA", "CD163", "TREM", "Fer", "CRP")

## Estimate full model
full_mod <- glm(sev.or.inte ~ VCAM + SDC + Ang + IL8 + IP10 + IL1RA + CD163 + TREM + Fer + CRP, 
                family=binomial, data=dat1, x=T, y=T)
full_est <- coef(full_mod)
full_se <- coef(summary(full_mod))[, "Std. Error"]
pred_name <- names(full_est)[-1]

## Selected model
sel_m <- dredge(full_mod, rank="AIC", m.lim=c(4,4))
sel_var <- matrix(0, ncol = length(pred), nrow = 1, dimnames = list(NULL, pred))
sel_var_tmp <- attr(get.models(sel_m, 1)[[1]]$terms, "term.labels")

for (j in 1:ncol(sel_var)) {
  sel_var[1,j] <- ifelse(names(sel_var[1,j]) %in% sel_var_tmp, 1, 0)
}

formula <- paste("sev.or.inte~", paste(names(sel_var[1,][sel_var[1,]==1]), collapse = "+"))
sel_mod <- glm(formula, data = dat1, family = binomial, x = T, y = T)
sel_est <- coef(sel_mod[1])[c("(Intercept)", pred_name)]
sel_est[is.na(sel_est)] <- 0
names(sel_est) <- c("(Intercept)", pred_name)
sel_se <- coef(summary(sel_mod))[, "Std. Error"][c("(Intercept)", pred_name)]
sel_se[is.na(sel_se)] <- 0

## Bootstrap
bootnum <- 1000
load("5b4.RData") # Load bootstrap results
boot_01 <- (boot_est != 0) * 1
boot_inclusion <- apply(boot_01, 2, function(x) sum(x) / length(x) * 100)

## Overview of estimates and measures
sqe <- (t(boot_est) - full_est) ^ 2
rmsd <- apply(sqe, 1, function(x) sqrt(mean(x)))
rmsdratio <- rmsd / full_se
boot_mean <- apply(boot_est, 2, mean)
boot_meanratio <- boot_mean / full_est
boot_relbias <- (boot_meanratio / (boot_inclusion / 100) - 1) * 100
boot_median <- apply(boot_est, 2, median)
boot_025per <- apply(boot_est, 2, function(x) quantile(x, 0.025))
boot_975per <- apply(boot_est, 2, function(x) quantile(x, 0.975))

overview <- round(cbind(full_est, full_se, boot_inclusion, sel_est, sel_se,  
                        rmsdratio, boot_relbias, boot_median, boot_025per, 
                        boot_975per), 4)
overview <- overview[order(overview[,"boot_inclusion"], decreasing=T),]
overview

## Model frequency
group_cols <- c(colnames(data.frame(boot_var))[1:length(pred)])
boot_modfreq <- data.frame(boot_var) %>%
  mutate(time = 1) %>%
  group_by_(.dots = group_cols) %>%
  summarise(aic_median = round(median(aic),1),
            aic_025 = round(quantile(aic, 0.025),1),
            aic_975 = round(quantile(aic, 0.975),1),
            count = sum(time)) %>%
  ungroup() %>%
  mutate(percent = count/bootnum * 100) %>%
  arrange(desc(count)) %>%
  slice(1:20) %>%
  as.data.frame(.)

out <- boot_modfreq
for (i in 1:(ncol(out)-5)) {
  out[,i] <- ifelse(out[,i]==0, NA, "+")
}

out

#----------------------------------------------------------------------------------------------------
# Best combination of 5 variables for children

pred <- c("VCAM", "SDC", "Ang", "IL8", "IP10", "IL1RA", "CD163", "TREM", "Fer", "CRP")

## Estimate full model
full_mod <- glm(sev.or.inte ~ VCAM + SDC + Ang + IL8 + IP10 + IL1RA + CD163 + TREM + Fer + CRP, 
                family=binomial, data=dat1, x=T, y=T)
full_est <- coef(full_mod)
full_se <- coef(summary(full_mod))[, "Std. Error"]
pred_name <- names(full_est)[-1]

## Selected model
sel_m <- dredge(full_mod, rank="AIC", m.lim=c(5,5))
sel_var <- matrix(0, ncol = length(pred), nrow = 1, dimnames = list(NULL, pred))
sel_var_tmp <- attr(get.models(sel_m, 1)[[1]]$terms, "term.labels")

for (j in 1:ncol(sel_var)) {
  sel_var[1,j] <- ifelse(names(sel_var[1,j]) %in% sel_var_tmp, 1, 0)
}

formula <- paste("sev.or.inte~", paste(names(sel_var[1,][sel_var[1,]==1]), collapse = "+"))
sel_mod <- glm(formula, data = dat1, family = binomial, x = T, y = T)
sel_est <- coef(sel_mod[1])[c("(Intercept)", pred_name)]
sel_est[is.na(sel_est)] <- 0
names(sel_est) <- c("(Intercept)", pred_name)
sel_se <- coef(summary(sel_mod))[, "Std. Error"][c("(Intercept)", pred_name)]
sel_se[is.na(sel_se)] <- 0

## Bootstrap
bootnum <- 1000
load("5b5.RData") # Load bootstrap results
boot_01 <- (boot_est != 0) * 1
boot_inclusion <- apply(boot_01, 2, function(x) sum(x) / length(x) * 100)

## Overview of estimates and measures
sqe <- (t(boot_est) - full_est) ^ 2
rmsd <- apply(sqe, 1, function(x) sqrt(mean(x)))
rmsdratio <- rmsd / full_se
boot_mean <- apply(boot_est, 2, mean)
boot_meanratio <- boot_mean / full_est
boot_relbias <- (boot_meanratio / (boot_inclusion / 100) - 1) * 100
boot_median <- apply(boot_est, 2, median)
boot_025per <- apply(boot_est, 2, function(x) quantile(x, 0.025, na.rm=T))
boot_975per <- apply(boot_est, 2, function(x) quantile(x, 0.975, na.rm=T))

overview <- round(cbind(full_est, full_se, boot_inclusion, sel_est, sel_se,  
                        rmsdratio, boot_relbias, boot_median, boot_025per, 
                        boot_975per), 4)
overview <- overview[order(overview[,"boot_inclusion"], decreasing=T),]
overview

## Model frequency
group_cols <- c(colnames(data.frame(boot_var))[1:length(pred)])
boot_modfreq <- data.frame(boot_var) %>%
  mutate(time = 1) %>%
  group_by_(.dots = group_cols) %>%
  summarise(aic_median = round(median(aic),1),
            aic_025 = round(quantile(aic, 0.025),1),
            aic_975 = round(quantile(aic, 0.975),1),
            count = sum(time)) %>%
  ungroup() %>%
  mutate(percent = count/bootnum * 100) %>%
  arrange(desc(count)) %>%
  slice(1:20) %>%
  as.data.frame(.)

out <- boot_modfreq
for (i in 1:(ncol(out)-5)) {
  out[,i] <- ifelse(out[,i]==0, NA, "+")
}

out


#----------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------
# SUMMARISE BOOTSTRAP RESULTS FOR ADULTS

#----------------------------------------------------------------------------------------------------
# Best combination of 2 variables for adults

pred <- c("VCAM", "SDC", "Ang", "IL8", "ns1(IP10)", "IL1RA", "CD163", "TREM", "Fer", "CRP")

## Estimate full model
full_mod <- glm(sev.or.inte ~ VCAM + SDC + Ang + IL8 + ns1(IP10) + IL1RA + CD163 + TREM + Fer + CRP, 
                family=binomial, data=dat2, x=T, y=T)
full_est <- coef(full_mod)
full_se <- coef(summary(full_mod))[, "Std. Error"]
pred_name <- names(full_est)[-1]

## Selected model
sel_m <- dredge(full_mod, rank="AIC", m.lim=c(2,2))
sel_var <- matrix(0, ncol = length(pred), nrow = 1, dimnames = list(NULL, pred))
sel_var_tmp <- attr(get.models(sel_m, 1)[[1]]$terms, "term.labels")

for (j in 1:ncol(sel_var)) {
  sel_var[1,j] <- ifelse(names(sel_var[1,j]) %in% sel_var_tmp, 1, 0)
}

formula <- paste("sev.or.inte~", paste(names(sel_var[1,][sel_var[1,]==1]), collapse = "+"))
sel_mod <- glm(formula, data = dat2, family = binomial, x = T, y = T)
sel_est <- coef(sel_mod[1])[c("(Intercept)", pred_name)]
sel_est[is.na(sel_est)] <- 0
names(sel_est) <- c("(Intercept)", pred_name)
sel_se <- coef(summary(sel_mod))[, "Std. Error"][c("(Intercept)", pred_name)]
sel_se[is.na(sel_se)] <- 0

## Bootstrap
bootnum <- 1000
load("5d2.RData")
boot_01 <- (boot_est != 0) * 1
boot_inclusion <- apply(boot_01, 2, function(x) sum(x) / length(x) * 100)

## Overview of estimates and measures
sqe <- (t(boot_est) - full_est) ^ 2
rmsd <- apply(sqe, 1, function(x) sqrt(mean(x)))
rmsdratio <- rmsd / full_se
boot_mean <- apply(boot_est, 2, mean)
boot_meanratio <- boot_mean / full_est
boot_relbias <- (boot_meanratio / (boot_inclusion / 100) - 1) * 100
boot_median <- apply(boot_est, 2, median)
boot_025per <- apply(boot_est, 2, function(x) quantile(x, 0.025))
boot_975per <- apply(boot_est, 2, function(x) quantile(x, 0.975))

overview <- round(cbind(full_est, full_se, boot_inclusion, sel_est, sel_se,  
                        rmsdratio, boot_relbias, boot_median, boot_025per, 
                        boot_975per), 4)
overview <- overview[order(overview[,"boot_inclusion"], decreasing=T),]
overview

## Model frequency
group_cols <- c(colnames(data.frame(boot_var))[1:length(pred)])
boot_modfreq <- data.frame(boot_var) %>%
  mutate(time = 1) %>%
  group_by_(.dots = group_cols) %>%
  summarise(aic_median = round(median(aic),1),
            aic_025 = round(quantile(aic, 0.025),1),
            aic_975 = round(quantile(aic, 0.975),1),
            count = sum(time)) %>%
  ungroup() %>%
  mutate(percent = count/bootnum * 100) %>%
  arrange(desc(count)) %>%
  slice(1:20) %>%
  as.data.frame(.)

out <- boot_modfreq
for (i in 1:(ncol(out)-5)) {
  out[,i] <- ifelse(out[,i]==0, NA, "+")
}

out

#----------------------------------------------------------------------------------------------------
# Best combination of 3 variables for adults

pred <- c("VCAM", "SDC", "Ang", "IL8", "ns1(IP10)", "IL1RA", "CD163", "TREM", "Fer", "CRP")

## Estimate full model
full_mod <- glm(sev.or.inte ~ VCAM + SDC + Ang + IL8 + ns1(IP10) + IL1RA + CD163 + TREM + Fer + CRP, 
                family=binomial, data=dat2, x=T, y=T)
full_est <- coef(full_mod)
full_se <- coef(summary(full_mod))[, "Std. Error"]
pred_name <- names(full_est)[-1]

## Selected model
sel_m <- dredge(full_mod, rank="AIC", m.lim=c(3,3))
sel_var <- matrix(0, ncol = length(pred), nrow = 1, dimnames = list(NULL, pred))
sel_var_tmp <- attr(get.models(sel_m, 1)[[1]]$terms, "term.labels")

for (j in 1:ncol(sel_var)) {
  sel_var[1,j] <- ifelse(names(sel_var[1,j]) %in% sel_var_tmp, 1, 0)
}

formula <- paste("sev.or.inte~", paste(names(sel_var[1,][sel_var[1,]==1]), collapse = "+"))
sel_mod <- glm(formula, data = dat2, family = binomial, x = T, y = T)
sel_est <- coef(sel_mod[1])[c("(Intercept)", pred_name)]
sel_est[is.na(sel_est)] <- 0
names(sel_est) <- c("(Intercept)", pred_name)
sel_se <- coef(summary(sel_mod))[, "Std. Error"][c("(Intercept)", pred_name)]
sel_se[is.na(sel_se)] <- 0

## Bootstrap
bootnum <- 1000
load("5d3.RData") # Load bootstrap results
boot_01 <- (boot_est != 0) * 1
boot_inclusion <- apply(boot_01, 2, function(x) sum(x) / length(x) * 100)

## Overview of estimates and measures
sqe <- (t(boot_est) - full_est) ^ 2
rmsd <- apply(sqe, 1, function(x) sqrt(mean(x)))
rmsdratio <- rmsd / full_se
boot_mean <- apply(boot_est, 2, mean)
boot_meanratio <- boot_mean / full_est
boot_relbias <- (boot_meanratio / (boot_inclusion / 100) - 1) * 100
boot_median <- apply(boot_est, 2, median)
boot_025per <- apply(boot_est, 2, function(x) quantile(x, 0.025))
boot_975per <- apply(boot_est, 2, function(x) quantile(x, 0.975))

overview <- round(cbind(full_est, full_se, boot_inclusion, sel_est, sel_se,  
                        rmsdratio, boot_relbias, boot_median, boot_025per, 
                        boot_975per), 4)
overview <- overview[order(overview[,"boot_inclusion"], decreasing=T),]
overview

## Model frequency
group_cols <- c(colnames(data.frame(boot_var))[1:length(pred)])
boot_modfreq <- data.frame(boot_var) %>%
  mutate(time = 1) %>%
  group_by_(.dots = group_cols) %>%
  summarise(aic_median = round(median(aic),1),
            aic_025 = round(quantile(aic, 0.025),1),
            aic_975 = round(quantile(aic, 0.975),1),
            count = sum(time)) %>%
  ungroup() %>%
  mutate(percent = count/bootnum * 100) %>%
  arrange(desc(count)) %>%
  slice(1:20) %>%
  as.data.frame(.)

out <- boot_modfreq
for (i in 1:(ncol(out)-5)) {
  out[,i] <- ifelse(out[,i]==0, NA, "+")
}

out

#----------------------------------------------------------------------------------------------------
# Best combination of 4 variables for adults

pred <- c("VCAM", "SDC", "Ang", "IL8", "ns1(IP10)", "IL1RA", "CD163", "TREM", "Fer", "CRP")

# Estimate full model ----------------------------------------------
full_mod <- glm(sev.or.inte ~ VCAM + SDC + Ang + IL8 + ns1(IP10) + IL1RA + CD163 + TREM + Fer + CRP, 
                family=binomial, data=dat2, x=T, y=T)
full_est <- coef(full_mod)
full_se <- coef(summary(full_mod))[, "Std. Error"]
pred_name <- names(full_est)[-1]

# Selected model ---------------------------------------------------
sel_m <- dredge(full_mod, rank="AIC", m.lim=c(4,4))
sel_var <- matrix(0, ncol = length(pred), nrow = 1, dimnames = list(NULL, pred))
sel_var_tmp <- attr(get.models(sel_m, 1)[[1]]$terms, "term.labels")

for (j in 1:ncol(sel_var)) {
  sel_var[1,j] <- ifelse(names(sel_var[1,j]) %in% sel_var_tmp, 1, 0)
}

formula <- paste("sev.or.inte~", paste(names(sel_var[1,][sel_var[1,]==1]), collapse = "+"))
sel_mod <- glm(formula, data = dat2, family = binomial, x = T, y = T)
sel_est <- coef(sel_mod[1])[c("(Intercept)", pred_name)]
sel_est[is.na(sel_est)] <- 0
names(sel_est) <- c("(Intercept)", pred_name)
sel_se <- coef(summary(sel_mod))[, "Std. Error"][c("(Intercept)", pred_name)]
sel_se[is.na(sel_se)] <- 0

# Bootstrap ----------------------------------------------------------
bootnum <- 1000
load("5d4.RData") # Load bootstrap results
boot_01 <- (boot_est != 0) * 1
boot_inclusion <- apply(boot_01, 2, function(x) sum(x) / length(x) * 100)

## Overview of estimates and measures
sqe <- (t(boot_est) - full_est) ^ 2
rmsd <- apply(sqe, 1, function(x) sqrt(mean(x)))
rmsdratio <- rmsd / full_se
boot_mean <- apply(boot_est, 2, mean)
boot_meanratio <- boot_mean / full_est
boot_relbias <- (boot_meanratio / (boot_inclusion / 100) - 1) * 100
boot_median <- apply(boot_est, 2, median)
boot_025per <- apply(boot_est, 2, function(x) quantile(x, 0.025, na.rm=T))
boot_975per <- apply(boot_est, 2, function(x) quantile(x, 0.975, na.rm=T))

overview <- round(cbind(full_est, full_se, boot_inclusion, sel_est, sel_se,  
                        rmsdratio, boot_relbias, boot_median, boot_025per, 
                        boot_975per), 4)
overview <- overview[order(overview[,"boot_inclusion"], decreasing=T),]
overview

## Model frequency
group_cols <- c(colnames(data.frame(boot_var))[1:length(pred)])
boot_modfreq <- data.frame(boot_var) %>%
  mutate(time = 1) %>%
  group_by_(.dots = group_cols) %>%
  summarise(aic_median = round(median(aic),1),
            aic_025 = round(quantile(aic, 0.025),1),
            aic_975 = round(quantile(aic, 0.975),1),
            count = sum(time)) %>%
  ungroup() %>%
  mutate(percent = count/bootnum * 100) %>%
  arrange(desc(count)) %>%
  slice(1:20) %>%
  as.data.frame(.)

out <- boot_modfreq
for (i in 1:(ncol(out)-5)) {
  out[,i] <- ifelse(out[,i]==0, NA, "+")
}

out

#----------------------------------------------------------------------------------------------------
# Best combination of 5 variables for adults

pred <- c("VCAM", "SDC", "Ang", "IL8", "ns1(IP10)", "IL1RA", "CD163", "TREM", "Fer", "CRP")

## Estimate full model
full_mod <- glm(sev.or.inte ~ VCAM + SDC + Ang + IL8 + ns1(IP10) + IL1RA + CD163 + TREM + Fer + CRP, 
                family=binomial, data=dat2, x=T, y=T)
full_est <- coef(full_mod)
full_se <- coef(summary(full_mod))[, "Std. Error"]
pred_name <- names(full_est)[-1]

## Selected model
sel_m <- dredge(full_mod, rank="AIC", m.lim=c(5,5))
sel_var <- matrix(0, ncol = length(pred), nrow = 1, dimnames = list(NULL, pred))
sel_var_tmp <- attr(get.models(sel_m, 1)[[1]]$terms, "term.labels")

for (j in 1:ncol(sel_var)) {
  sel_var[1,j] <- ifelse(names(sel_var[1,j]) %in% sel_var_tmp, 1, 0)
}

formula <- paste("sev.or.inte~", paste(names(sel_var[1,][sel_var[1,]==1]), collapse = "+"))
sel_mod <- glm(formula, data = dat2, family = binomial, x = T, y = T)
sel_est <- coef(sel_mod[1])[c("(Intercept)", pred_name)]
sel_est[is.na(sel_est)] <- 0
names(sel_est) <- c("(Intercept)", pred_name)
sel_se <- coef(summary(sel_mod))[, "Std. Error"][c("(Intercept)", pred_name)]
sel_se[is.na(sel_se)] <- 0

## Bootstrap
bootnum <- 1000
load("5d5.RData") # Load bootstrap results
boot_01 <- (boot_est != 0) * 1
boot_inclusion <- apply(boot_01, 2, function(x) sum(x) / length(x) * 100)

## Overview of estimates and measures
sqe <- (t(boot_est) - full_est) ^ 2
rmsd <- apply(sqe, 1, function(x) sqrt(mean(x)))
rmsdratio <- rmsd / full_se
boot_mean <- apply(boot_est, 2, mean)
boot_meanratio <- boot_mean / full_est
boot_relbias <- (boot_meanratio / (boot_inclusion / 100) - 1) * 100
boot_median <- apply(boot_est, 2, median)
boot_025per <- apply(boot_est, 2, function(x) quantile(x, 0.025, na.rm=T))
boot_975per <- apply(boot_est, 2, function(x) quantile(x, 0.975, na.rm=T))

overview <- round(cbind(full_est, full_se, boot_inclusion, sel_est, sel_se,  
                        rmsdratio, boot_relbias, boot_median, boot_025per, 
                        boot_975per), 4)
overview <- overview[order(overview[,"boot_inclusion"], decreasing=T),]
overview

## Model frequency
group_cols <- c(colnames(data.frame(boot_var))[1:length(pred)])
boot_modfreq <- data.frame(boot_var) %>%
  mutate(time = 1) %>%
  group_by_(.dots = group_cols) %>%
  summarise(aic_median = round(median(aic),1),
            aic_025 = round(quantile(aic, 0.025),1),
            aic_975 = round(quantile(aic, 0.975),1),
            count = sum(time)) %>%
  ungroup() %>%
  mutate(percent = count/bootnum * 100) %>%
  arrange(desc(count)) %>%
  slice(1:20) %>%
  as.data.frame(.)

out <- boot_modfreq
for (i in 1:(ncol(out)-5)) {
  out[,i] <- ifelse(out[,i]==0, NA, "+")
}

out
