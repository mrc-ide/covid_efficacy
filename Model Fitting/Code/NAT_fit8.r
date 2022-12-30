library(rstudioapi)
library(drjacoby)
library(matrixStats)
library(dplyr)
library(ggplot2)


setwd(dirname(getActiveDocumentContext()$path))
cores <- 8

# Parameters & log like

source("NAT_fit_priors_inits_loglik_simp.r")
source("NAT_fit_main_fn_new.r")

### model version

misc <- list(fit=TRUE,AZ_PD2=FALSE,SingleBoost=FALSE,NewDecay=TRUE,AdditiveBoost=FALSE,AltSevMod=FALSE)

set.seed(1790917)
run_name <- paste0("UKHSA_v6_65+_20220702_AZPD2=",misc$AZ_PD2,"_SB=",misc$SingleBoost,"_NewDecay=", misc$NewDecay,"_AddBst=",misc$AdditiveBoost,"_AltSev=",misc$AltSevMod)
data_file <- "UKHSA_VE_Jun22_65+.csv"
r_fitmodel(run_name, data_file, misc)


set.seed(1790917)
run_name <- paste0("UKHSA_v6_20220702_AZPD2=",misc$AZ_PD2,"_SB=",misc$SingleBoost,"_NewDecay=",misc$NewDecay,"_AddBst=",misc$AdditiveBoost,"_AltSev=",misc$AltSevMod)
data_file <- "UKHSA_VE_Jun22.csv"
r_fitmodel(run_name, data_file, misc)


misc <- list(fit=TRUE,AZ_PD2=FALSE,SingleBoost=FALSE,NewDecay=TRUE,AdditiveBoost=TRUE,AltSevMod=FALSE)

set.seed(1790917)
run_name <- paste0("UKHSA_v6_65+_20220702_AZPD2=",misc$AZ_PD2,"_SB=",misc$SingleBoost,"_NewDecay=",misc$NewDecay,"_AddBst=",misc$AdditiveBoost,"_AltSev=",misc$AltSevMod)
data_file <- "UKHSA_VE_Jun22_65+.csv"
r_fitmodel(run_name, data_file, misc)

set.seed(1790917)
run_name <- paste0("UKHSA_v6_20220702_AZPD2=",misc$AZ_PD2,"_SB=",misc$SingleBoost,"_NewDecay=",misc$NewDecay,"_AddBst=",misc$AdditiveBoost,"_AltSev=",misc$AltSevMod)
data_file <- "UKHSA_VE_Jun22.csv"
r_fitmodel(run_name, data_file, misc)

misc <- list(fit=TRUE,AZ_PD2=FALSE,SingleBoost=FALSE,NewDecay=TRUE,AdditiveBoost=FALSE,AltSevMod=TRUE)

set.seed(1790917)
run_name <- paste0("UKHSA_v6_65+_20220702_AZPD2=",misc$AZ_PD2,"_SB=",misc$SingleBoost,"_NewDecay=",misc$NewDecay,"_AddBst=",misc$AdditiveBoost,"_AltSev=",misc$AltSevMod)
data_file <- "UKHSA_VE_Jun22_65+.csv"
r_fitmodel(run_name, data_file, misc)

set.seed(1790917)
run_name <- paste0("UKHSA_v6_20220702_AZPD2=",misc$AZ_PD2,"_SB=",misc$SingleBoost,"_NewDecay=",misc$NewDecay,"_AddBst=",misc$AdditiveBoost,"_AltSev=",misc$AltSevMod)
data_file <- "UKHSA_VE_Jun22.csv"
r_fitmodel(run_name, data_file, misc)
