
r_fitmodel<- function(run_name, data_file, misc) {
# time offsets by dose
t_offset=c(21,14,14)

# load data

data_all <- read.csv(data_file)
data_all$t_min <- data_all$t_min-t_offset[data_all$dose]
data_all$t_max <- data_all$t_max-t_offset[data_all$dose]
data_all$t <- floor((data_all$t_min + data_all$t_max)/2)
data_all$RR <- 1- data_all$VE/100

data_all$RR[data_all$RR>1]=1 ## can't fit > 1
data_fit <- data_all[!is.na(data_all$VE) & data_all$t_min>=0 & data_all$vaccine_num<4,]
data_fit_label <- data_fit$Label
data_fit <- data_fit[,c("variant","vaccine_num","booster_num","dose","endpoint","t","N","RR","L95","U95")]

data_fit$vacc_lookup = ifelse(data_fit$dose<3,
                              (data_fit$dose-1)*3+data_fit$vaccine_num,
                              6+(data_fit$booster_num-1)*3+data_fit$vaccine_num)
data_fit$N[data_fit$N>3000]=3000

# MCMC --------------------------------------------------------------------

#cores <- parallel::detectCores()

if(misc$fit==TRUE) {

  cl <- parallel::makeCluster(cores)
  
  # run MCMC
  mcmc <- run_mcmc(data = data_fit,
                   df_params = df_params,
                   misc = misc,
                   loglike = r_loglike,
                   logprior = r_logprior,
                   burnin =  50000,
                   samples = 50000,
                   rungs = cores,
                   chains = 8,
                   alpha = 5.0,
                   cluster = cl,
                   pb_markdown = TRUE)
  
  parallel::stopCluster(cl)
  
  save(mcmc,file=paste0(run_name,"_mcmc_chain.Rdata"))
  
} else
  load(file=paste0(run_name,"_mcmc_chain.Rdata"))

#plot_par(mcmc,phase="burnin")
#plot_par(mcmc)
#plot_cor(mcmc, "k", "fold_red_AZ")

print("ESS:")
print(mcmc$diagnostics$ess)
print("Rhat:")
print(mcmc$diagnostics$rhat)
misc$fit=FALSE

#calculate parameter median and 95% CrI
mcmc.samp <- mcmc$output[mcmc$output$phase=="sampling",4:(3+nrow(df_params))]
mcmc.loglik <- mcmc$output$loglikelihood[mcmc$output$phase=="sampling"]
mode.index <- which.max(mcmc.loglik)  ## modal sample
param.est <- as.data.frame(colQuantiles(as.matrix(mcmc.samp),probs=c(0.5,0.025,0.975)))
param.est$variable=rownames(param.est)
param.est$mode=as.numeric(mcmc.samp[mode.index,])
write.csv(param.est,paste0(run_name,"_param_est.csv"))

## calculated predicted values for data
mcmc.samp <- mcmc$output[mcmc$output$phase=="sampling" & mcmc$output$iteration %% 10==0,4:(3+nrow(df_params))]
post.pred <- matrix(ncol=nrow(mcmc.samp),nrow=nrow(data_fit))
for(i in 1:nrow(mcmc.samp)) {
  temp <- r_loglike(mcmc.samp[i,], data_fit, misc)
  post.pred[,i] <- temp[[1]]
}
post.pred.q <- rowQuantiles(post.pred,probs=c(0.5,0.025,0.975))
colnames(post.pred.q) <- c("predicted","l95","u95")  

## and add to data table
data_fit$Label <- data_fit_label
data_fit_pred <- data_fit
data_fit_pred$VE <- 1-data_fit_pred$RR
data_fit_pred <- cbind(data_fit_pred,post.pred.q)
data_fit_pred$dose <- as.factor(data_fit_pred$dose)
variants <- c("Delta","Omicron")
data_fit_pred$variant <- variants[data_fit_pred$variant]
vaccines <- c("AZ","Pfizer","MD")
data_fit_pred$vaccine <- vaccines[data_fit_pred$vaccine_num]
outcomes <- c("Infection","Hospitalisation","Death")
data_fit_pred$outcome <- outcomes[data_fit_pred$endpoint]
data_fit_pred$L95=data_fit_pred$L95/100
data_fit_pred$U95=data_fit_pred$U95/100
data_fit_pred$L95[data_fit_pred$L95<0]=0

write.csv(data_fit_pred,paste0(run_name,"_pred_data.csv"))

data_fit_pred1 <- data_fit_pred[data_fit_pred$variant=="Delta",]

## plot
fit.plot <- data_fit_pred1 %>% ggplot( aes(x=t, y=predicted, ymin=l95, ymax=u95, group=dose, colour=dose,fill=dose)) +
  facet_grid(rows=vars(Label),cols=vars(outcome)) +
  geom_pointrange(position = "jitter") + geom_pointrange(aes(y=VE,ymin=L95,ymax=U95),shape=0,size=0.5) +
  theme_bw() +
  scale_x_continuous(breaks = seq(0, 180, by = 30)) +
  scale_y_continuous(breaks = seq(0.0, 1, by = 0.1)) +
  theme(strip.background = element_blank(),
        panel.border = element_blank())

fit.plot
ggsave(paste0(run_name,"_fit_delta.pdf"),fit.plot,width=8,height=14)


data_fit_pred2 <- data_fit_pred[data_fit_pred$variant=="Omicron",]
## plot
fit.plot2 <- data_fit_pred2 %>% ggplot( aes(x=t, y=predicted, ymin=l95, ymax=u95, group=dose, colour=dose,fill=dose)) +
  facet_grid(rows=vars(Label),cols=vars(outcome)) +
  geom_pointrange(position = "jitter") + geom_pointrange(aes(y=VE,ymin=L95,ymax=U95),shape=0,size=0.5) +
  theme_bw() +
  scale_x_continuous(breaks = seq(0, 180, by = 30)) +
  scale_y_continuous(breaks = seq(0.0, 1, by = 0.1)) +
  theme(strip.background = element_blank(),
        panel.border = element_blank())

fit.plot2
ggsave(paste0(run_name,"_fit_omicron.pdf"),fit.plot2,width=8,height=14)


## now calculate PD1, PD2, PD3 profiles (median and CrI)
data_profile1 <- data.frame(variant=rep(1,730*3*3),
                          vaccine_num=c(rep(1,730*3),rep(2,730*3),rep(3,730*3)),
                          booster_num=c(rep(1,730*3),rep(2,730*3),rep(3,730*3)),
                          dose=rep(c(rep(1,730),rep(2,730),rep(3,730)),3),
                          endpoint=rep(1,730*3*3),
                          t=rep(1:730,3*3))
data_profile1$vacc_lookup = ifelse(data_profile1$dose<3,
                              (data_profile1$dose-1)*3+data_profile1$vaccine_num,
                              6+(data_profile1$booster_num-1)*3+data_profile1$vaccine_num)
data_profile <- data_profile1
data_profile1$endpoint <- 2
data_profile <- rbind(data_profile,data_profile1)
data_profile1$endpoint <- 3
data_profile <- rbind(data_profile,data_profile1)
data_profile1 <- data_profile
data_profile1$variant <-2
data_profile <- rbind(data_profile,data_profile1)

mcmc.samp <- mcmc$output[mcmc$output$phase=="sampling" & mcmc$output$iteration %% 100==0,4:(3+nrow(df_params))]
mcmc.loglik <- mcmc$output$loglikelihood[mcmc$output$phase=="sampling" & mcmc$output$iteration %% 100==0]
mode.index <- which.max(mcmc.loglik)  ## modal sample
post.pred <- matrix(ncol=nrow(mcmc.samp),nrow=nrow(data_profile))
post.drvec <- post.pred
for(i in 1:nrow(mcmc.samp)) {
  temp <- r_loglike(mcmc.samp[i,], data_profile, misc)
  post.pred[,i] <- temp[[1]]
  post.drvec[,i] <- temp[[2]]
}
post.pred.q <- rowQuantiles(post.pred,probs=c(0.5,0.025,0.975))
post.drvec.q <- rowQuantiles(post.drvec,probs=c(0.5,0.025,0.975))
post.pred.mode <-  post.pred[, mode.index]
post.drvec.mode <-  post.drvec[, mode.index]
post.pred.q <- cbind(post.pred.q, post.pred.mode)
post.drvec.q <- cbind(post.drvec.q, post.drvec.mode)
colnames(post.pred.q) <- c("Median.VE","VE.l95","VE.u95","VE.mode")
colnames(post.drvec.q) <- c("Median.drvec","drvec.l95","drvec.u95","drvec.mode")


variants <- c("Delta","Omicron")
data_profile$variant <- variants[data_profile$variant]
vaccines <- c("AZ","Pfizer")
data_profile$vaccine <- vaccines[data_profile$vaccine_num]
outcomes <- c("Infection","Hospitalisation","Death")
data_profile$outcome <- outcomes[data_profile$endpoint]
data_profile <- cbind(data_profile,post.pred.q,post.drvec.q)
write.csv(data_profile,paste0(run_name,"_pred_profile.csv"))

}




