## function to create titre and efficacy profiles



draw <-
  function(mu_ab_d1 = mu_ab_d1,
           mu_ab_d2 = mu_ab_d2,
           mu_ab_d3,
           std10 = std10,
           ab_50 = ab_50,
           ab_50_severe = ab_50_severe,
           ab_50_death = ab_50_death,
           dr_s = dr_s,
           dr_l = dr_l,
           period_s = period_s,
           t = t,
           t_d2 = t_d2,
           t_d3,
           k = k
  ){
  
  # decay rates over time
  ## simple biphasic decay on natural log scale
  
  denom=log(exp(dr_l*period_s)+exp(dr_s*period_s))
  cum_dr_vec=log(exp(dr_s*t+dr_l*period_s)+exp(dr_l*t+dr_s*period_s))-denom
  dr_vec=c(0,diff(cum_dr_vec,1)) 
 
  
 # dr_vec_sub <- c(rep(dr_s, round(period_s)-1),
 #                 seq(dr_s, dr_l, length.out = round(t_period_l)-round(period_s)+1),
 #                 rep(dr_l, length(t) - round(t_period_l)))
  
 # dr_vec <- rep(0, length(t))
 # dr_vec[2:t_d2] <- dr_vec_sub[1:(t_d2-1)]
 # dr_vec[(2+t_d2):(t_d3+t_d2)] <- dr_vec_sub[1:((t_d3+t_d2)-(t_d2+1))]
 # dr_vec[(2+t_d2+t_d3):length(t)] <- dr_vec_sub[1:(length(t) - (t_d2+t_d3+1)) ]
  
  # nab titre distribution for a given vaccine product: draw from a log-normal distribution
  # here the first titre is dependent on the second but I suppose could do it the other way around
  
#  z2 <- rnorm(1, log10(mu_ab_d2), std10)
 # z3 <- log10(10^z2 * (mu_ab_d3/mu_ab_d2))
#  z1 <- log10(10^z2 * (mu_ab_d1/mu_ab_d2))
  
  z2 <- log10(mu_ab_d2)
  z3 <- log10(mu_ab_d3)
  z1 <- log10(mu_ab_d1)
  
  
  
  # initiate titre vector
  nt <- rep(0, length(t))
  nt[1] <- log(10^z1)
  nt[t_d2+1] <- log(10^z2)
  nt[t_d2+t_d3+1] <- log(10^z3)
  
  # decay antibodies over time on natural log scale
  for (i in (2:t_d2)){
    nt[i] <- nt[i-1] + dr_vec[i]
  }
  for (i in ((t_d2+2):(t_d2+t_d3))){
    nt[i] <- nt[i-1] + dr_vec[i-t_d2-1]
  }
  for (i in ((t_d2+t_d3+2):length(t))){
    nt[i] <- nt[i-1] + dr_vec[i-(t_d2+t_d3)-1]
  }
  
  nt <- exp(nt) # return to linear scale
  
  # relate titre to efficacy over time - using log-10 parameters
  ef_infection <- 1 / (1 + exp(-k * (log10(nt) - log10(ab_50))))
  ef_severe <- 1 / (1 + exp(-k * (log10(nt) - log10(ab_50_severe))))
  ef_death <- 1 / (1 + exp(-k * (log10(nt) - log10(ab_50_death))))
  
  # output list containing titre and vaccine efficacy
  ret <- list(titre = nt, VE = ef_infection, VE_severe = ef_severe, VE_death = ef_death)
  return(ret)
  }
