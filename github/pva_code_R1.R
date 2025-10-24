
# Viability of an endangered songbird population is most affected by adult survival and not brood parasitism
# SG English, A Khan, AM Bezener, M Bieber, TR Forrester, T Luszcz, K Mancuso, R McKibbin, CA Bishop
### LIBRARIES ###############################################################################################

library(nimble)

### LOAD DATA ###############################################################################################

load("./pva_input.rda")

### MODEL SETUP #############################################################################################

dzip <- nimbleFunction(
  run = function(x = integer(), lambda = double(), zeroProb = double(), log = logical(0, default = 0)) {
    returnType(double())
    ## First handle non-zero data
    if (x != 0) {
      ## return the log probability if log = TRUE
      if (log) return(dpois(x, lambda, log = TRUE) + log(1 - zeroProb))
      ## or the probability if log = FALSE
      else return((1 - zeroProb) * dpois(x, lambda, log = FALSE))
    }
    ## From here down we know x is 0
    totalProbZero <- zeroProb + (1 - zeroProb) * dpois(0, lambda, log = FALSE)
    if (log) return(log(totalProbZero))
    return(totalProbZero)
  })

rzip <- nimbleFunction(
  run = function(n = integer(), lambda = double(), zeroProb = double()) {
    returnType(integer())
    isStructuralZero <- rbinom(1, prob = zeroProb, size = 1)
    if (isStructuralZero) return(0)
    return(rpois(1, lambda))
  })

registerDistributions(list(
  dzip = list(
    BUGSdist = "dzip(lambda, zeroProb)",
    discrete = T)
))

### MODEL CODE ##############################################################################################

PVA <- nimbleCode({
  ##### SURVEY MODEL ######################################################################################
  # Poisson likelihood observation model
  for (t in 1:n.occ){
    count[t] ~ dpois(survey[t])
    count_sim[t] ~ dpois(survey[t])
    
    log(survey[t]) <- log(N_tot[1,t] + ext[1,t]) + eps_ds[t]
    eps_ds[t] ~ dnorm(0, sd = sigma_ds)
  }
  
  sigma_ds ~ dgamma(2,4)
  
  # GoF for survey data: Freeman-Tukey test statistics
  ### Freeman–Tukey residuals
  R_srv_r <- sum(pow(sqrt(count[1:n.occ]) - sqrt(survey[1:n.occ]), 2))
  R_srv_p <- sum(pow(sqrt(count_sim[1:n.occ]) - sqrt(survey[1:n.occ]), 2))
  # Bayesian p-value
  Bp_survey <- R_srv_p > R_srv_r
  
  ##### SURVIVAL MODEL ####################################################################################
  # Priors
  beta_phia_int ~ dnorm(0.5, sd = 0.25) ### source of var
  beta_phij_int ~ dnorm(-0.5, sd = 0.25) ### source of var
  beta_p_int ~ dnorm(1, sd = 0.5)
  beta_p_pas ~ dnorm(0, sd = 1)
  
  sigma_p ~ dgamma(2,4)
  
  sigma_tj ~ dgamma(16,40)       ### source of var
  sigma_ta ~ dgamma(16,40)       ### source of var
  
  # Linear models
  for (t in 1:(n.occ+n.prj-1)){ # Here we extend the loop to n.prj more years 
    logit(phi_j[t]) <- beta_phij_int + eps_j[t]
    logit(phi_a[t]) <- beta_phia_int + eps_a[t]
    
    eps_j[t] ~ dnorm(0, sd = sigma_tj)
    eps_a[t] ~ dnorm(0, sd = sigma_ta)
  }
  
  # CJS model with multinomial likelihood
  # Define the multinomial likelihood
  for (t in 1:(n.occ-1)){
    marr.j[t,1:n.occ] ~ dmulti(pr_j[t, 1:n.occ], rel.j[t])
    eps_p[t] ~ dnorm(0, sd = sigma_p)
    
    # Define the cell probabilities of the m-arrays
    # Main diagonal
    logit(p[t]) <- beta_p_int +
      beta_p_pas * passive[t] +
      eps_p[t]
    q[t] <- 1-p[t]
    pr_j[t,t] <- phi_j[t] * p[t]
  }
  for (t in 1:(n.occ-2)){
    # Above main diagonal
    for (j in (t+1):(n.occ-1)){
      pr_j[t,j] <- phi_j[t] * prod(phi_a[(t+1):j]) * prod(q[t:(j-1)]) * p[j]
    } #j
  } #t
  for (t in 2:(n.occ-1)){
    # Below main diagonal
    for (j in 1:(t-1)){
      pr_j[t,j] <- 0
    } #j
  } #t
  # Last column: probability of non-recapture
  for (t in 1:(n.occ-1)){
    pr_j[t,n.occ] <- 1-sum(pr_j[t,1:(n.occ-1)])
  }
  
  # GoF for CMR data: Freeman-Tukey test statistics
  for (t in 1:(n.occ-1)){
    # Simulated m-arrays
    marr.j_sim[t,1:n.occ] ~ dmulti(pr_j[t,1:n.occ], rel.j[t])
    # Expected values and test statistics
    for (j in 1:n.occ){
      marr.j_r[t,j] <- pr_j[t,j] * rel.j[t]
      Exp.j_r[t,j] <- pow(sqrt(marr.j[t,j]) - sqrt(marr.j_r[t,j]), 2)
      Exp.j_p[t,j] <- pow(sqrt(marr.j_sim[t,j]) - sqrt(marr.j_r[t,j]), 2)
    } #j
  } #t
  R_cjs_r <- sum(Exp.j_r[1:(n.occ-1), 1:n.occ])
  R_cjs_p <- sum(Exp.j_p[1:(n.occ-1), 1:n.occ])
  # Bayesian p-value
  Bp_CJS <- R_cjs_p > R_cjs_r
  
  ##### PRODUCTIVITY MODEL ################################################################################
  # Priors
  beta_f_int ~ dnorm(0, sd = 1)
  beta_f_par ~ dnorm(0, sd = 1)
  theta_f ~ dunif(0,1)
  
  beta_pi_int ~ dnorm(0, sd = 1)
  beta_nu_int ~ dnorm(0, sd = 1)
  
  sigma_pi ~ dgamma(16,40)
  sigma_nu ~ dgamma(16,40)
  sigma_tf ~ dgamma(16,40)
  
  # fertility parameters
  for (t in 1:(n.occ+n.prj)) {
    logit(pi[t]) <- beta_pi_int + eps_pi[t]
    logit(nu[t]) <- beta_nu_int + eps_nu[t]
    
    eps_pi[t] ~ dnorm(0, sd = sigma_pi)
    eps_nu[t] ~ dnorm(0, sd = sigma_nu)
    eps_f[t] ~ dnorm(0, sd = sigma_tf)
    
    for (s in 1:n.scn) {
      nu_unbounded[s,t] <- (1 - nu[t] * bfr_mat[s,t])
      # reproductive success rate is >= 0 == 1 (valid) ; success rate is < 0 == 0 (invalid)
      nu_zero[s,t] <- step(nu_unbounded[s,t])
      
      f[s,t] <- exp(beta_f_int + beta_f_par * (pi[t]*par_mat[s,t])+eps_f[t]) * nu_unbounded[s,t]*nu_zero[s,t]
    }
  }
  # Models
  for (i in 1:n.rep) {
    # Probability of reproductive failure
    rep_fail[i] ~ dbern(nu[rep_yr[i]])
  }
  
  for(i in 1:n.nests) {
    # Probability of parasitism
    parasi_i[i] ~ dbern(pi[yr.n[i]])
    # Linear model of reproductive output for nest n (for both sexes)
    log(lam[i]) <- beta_f_int + beta_f_par * parasi_i[i] + eps_f[yr.n[i]]
    # Poisson likelihood model of productivity
    fled_i[i] ~ dzip(lam[i], theta_f)
    fled_i_sim[i] ~ dzip(lam[i], theta_f)
  } # n
  # GoF for clutch size data: Freeman-Tukey test statistics
  R_fled_r <- sum(pow(sqrt(fled_i[1:n.nests]) - sqrt(lam[1:n.nests]), 2))
  R_fled_p <- sum(pow(sqrt(fled_i_sim[1:n.nests]) - sqrt(lam[1:n.nests]), 2))
  # Bayesian p-value
  Bp_fledge <- R_fled_p > R_fled_r
  
  ##### IMMIGRATION MODEL #################################################################################
  # Priors
  beta_om_int ~ dnorm(0, sd = 0.1)
  beta_om_ddi ~ dnorm(0, sd = 0.1)
  sigma_ti ~ dgamma(16,40)
  
  for (t in 1:(n.occ + n.prj - 1)) {
    eps_ti[t] ~ dnorm(0, sd = sigma_ti)
    
    for (s in 1:n.scn) {
      log(om[s,t]) <- 
        beta_om_int + (beta_om_ddi * ((log(N_tot[s,t] + ext[s, t]) - 4.97)/0.56)) + eps_ti[t]
    }
  }
  
  ##### POPULATION MODEL ##################################################################################
  # Process model: our model of population dynamics
  # Binomial and Poisson likelihood models for population size
  N_j[1,1] ~ dpois(Ninit[1])
  N_a[1,1] ~ dpois(Ninit[2])
  N_i[1,1] ~ dpois(Ninit[3])
  
  # Total population size
  for (t in 1:(n.occ+n.prj)){
    for (s in 1:n.scn) {
      N_tot[s,t] <- N_j[s,t] + N_a[s,t] + N_i[s,t]
    }
  }
  
  ##### PROJECTION SCENARIOS ##############################################################################
  # realized scenario
  for (t in 1:n.occ-1){
    N_j[1,t+1] ~ dpois(phi_j[t] * f[1,t]/2 * N_tot[1,t])
    N_a[1,t+1] ~ dbin(phi_a[t], N_tot[1,t])
    N_i[1,t+1] ~ dpois(om[1,t])
  }
  # The past is the same for every future scenario
  for (t in 1:n.occ) {
    for (s in 2:n.scn) {
      N_j[s,t] <- N_j[1,t]
      N_a[s,t] <- N_a[1,t]
      N_i[s,t] <- N_i[1,t]
    }
  }
  # projection scenarios
  for (t in n.occ:(n.occ+n.prj-1)){
    for (s in 1:n.scn) {
      N_j[s,t+1] ~ dpois(phi_j[t] * f[s,t]/2 * N_tot[s,t])
      N_a[s,t+1] ~ dbin(phi_a[t] * srv_mat[s,t] , N_tot[s,t])
      N_i[s,t+1] ~ dpois(om[s,t])
    }
  }
  
  ##### DERIVED PARAMETERS ################################################################################
  
  for (s in 1:n.scn) {
    proj_gr[s] <- pow(N_tot[s,n.occ+n.prj] / N_tot[s,n.occ], 1/(n.prj-1))
    for (t in 1:(n.occ + n.prj - 1)){
      # has the population been extirpated? 1 = Yes ; 0 = No
      ext[s,t] <- equals(N_tot[s, t], 0)
      
      annual_growth_rate[s,t] <- N_tot[s, t + 1] / (N_tot[s, t] + ext[s, t])
    }
    for (t in 1:n.prj) {
      Pr_ext[s, t] <- (N_tot[s, (n.occ + t)] == 0)
      Pr_qext[s, t] <- (N_tot[s, (n.occ + t)] < 20)
      Pr_K[s, t] <- (N_tot[s, (n.occ + t)] > 1000)
    }
  }
  
  # Mean demographic rates
  phia_bar <- mean(phi_a[1:(n.occ-1)])
  phij_bar <- mean(phi_j[1:(n.occ-1)])
  fbar <- mean(f[1,1:n.occ]/2)
  realized_gr <- pow(N_tot[1,n.occ] / N_tot[1,1], 1/(n.occ-1))
  
  # Partial derivatives are calculated by implicit differentiation of matrix model parameters, as described in:
  # # Oli MK, Zinner B. Partial life-cycle analysis: a model for birth-pulse populations.
  # # Ecology 82, 1180–1190 (2001).
  
  ##### MATRIX MODEL DERIVATIVES
  sen_alpha <- fbar*
    exp(log(phij_bar)*alpha)*(realized_gr*(-log(phij_bar)*realized_gr*exp(-log(
      realized_gr)*alpha+log(realized_gr)*omega)+realized_gr*log(realized_gr)*
        exp(-log(realized_gr)*alpha+log(realized_gr)*omega)+log(phij_bar)*phia_bar*
        exp(-log(realized_gr)*alpha+log(realized_gr)*omega)-phia_bar*log(realized_gr)*
        exp(-log(realized_gr)*alpha+log(realized_gr)*omega)-log(phij_bar)*phij_bar*
        exp(-log(realized_gr)*alpha+log(realized_gr)*omega)+log(realized_gr)*phij_bar*
        exp(-log(realized_gr)*alpha+log(realized_gr)*omega)+log(phij_bar)*
        exp(log(phia_bar)*omega-log(phia_bar)*alpha)*phij_bar-
        exp(log(phia_bar)*omega-log(phia_bar)*alpha)*log(phia_bar) * phij_bar)
    )/(
      (-fbar * realized_gr*alpha*
         exp(log(phij_bar)*alpha-log(realized_gr)*alpha+log(realized_gr)*omega)+fbar*phia_bar*alpha*
         exp(log(phij_bar)*alpha-log(realized_gr)*alpha+log(realized_gr)*omega)+fbar*phia_bar*
         exp(log(phij_bar)*alpha-log(realized_gr)*alpha+log(realized_gr)*omega)-fbar*alpha*phij_bar*
         exp(log(phij_bar)*alpha-log(realized_gr)*alpha+log(realized_gr)*omega)-fbar*phij_bar*
         exp(log(phij_bar)*alpha-log(realized_gr)*alpha+log(realized_gr)*omega)+fbar*omega*phij_bar*
         exp(log(phij_bar)*alpha+log(phia_bar)*omega-log(phia_bar)*alpha)+fbar*phij_bar*
         exp(log(phij_bar)*alpha+log(phia_bar)*omega-log(phia_bar)*alpha)-phia_bar*phij_bar*
         exp(log(realized_gr)*omega))
    )
  
  sen_omega <- fbar*realized_gr*phij_bar*log(phia_bar/realized_gr)/
    (-fbar*realized_gr*alpha*
       exp(log(phij_bar)*alpha-log(realized_gr)*alpha+log(realized_gr)*omega)+fbar*phia_bar*alpha*
       exp(log(phij_bar)*alpha-log(realized_gr)*alpha+log(realized_gr)*omega)+fbar*phia_bar*
       exp(log(phij_bar)*alpha-log(realized_gr)*alpha+log(realized_gr)*omega)-fbar*alpha*phij_bar*
       exp(log(phij_bar)*alpha-log(realized_gr)*alpha+log(realized_gr)*omega)-fbar*phij_bar*
       exp(log(phij_bar)*alpha-log(realized_gr)*alpha+log(realized_gr)*omega)+fbar*omega*phij_bar*
       exp(log(phij_bar)*alpha+log(phia_bar)*omega-log(phia_bar)*alpha)+ fbar*phij_bar*
       exp(log(phij_bar)*alpha+log(phia_bar)*omega-log(phia_bar)*alpha)-phia_bar*phij_bar*
       exp(log(realized_gr)*omega))*exp(log(phij_bar)*alpha+log(phia_bar)*omega-log(phia_bar)*alpha)
  
  sen_phij <- fbar*exp(log(phij_bar)*alpha)*realized_gr*
    (realized_gr*alpha*
       exp(-log(realized_gr)*alpha+log(realized_gr)*omega)-realized_gr*
       exp(-log(realized_gr)*alpha+log(realized_gr)*omega)-phia_bar*alpha*
       exp(-log(realized_gr)*alpha+log(realized_gr)*omega)+phia_bar*
       exp(-log(realized_gr)*alpha+log(realized_gr)*omega)+alpha*phij_bar*
       exp(-log(realized_gr)*alpha+log(realized_gr)*omega)-alpha*
       exp(log(phia_bar) * omega-log(phia_bar) * alpha) * phij_bar
    )/(
      phij_bar * 
        (fbar*realized_gr*alpha*
           exp(log(phij_bar)*alpha-log(realized_gr)*alpha+log(realized_gr)*omega)-fbar*phia_bar*alpha*
           exp(log(phij_bar)*alpha-log(realized_gr)*alpha+log(realized_gr)*omega)-fbar*phia_bar*
           exp(log(phij_bar)*alpha-log(realized_gr)*alpha+log(realized_gr)*omega)+fbar*alpha*phij_bar*
           exp(log(phij_bar)*alpha-log(realized_gr)*alpha+log(realized_gr)*omega)+fbar*phij_bar*
           exp(log(phij_bar)*alpha-log(realized_gr)*alpha+log(realized_gr)*omega)-fbar*omega*phij_bar*
           exp(log(phij_bar)*alpha+log(phia_bar)*omega-log(phia_bar)*alpha)-fbar*phij_bar*
           exp(log(phij_bar)*alpha+log(phia_bar)*omega-log(phia_bar)*alpha)+phia_bar*phij_bar*
           exp(log(realized_gr)*omega)))
  
  sen_phia <- -realized_gr*(
    fbar*phia_bar*
      exp(log(phij_bar)*alpha-log(realized_gr)*alpha+log(realized_gr)*omega)+fbar*omega*phij_bar*
      exp(log(phij_bar)*alpha+log(phia_bar)*omega-log(phia_bar)*alpha)-fbar*phij_bar*alpha*
      exp(log(phij_bar)*alpha+log(phia_bar)*omega-log(phia_bar)*alpha)-phia_bar*phij_bar*
      exp(log(realized_gr)*omega)
  )/(
    phia_bar *
      (fbar * realized_gr*alpha*
         exp(log(phij_bar)*alpha-log(realized_gr)*alpha+log(realized_gr)*omega)-fbar*phia_bar*alpha*
         exp(log(phij_bar)*alpha-log(realized_gr)*alpha+log(realized_gr)*omega)-fbar*phia_bar*
         exp(log(phij_bar)*alpha-log(realized_gr)*alpha+log(realized_gr)*omega)+fbar*alpha*phij_bar*
         exp(log(phij_bar)*alpha-log(realized_gr)*alpha+log(realized_gr)*omega)+fbar*phij_bar*
         exp(log(phij_bar)*alpha-log(realized_gr)*alpha+log(realized_gr)*omega)-fbar*omega*phij_bar*
         exp(log(phij_bar)*alpha+log(phia_bar)*omega-log(phia_bar)*alpha)-fbar*phij_bar*
         exp(log(phij_bar)*alpha+log(phia_bar)*omega-log(phia_bar)*alpha)+phia_bar*phij_bar*
         exp(log(realized_gr) * omega)))
  
  sen_fer <- (exp(log(phij_bar)*alpha)*realized_gr*(
    realized_gr*
      exp(-log(realized_gr)*alpha+log(realized_gr)*omega)-phia_bar*
      exp(-log(realized_gr)*alpha+log(realized_gr)*omega)+phij_bar*
      exp(-log(realized_gr)*alpha+log(realized_gr)*omega)-
      exp(log(phia_bar)*omega-log(phia_bar)*alpha)*phij_bar)
  )/(
    (fbar*realized_gr*alpha*
       exp(log(phij_bar)*alpha-log(realized_gr)*alpha+log(realized_gr)*omega)-fbar*phia_bar*alpha*
       exp(log(phij_bar)*alpha-log(realized_gr)*alpha+log(realized_gr)*omega)-fbar*phia_bar*
       exp(log(phij_bar)*alpha-log(realized_gr)*alpha+log(realized_gr)*omega)+fbar*alpha*phij_bar*
       exp(log(phij_bar)*alpha-log(realized_gr)*alpha+log(realized_gr)*omega)+fbar*phij_bar*
       exp(log(phij_bar)*alpha-log(realized_gr)*alpha+log(realized_gr)*omega)-fbar*omega*phij_bar*
       exp(log(phij_bar)*alpha+log(phia_bar)*omega-log(phia_bar)*alpha)-fbar*phij_bar*
       exp(log(phij_bar)*alpha+log(phia_bar)*omega-log(phia_bar)*alpha)+phia_bar*phij_bar*
       exp(log(realized_gr)*omega)))
  
  ela_omega <- sen_omega*omega/realized_gr
  ela_alpha <- sen_alpha*alpha/realized_gr
  ela_phij <- sen_phij*phij_bar/realized_gr
  ela_phia <- sen_phia*phia_bar/realized_gr
  ela_fer <- sen_fer*fbar/realized_gr
})

### PARAMS ##################################################################################################

# Parameters monitored
params <- c(
  "beta_phij_int","phi_j","phij_bar","beta_phia_int","phi_a","phia_bar","beta_om_int","beta_om_ddi","om",
  "beta_f_int","beta_f_par","theta_f","fbar","beta_pi_int","sigma_pi","pi","beta_p_int","beta_p_pas",
  "sigma_p","nu","sigma_nu","beta_nu_int","p","f","N_tot","sigma_ds","annual_growth_rate","realized_gr",
  "Pr_ext","Pr_qext","Pr_K","proj_gr","sen_omega","sen_alpha","sen_phij","sen_phia","sen_fer","ela_omega",
  "ela_alpha","ela_phij","ela_phia","ela_fer","Bp_survey","R_srv_p","R_srv_r","Bp_CJS","R_cjs_p","R_cjs_r",
  "Bp_fledge","R_fled_p","R_fled_r","sigma_ta","sigma_tj","sigma_ti","sigma_tf")

### RUN MODEL ###############################################################################################

samples <- nimbleMCMC(
  code = PVA,
  inits = ybch_ini,
  data = ybch_dat,
  constants = ybch_con,
  monitors = params,
  nchains = 1,
  nburnin = 500000,
  niter = 1500000,
  thin = 100,
  samplesAsCodaMCMC = T
)

### SAVE DATA ###############################################################################################

idx <- length(list.files("./", pattern = "pva_chain"))+1
assign(paste0("pva_chain", idx), value = samples)

save(list = paste0("pva_chain", idx), file = paste0("./pva_chain", idx, ".rda"))
