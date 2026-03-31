# BPOM
# Preparation ----
rm(list = ls())

library(arm)
library(boot)
library(brms)
library(cmdstanr)
library(gsDesign)
library(lme4)
library(magrittr)
library(marginaleffects)
library(parallel)
library(parallelly)
library(readxl)
library(rstanarm)
library(tidybayes)
library(tidyverse)



# Environment ----
env <- ifelse(str_detect(getwd(), 'C|D:/'), 1, 2)
t0 <- Sys.time()

set.seed(123)
current_queue <- Sys.getenv("PBS_QUEUE")
cat("Currently running in queue:", current_queue, "\n")

if (env == 1){
  setwd('D:/NUS Dropbox/Xiangyuan Huang/github/BPOM')
} else if(env == 2){
  setwd('BPOM')
}

if (sum(str_detect(dir(), 'cmdstan')) == 0){
  install_cmdstan(dir = 'cmdstan', cores = 2, overwrite = T)
} else {
  set_cmdstan_path(dir('cmdstan', full.names = T)[1])
  check_cmdstan_toolchain()
}



# Pattern ----
treatmentCombination <- read_excel(
  'input/Antibiotic list_24Feb2026.xlsx', sheet = 1) 

treatmentCombination %<>% rename(intervention = 'Intervention arms',
                                 aki = 'Acute kidney injury') %>%
  filter(aki != 'Yes', !(Country %in% c('Australia', 'Spain'))) %>%
  mutate(combination = row_number(), 
         Organism = toupper(Organism),
         Organism = ifelse(Organism == 'CRE',
                           ifelse(str_detect(`Resistance gene`, 'A'),
                                  paste0(Organism, '_A'),
                                  paste0(Organism, '_B')),
                           Organism))
maxTreat <- str_extract_all(treatmentCombination$intervention,'[0-9]+') %>% 
  unlist %>% as.numeric %>% unique %>% length

patternS <- treatmentCombination %>%
  mutate(arm = str_split(intervention, "\\r?\\n")) %>%
  unnest(arm) %>%
  mutate(arm = str_remove(arm, "^\\d+\\.\\s*"),
         arm = str_trim(arm),
         groupT = ifelse(str_detect(arm, 'Colistin'), 'A',
                         ifelse(str_detect(arm, 'Meropenem'),
                                'B', 'C')),
         prob = recode(Organism, 
                       'CRAB' = 0.5,
                       'CRE_A' = 0.15,
                       'CRE_B' = 0.15,
                       'CRPA' = 0.2)) %>%
  group_by(combination) %>%
  mutate(nInt = n(), 
         nTreat = min(sum(groupT == 'A'), 1) + 
           min(sum(groupT == 'B'), 1) + 
           min(sum(groupT == 'C'), 1),
         cAB = min(sum(groupT == 'A'), 1) + min(sum(groupT == 'B'), 1),
         cAC = min(sum(groupT == 'A'), 1) + min(sum(groupT == 'C'), 1),
         cBC = min(sum(groupT == 'B'), 1) + min(sum(groupT == 'C'), 1),
         proportion = ifelse(Country == 'Indonesia', 0.2, 0.8)) %>%
  filter(nTreat > 1) 



# Scenario ----
if (env == 1){
  n_cores <- 1
  cat('Available CPUs: ', n_cores)
  N_iteration <- 20
  sampleSize <- seq(400, 3000, 200)
  samplePrior <- 3000
  p1 <- c(0.31, 0.36, 0.51)  # polymyxin
  p2 <- c(0.41)                          # carbapenem
  p3 <- c(0.21, 0.26, 0.31, 0.36)        # cefta
  # primary analysis: polymyxin versus carbapenem(BAT), non-inferirority, 10% margin
  # secondary analysis: cefta versus polymyxin, non-inferiority, 10% margin
  maxit <- 1000
  priorSet <- c('skeptical', 'optimistic', 'noninformative')
  nInterim <- 3
} else if (env == 2){
  n_cores <- print(availableCores()) - 2
  cat('Available CPUs: ', n_cores, '\n')
  N_iteration <- 30
  sampleSize <- seq(800, 2000, 300)
  samplePVB <- sampleSize*0.2*0.72
  samplePrior <- 3000
  p1 <- c(0.31, 0.36, 0.41)  # polymyxin
  p2 <- c(0.41)                          # carbapenem
  p3 <- c(0.26)        # cefta
  # primary analysis: polymyxin versus carbapenem(BAT), non-inferirority, 10% margin
  # secondary analysis: cefta versus polymyxin, non-inferiority, 10% margin
  maxit <- 1000
  priorSet <- c('skeptical')
  # fixed: sticking to the same priors; dynamic: updating prior every interim analysis
  nInterim <- 3
}

scenario <- expand.grid(sampleSize = sampleSize,
                        samplePrior = samplePrior,
                        p1 = p1,
                        p2 = p2,
                        p3 = p3,
                        priorSet = priorSet) %>%
  mutate(scenarioID = paste0('scenario', row_number())) %>%
  filter(!((p1 == 0.51)&(p3 %in% c(0.26, 0.31, 0.36))))



# Result export ----
new_folder <- paste0('results/run', format(Sys.time(), "%Y-%m-%d_%H-%M-%S"))
dir.create(new_folder)
setwd(new_folder)
save(scenario, file = 'cenario.RData')



# Interim analysis O-Brien Fleming
design <- gsDesign(k = 3, test.type = 1, alpha = 0.05, beta = 0.2,
                   timing = c(0.33, 0.66, 1.0), sfu = 'OF', sfupar = 0)
z_Val <- -design$upper$bound
cat('Z values:', z_Val, '\n')
p_nomi <- 1 - 2*pnorm(z_Val)
cat('P values:', p_nomi, '\n')



# Function mytrycatch ----
myTryCatch <- function(expr) {
  warn_list <- character()
  err_msg <- NULL

  value <- tryCatch(
    withCallingHandlers(
      expr,
      warning = function(w) {
        warn_list <<- c(warn_list, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      err_msg <<- conditionMessage(e)
      return(NULL)
    }
  )

  # Deep-nested warnings for lme4/glmer models
  if (!is.null(value)) {
    if (inherits(value, "glmerMod")) {
      opt_msg <- value@optinfo$conv$lme4$message
      if (!is.null(opt_msg) && !isSingular(value)) {
        warn_list <- c(warn_list, opt_msg)
      }
      if (length(value@optinfo$warnings) > 0) {
        warn_list <- c(warn_list, unlist(value@optinfo$warnings))
      }
    }
  }

  list(
    value = value,
    warn  = if (length(warn_list) > 0) unique(warn_list) else NULL,
    error = err_msg
  )
}

myTryCatch_brms <- function(expr) {
  warn_list <- character()
  err_msg <- NULL
  
  # 1. Standard R-level Catching
  value <- tryCatch(
    withCallingHandlers(
      expr,
      warning = function(w) {
        warn_list <<- c(warn_list, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) {
      err_msg <<- conditionMessage(e)
      return(NULL)
    }
  )
  
  # 2. Deep-nested Bayesian Diagnostics for brms
  if (!is.null(value) && inherits(value, "brmsfit")) {
    
    # Extract Stan sampler diagnostics (Divergences)
    # This identifies if the MCMC failed to explore the posterior properly
    sampler_params <- rstan::get_sampler_params(value$fit, inc_warmup = FALSE)
    n_divergent <- sum(sapply(sampler_params, function(x) sum(x[, "divergent__"])))
    
    if (n_divergent > 0) {
      warn_list <- c(warn_list, paste("There were", n_divergent, "divergent transitions after warmup."))
    }
    
    # Check R-hat (Convergence)
    # R-hat > 1.05 or 1.1 suggests the chains haven't converged
    rhats <- brms::rhat(value)
    if (any(rhats > 1.05, na.rm = TRUE)) {
      warn_list <- c(warn_list, "R-hat values indicate lack of convergence (R-hat > 1.05).")
    }
    
    # Check Effective Sample Size (ESS)
    # Low ESS suggests the "Full Posterior" estimate is noisy
    ess_vals <- brms::neff_ratio(value)
    if (any(ess_vals < 0.1, na.rm = TRUE)) {
      warn_list <- c(warn_list, "Bulk or Tail ESS ratios are low (< 0.1).")
    }
  }
  
  list(
    value = value,
    warn  = if (length(warn_list) > 0) unique(warn_list) else NULL,
    error = err_msg
  )
}



# Simulation engine ----
get_safe_res <- function(model_obj, model_label, interim_val,
                         scenario_row, group_label = 'B', dfInterim = dfInterim) {
  if (is.null(model_obj$value)) return(NULL)

  tryCatch({
    summ <- summary(model_obj$value)
    coefs <- as.data.frame(summ$coefficients) %>%
      rownames_to_column("term") %>%
      rename(Estimate = 'Estimate',
        zVal = 'z value', std = 'Std. Error', pVal = 'Pr(>|z|)') %>%
      mutate(
        model = model_label,
        interim = interim_val,
        scenarioID = scenario_row$scenarioID,
        interve = sum(dfInterim$groupT == group_label),
        control = sum(dfInterim$groupT == 'A'),
        conf_low = NA,
        conf_high = NA
      )
    return(coefs)
  }, error = function(e) return(NULL))
}

posteriorExtract <- function(modelBRM){
  post_summary <- as_draws_df(modelBRM) %>% as.data.frame %>% 
    dplyr::select(starts_with('b_')) %>%
    pivot_longer(cols = starts_with('b_'), 
                 names_to = 'rowname', values_to = 'value') %>%
    mutate(rowname = str_replace(rowname, 'b_', '')) %>%
    group_by(rowname) %>% 
    summarise(Estimate = mean(value), Est.Error = sd(value))
  return(post_summary)
}

coefToPrior <- function(coef_df) {
  prior_list <- list()
  
  for(i in 1:nrow(coef_df)) {
    dist_str <- paste0("normal(", coef_df$Estimate[i], ", ", coef_df$Est.Error[i], ")")
    
    if(coef_df$rowname[i] == 'Intercept') {
      # Use prior_string to handle the dynamic variables
      prior_list[[i]] <- prior_string(dist_str, class = "Intercept")
    } else {
      prior_list[[i]] <- prior_string(dist_str, class = "b", coef = coef_df$rowname[i])
    }
  }
  final_prior <- do.call(c, prior_list)
  return(final_prior)
}


simu <- function(iteration, scenario){
  for (i in 1:nrow(scenario)){
    cat('Scenario:', scenario$scenarioID[i])
    
    sampleSize <- scenario$sampleSize[i]
    samplePrior <- scenario$samplePrior[i]
    priorScale <- scenario$priorScale[i]
    
    
    # Treatment allocation
    pattern <- patternS %>%
      crossing(priordata = c('trial', 'prior')) %>% 
      mutate(countryAllocation = ifelse(priordata == 'trial',
                                        proportion*sampleSize,
                                        proportion*samplePrior),
             organismAllocation = countryAllocation*prob)

    pattern <- pattern[rep(1:nrow(pattern), each = max(sampleSize, samplePrior)),]  
    
    pattern %<>% arrange(priordata, proportion, Organism) %>% 
      group_by(priordata, proportion, Organism) %>%
      slice_sample(prop = 1) %>%
      filter(row_number() <= organismAllocation) %>% 
      group_by(priordata, combination) %>%
      mutate(errO = rnorm(1, 0, 0.01),
             errO = first(errO)) %>% ungroup %>%
      mutate(mortality = if_else(groupT == 'A', scenario$p1[i], 
                                 ifelse(groupT == 'B', scenario$p2[i], 
                                        scenario$p3[i])),
        mortality = mortality + errO) %>%
      group_by(priordata, proportion)  %>%
      mutate(slice = ntile(runif(n(), 0, 1), nInterim))
    pattern$death <- rbinom(n = nrow(pattern), 
                            size = 1, prob = pattern$mortality)
      
    
    # Randomized outcome event
    dfRCT <- pattern %>% filter(priordata == 'trial')
    dfRCTPrior <- pattern %>% filter(priordata == 'prior')

    dfRCTAB <- dfRCT %>% filter(cAB > 1, groupT %in% c('A', 'B'))
    dfRCTPriorAB <- dfRCTPrior %>% filter(cAB > 1, groupT %in% c('A', 'B'))
    
    # Priors
    model <- myTryCatch_brms(brm(
      data = dfRCTAB, family = bernoulli('logit'),
      death ~ groupT + Organism, backend = 'cmdstanr', silent = 1,
      chains = 1, iter = 1000, warmup = 500, refresh = 0))
    modelPrior <- myTryCatch_brms(brm(
      data = dfRCTPriorAB, family = bernoulli('logit'),
      death ~ groupT + Organism, backend = 'cmdstanr', silent = 1,
      chains = 1, iter = 1000, warmup = 500, refresh = 0))

    
    if (!(is.null(model$value)|(is.null(modelPrior$value)))){
      model_coef <- posteriorExtract(model$value)
      model_coefPrior <- posteriorExtract(modelPrior$value)
      model_coef %<>% rows_update(
        model_coefPrior, by = 'rowname', unmatched = "ignore")
      
      if(scenario$priorSet[i] == 'skeptical'){
        model_coef$Estimate[str_detect(model_coef$rowname, 'groupT')] <- 0
        model_coef$Est.Error <- 1
      } else if (scenario$priorSet[i] == 'optimistic'){
        model_coef$Est.Error <- 0.5
      } else if (scenario$priorSet[i] == 'noninformative'){
        model_coef$Estimate[str_detect(model_coef$rowname, 'groupT')] <- 0
        model_coef$Est.Error <- 10
      } 
      current_coef <- model_coef
    } else {break}
    original_prior <- coefToPrior(current_coef)

    
    for (interim in 1:nInterim){
      current_prior <- coefToPrior(current_coef)
      
      dfInterimA <- dfRCTAB %>% filter(slice == interim)
      model1 <- myTryCatch_brms(brm(
          data = dfInterimA, family = bernoulli('logit'), prior = current_prior,
          death ~ groupT + Organism, backend = 'cmdstanr', silent = 1,
          chains = 4, iter = 1000, warmup = 500, refresh = 0))
      
      dfInterimB <- dfRCTAB %>% filter(slice <= interim)
      model2 <- myTryCatch_brms(brm(
        data = dfInterimB, family = bernoulli('logit'), prior = original_prior,
        death ~ groupT + Organism, backend = 'cmdstanr', silent = 1,
        chains = 4, iter = 1000, warmup = 500, refresh = 0))
      
      dfCVP <- dfRCT %>% filter(cAC > 1, groupT %in% c('A', 'C'), slice <= interim)
      model3 <- myTryCatch(glm(data = dfCVP, death ~ groupT + as.factor(combination), 
                               family = binomial('identity')))
      res3 <- get_safe_res(model3, 'CVP', interim, scenario[i, ], 
                           group_label = 'C', dfInterim = dfCVP)
    
      
      if (!(is.null(model1$value)|(is.null(model2$value)))){
        # Update prior
        model1_coef <- fixef(model1$value) %>% as.data.frame %>% 
          rownames_to_column() %>% 
          dplyr::select(rowname, Estimate, Est.Error)
        current_coef %<>% rows_update(model1_coef, by = 'rowname')
        
        
        # posterior extract
        rd_draws1 <- avg_comparisons(model1$value, 
                                     variables = 'groupT',
                                    comparison = 'difference',
                                    type = "response")
        rd_vector1 <- posterior_draws(rd_draws1)$draw
        save(rd_vector1, file = 
               paste0('pos-', scenario$scenarioID[i], '-m1-intrim', 
                      interim, '-iteration', iteration, '.RData'))

        rd_draws2 <- avg_comparisons(model2$value, 
                                     variables = 'groupT',
                                     comparison = 'difference',
                                     type = "response")
        rd_vector2 <- posterior_draws(rd_draws2)$draw
        save(rd_vector2, file = 
               paste0('pos-', scenario$scenarioID[i], '-m2-intrim', 
                      interim,  '-iteration', iteration, '.RData'))

        # Extract results
        # res1 <- get_safe_res(model1, 'dynamic', interim, scenario[i,], 
        #                      group_label = 'B', dfInterim = dfInterimA)
        # res2 <- get_safe_res(model2, 'fixed', interim, scenario[i,], 
        #                      group_label = 'B', dfInterim = dfInterimB)
        
        # ceft versus polymyxin
       
        if(!is.null(res3)){
          res3 <- res3[str_detect(res3$term, 'groupT'), ]
        
        # res <- rbind(res1, res2, res3) %>%
        res <- res3 %>%
          mutate(nAIND = sum(dfInterimB$groupT == 'A'), 
                 nBIND = sum(dfInterimB$groupT == 'B'),
                 nAGLO = sum(dfCVP$groupT == 'A'),
                 nCGLO = sum(dfCVP$groupT == 'C'))
        assign(paste0(scenario$scenarioID[i], '-', interim), 
               res)
        }
      } else {next}
    }
  }
  resCoef <- ls(pattern = 'scenario[0-9]')
  if (length(resCoef) > 0) {
    resCoef <- mget(resCoef)
    resCoef <- do.call(rbind, resCoef)
    rm(list = ls(pattern = 'scenario[0-9]'))
  } else {
    return(NULL) # Or an empty data frame with your expected columns
  }
}



# Simulation ----
results <- mclapply(1:N_iteration, function(iteration) {
  cat('Iteration: ', iteration, '\n')
  out1 <- simu(iteration, scenario)
  out1$iteration_id <- iteration

  return(list(out1 = out1))
}, 
mc.cores = n_cores)
t1 <- Sys.time()

save(results, file = 'final.RData')

cat('Number of scenarios: ', nrow(scenario), '\n')
cat('Number of iterations: ', N_iteration, '\n')
print(t1 - t0)

q()



# Windows parallelization ----
# N_iteration <- 
cat("\nStarting Parallel Simulation: ", N_iteration, " iterations per scenario.\n")
plan(multisession, workers = max(1, n_cores))

t0 <- Sys.time()
all_results_list <- list()
scenario_iterations <- future_lapply(1:N_iteration, function(x) {
  # Call your existing simu function for just this one scenario row
  cat(x, ":", Sys.time())
  res <- simu(scenario)
  
  # Add an iteration ID column so you can distinguish them later
  if (!is.null(res)) res$iteration <- x
  return(res)
}, future.seed = TRUE) # Essential for valid random numbers in parallel

result <- do.call(rbind, scenario_iterations)
t1 <- Sys.time()
t1 - t0

# Load data ----
load('results/final.RData')
load('results/scenario.RData')



# Post-processing ----
final_data <- results %>%
  keep(~ !inherits(.x, "try-error")) %>%
  map_dfr(~ .x$out1)

resultSum <- final_data %>% filter(str_detect(term, 'groupT')) %>%
  mutate(zThreshold = ifelse(interim == 1, z_Val[1], 
                             ifelse(interim == 2, z_Val[2],
                                    ifelse(interim == 3, z_Val[3], NA)))) %>%
  filter(!is.na(zThreshold)) %>%
  group_by(model, scenarioID, iteration_id) %>%
  mutate(count = n()) %>% filter(count == 3)
# save(resultSum, file = 'results/resultSum.RData')
# load('results/resultSum.RData')
resultSum %<>%
  mutate(LL = Estimate + zThreshold*std,  # Z value as negative number
         UL = Estimate - zThreshold*std,
         pass = ifelse(model == 'CVP', UL <= 0.1, LL >= -0.1)) %>%
  summarise(positive = sum(pass) > 0) %>%
  group_by(scenarioID, model) %>%
  summarise(nPOS = sum(positive),
            N_iteration = n(),
            power = 100*nPOS/N_iteration) %>%
  left_join(scenario, by = 'scenarioID')



# polymyxin
PVB <- resultSum %>% filter(model != 'CVP') %>%
  arrange(model, priorSet, p2, p1, sampleSize) %>% 
  group_by(model, priorSet, p2, p1, sampleSize) %>%
  summarise(N_iteration = sum(N_iteration),
            nPOS = sum(nPOS),
            power = 100*nPOS/N_iteration)

PVB %>% mutate(mortalityPoly = paste0('Poly - BAT = ', 100*(p1 - p2), '%')) %>%
  ggplot(data = ., aes(x = sampleSize, y = power)) + 
  geom_point(aes(shape = priorSet, col = priorSet)) +
  geom_hline(yintercept = c(5, 80), linetype = 2) +
  facet_wrap(~mortalityPoly+model, ncol = 2) 
  
# sample size of polymyxin+BAT and power
PVB %>% mutate(mortalityPoly = paste0('Poly - BAT = ', 100*(p1 - p2), '%'),
               sampleSizePB = sampleSize*0.2*0.72) %>%
  ggplot(data = ., aes(x = sampleSizePB, y = power)) + 
  geom_point(aes(shape = priorSet, col = priorSet)) +
  geom_hline(yintercept = c(5, 80), linetype = 2) +
  facet_wrap(~mortalityPoly+model, ncol = 2) 


# sample size and power for ceft versus power 
CVP <- resultSum %>% filter(model == 'CVP') %>%
  arrange(model, priorSet, p1, p3, sampleSize) %>% 
  group_by(model, priorSet, p1, p3, sampleSize) %>%
  summarise(N_iteration = sum(N_iteration),
            nPOS = sum(nPOS),
            power = 100*nPOS/N_iteration)

CVP %>% mutate(mortalityCeft = paste0('Ceft - BAT = ', 100*(p3 - p1), '%')) %>%
  ggplot(data = ., aes(x = sampleSize, y = power)) + 
  geom_point(aes(shape = priorSet, col = priorSet)) +
  geom_hline(yintercept = c(5, 80), linetype = 2) +
  facet_wrap(~mortalityPoly+model, ncol = 2) 