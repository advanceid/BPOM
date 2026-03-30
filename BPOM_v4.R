# BPOM
# Preparation ----
rm(list = ls())

library(arm)
library(boot)
# library(future.apply)
library(gsDesign)
library(lme4)
library(magrittr)
library(marginaleffects)
library(parallel)
library(parallelly)
library(readxl)
library(rstanarm)
library(tidyverse)



# Environment ----
env <- ifelse(str_detect(getwd(), 'C|D:/'), 1, 2)
t0 <- Sys.time()

current_queue <- Sys.getenv("PBS_QUEUE")
cat("Currently running in queue:", current_queue, "\n")

if (env == 1){
  setwd('D:/NUS Dropbox/Xiangyuan Huang/github/BPOM')
} else if(env == 2){
  setwd('BPOM')
}

treatmentCombination <- read_excel('input/Antibiotic list_24Feb2026.xlsx', sheet = 1) 

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
maxTreat <- str_extract_all(treatmentCombination$intervention, 
                            '[0-9]+') %>% unlist %>% as.numeric %>% unique %>% length

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
  interimMethod <- c('fixed', 'dynamic')   
  nInterim <- 3
} else if (env == 2){
  n_cores <- print(availableCores()) - 2
  cat('Available CPUs: ', n_cores, '\n')
  N_iteration <- 1500
  sampleSize <- seq(400, 3000, 250)
  samplePrior <- 3000
  p1 <- c(0.31, 0.36, 0.51)  # polymyxin
  p2 <- c(0.41)                          # carbapenem
  p3 <- c(0.21, 0.26, 0.31, 0.36)        # cefta
  # primary analysis: polymyxin versus carbapenem(BAT), non-inferirority, 10% margin
  # secondary analysis: cefta versus polymyxin, non-inferiority, 10% margin
  maxit <- 1000
  priorSet <- c('skeptical', 'optimistic', 'noninformative')
  interimMethod <- c('fixed', 'dynamic')   
  # fixed: sticking to the same priors; dynamic: updating prior every interim analysis
  nInterim <- 3
}

scenario <- expand.grid(sampleSize = sampleSize,
                        samplePrior = samplePrior,
                        p1 = p1,
                        p2 = p2,
                        p3 = p3,
                        priorSet = priorSet) %>%
  mutate(scenarioID = paste0('scenario', row_number()))
save(scenario, file = 'results/scenario.RData')

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



simu <- function(scenario){
  for (i in 1:nrow(scenario)){
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
    model <- myTryCatch(
      glm(data = dfRCTAB,
      death ~ groupT + Organism, family = binomial('logit')))
    modelPrior <- myTryCatch(
      glm(data = dfRCTPriorAB, 
          death ~ groupT + Organism, family = binomial('identity')))
    
    if (!(is.null(model$value)|(is.null(modelPrior$value)))){
      model_coef <- summary(model$value)$coefficients %>% 
        as.data.frame %>% rownames_to_column()
      model_coefPrior <- summary(modelPrior$value)$coefficients %>% 
        as.data.frame %>% rownames_to_column()
      
      model_coef %<>% rows_update(model_coefPrior, by = 'rowname')
      if(scenario$priorSet[i] == 'skeptical'){
        model_coef$Estimate[str_detect(model_coef$rowname, 'groupT')] <- 0
        model_coef$`Std. Error` <- 1
      } else if (scenario$priorSet[i] == 'optimistic'){
        model_coef$`Std. Error` <- 0.5
      } else if (scenario$priorSet[i] == 'noninformative'){
        model_coef$Estimate[str_detect(model_coef$rowname, 'groupT')] <- 0
        model_coef$`Std. Error` <- 10
      } 
      current_prior <- model_coef
      } else {break}
    

    for (interim in 1:nInterim){
      dfInterimA <- dfRCTAB %>% filter(slice == interim)
      model1 <- myTryCatch(bayesglm(
        death ~ groupT + Organism, 
        data = dfInterimA, 
        maxit = maxit,
        prior.mean = current_prior$Estimate[-1], 
        prior.scale = current_prior$`Std. Error`[-1],
        prior.mean.for.intercept = current_prior$Estimate[1],
        prior.scale.for.intercept = current_prior$`Std. Error`[1],
        family = binomial(link = "identity")))
      
      dfInterimB <- dfRCTAB %>% filter(slice <= interim)
      model2 <- myTryCatch(bayesglm(
        death ~ groupT + Organism, 
        data = dfInterimB, 
        maxit = maxit,
        prior.mean = model_coef$Estimate[-1], 
        prior.scale = model_coef$`Std. Error`[-1],
        prior.mean.for.intercept = model_coef$Estimate[1],
        prior.scale.for.intercept = model_coef$`Std. Error`[1],
        family = binomial(link = "identity")))
    
      
      # Cefta-polymyxin
      if (!(is.null(model1$value)|(is.null(model2$value)))){
        new_stats <- summary(model1$value)$coefficients %>%
          as.data.frame() %>%
          rownames_to_column("rowname") %>%
          select(rowname, Estimate, `Std. Error`)
        current_prior %<>% rows_update(new_stats, by = 'rowname')

        
        # Secondary analysis
        dfCVP <- dfRCT %>% filter(cAC > 1, groupT %in% c('A', 'C'), slice <= interim)
        model3 <- myTryCatch(glm(data = dfCVP, death ~ groupT + as.factor(combination), 
                                 family = binomial('identity')))
        
        # Extract results
        res1 <- get_safe_res(model1, 'dynamic', interim, scenario[i,], 
                             group_label = 'B', dfInterim = dfInterimA)
        res2 <- get_safe_res(model2, 'fixed', interim, scenario[i,], 
                             group_label = 'B', dfInterim = dfInterimB)
        res3 <- get_safe_res(model3, 'CVP', interim, scenario[i, ], 
                             group_label = 'C', dfInterim = dfCVP)
        if(!is.null(res3)) res3 <- res3[str_detect(res3$term, 'groupT'), ]
        
        res <- rbind(res1, res2, res3) %>%
          mutate(nA = sum(dfInterimB$groupT == 'A'), 
                 nB = sum(dfInterimB$groupT == 'B'),
                 nC = sum(dfCVP$groupT == 'C'))
        assign(paste0(scenario$scenarioID[i], '-', interim), 
               res)
      } else {next}
    }
  }
  resCoef <- ls(pattern = 'scenario[0-9]')
  resCoef <- mget(resCoef)
  resCoef <- do.call(rbind, resCoef)
  rm(list = ls(pattern = 'scenario[0-9]'))
  return(resCoef)
}



# Simulation ----
results <- mclapply(1:N_iteration, function(i) {
  cat('Iteration: ', i, '\n')
  out1 <- simu(scenario)
  out1$iteration_id <- i

  return(list(out1 = out1))
}, 
mc.cores = n_cores)
t1 <- Sys.time()

save(results, file = 'results/final.RData')

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
final_data <- results %>% map_dfr(~ .x$out1)

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