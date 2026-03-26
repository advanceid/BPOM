# BPOM
# Preparation ----
rm(list = ls())

library(arm)
library(boot)
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
  filter(aki != 'Yes') %>%
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
  N_iteration <- 20
  sampleSize <- seq(200, 800, 40)
  samplePrior <- c(30, 100, 400)
  p1 <- 0.4
  p2 <- c(0.5, 0.45, 0.4, 0.38, 0.37, 0.36, 0.35)
  p3 <- 0.4
  maxit <- 1000
  nInterim <- 3
} else if (env == 2){
  n_cores <- print(availableCores()) - 2
  cat('Available CPUs: ', n_cores)
  N_iteration <- 1500
  sampleSize <- seq(1000, 3000, 200)
  samplePrior <- 3000
  p1 <- c(0.41, 0.42, 0.46, 0.51, 0.61)  # polymyxin
  p2 <- c(0.51)                          # carbapenem
  p3 <- c(0.31, 0.36, 0.41, 0.51)        # cefta
  # primary analysis: polymyxin versus carbapenem(BAT), non-inferirority, 10% margin
  # secondary analysis: cefta versus polymyxin, non-inferiority, 10% margin
  maxit <- 1000
  priorScale <- c(1, 0.5, 0.1, 0.05)
  nInterim <- 3
}

scenario <- expand.grid(sampleSize = sampleSize,
                        samplePrior = samplePrior,
                        p1 = p1,
                        p2 = p2,
                        p3 = p3,
                        priorScale = priorScale) %>%
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
gcomp_boot_fn <- function(data, indices) {
  d <- data[indices, ]
  
  fit <- tryCatch({
    bayesglm(
      death ~ factor(groupT) + factor(Organism), 
      data = d, 
      family = binomial(link = "logit"),
      prior.mean = model_coef[2:length(model_coef)],
      prior.scale = priorScale,
      prior.mean.for.intercept = model_coef[1],
      prior.scale.for.intercept = priorScale,
      maxit = maxit
    )
  }, error = function(e) return(NULL))
  
  if (is.null(fit)) return(NA)
  
  levs_groupT <- levels(as.factor(d$groupT))
  levs_org    <- levels(as.factor(d$Organism))
  
  df_a1 <- d %>% mutate(groupT = factor(levs_groupT[2], levels = levs_groupT),
                        Organism = factor(Organism, levels = levs_org))
  
  df_a0 <- d %>% mutate(groupT = factor(levs_groupT[1], levels = levs_groupT),
                        Organism = factor(Organism, levels = levs_org))
  
  # 4. Predict
  p1 <- mean(predict(fit, newdata = df_a1, type = "response"))
  p0 <- mean(predict(fit, newdata = df_a0, type = "response"))
  
  return(p1 - p0)
}

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
    # i <- 10
    sampleSize <- scenario$sampleSize[i]
    samplePrior <- scenario$samplePrior[i]
    priorScale <- scenario$priorScale[i]
    
    # Treatment allocation
    pattern <- patternS
    pattern %<>% mutate(countryAllocation = proportion*sampleSize,
                        organismAllocation = prob*countryAllocation,
                        countryAllocationPrior = proportion*samplePrior,
                        organismAllocationPrior = prob*countryAllocationPrior)
    pattern <- pattern[rep(1:nrow(pattern), each = sampleSize),] %>% 
      arrange(proportion, Organism) %>% group_by(proportion, Organism) %>%
      slice_sample(prop = 1) 
    
    # Randomized outcome event
    dfRCT <- pattern %>%
      filter(row_number() <= organismAllocation) %>%
      group_by(Organism) %>%
      mutate(errO = rnorm(1, 0, 0.01),
             errO = first(errO)) %>% ungroup %>%
      mutate(mortality = if_else(
        groupT == 'A', scenario$p1[i], ifelse(
          groupT == 'B', scenario$p2[i], scenario$p3[i])),
        mortality = mortality + errO) %>%
      group_by(proportion)  %>%
      mutate(slice = ntile(runif(n(), 0, 1), nInterim))
    dfRCT$death <- rbinom(n = nrow(dfRCT), size = 1, prob = dfRCT$mortality)
    
    dfRCTPrior <- pattern %>%
      filter(row_number() <= organismAllocationPrior) %>%
      group_by(Organism) %>%
      mutate(errO = rnorm(1, 0, 0.01),
             errO = first(errO)) %>% ungroup %>%
      mutate(mortality = if_else(
        groupT == 'A', scenario$p1[i], ifelse(
          groupT == 'B', scenario$p2[i], scenario$p3[i])),
        mortality = mortality + errO, 
        slice = ntile(row_number(), nInterim))
    dfRCTPrior$death <- rbinom(n = nrow(dfRCTPrior), 
                               size = 1, prob = dfRCTPrior$mortality)
    
    dfRCTAB <- dfRCT %>% filter(cAB > 1, groupT %in% c('A', 'B'))
    dfRCTPriorAB <- dfRCTPrior %>% filter(cAB > 1, groupT %in% c('A', 'B'))
    
    
    for (interim in 1:nInterim){
      # Risk difference model 
      model <- myTryCatch(
        glm(death ~ groupT + Organism, family = binomial('identity'), 
            data = dfRCTAB,
            mustart = rep(mean(dfRCTAB$death), nrow(dfRCTAB))))
      model_coef <- coefficients(model[[1]])
      
      modelPrior <- myTryCatch(
        glm(data = dfRCTPriorAB, 
            death ~ groupT + Organism, family = binomial('identity')))
      modelPrior_coef <- coefficients(modelPrior[[1]])
      
      model_coef[names(model_coef)] <- 0
      modelPrior_coef <- modelPrior_coef[
        names(modelPrior_coef) %in% names(model_coef)]
      model_coef[names(modelPrior_coef)] <- modelPrior_coef
      model_coef[is.na(model_coef)] <- 0
      
      dfInterim <- dfRCTAB %>% filter(slice <= interim)
      model1 <- myTryCatch(bayesglm(
        death ~ groupT + Organism, 
        data = dfInterim, 
        maxit = maxit,
        prior.mean = model_coef[2:length(model_coef)],
        prior.scale = priorScale,
        prior.mean.for.intercept = model_coef[1],
        prior.scale.for.intercept = priorScale,
        family = binomial(link = "identity")))
      model2 <- myTryCatch(bayesglm(
        death ~ groupT + Organism, 
        data = dfInterim, 
        maxit = maxit,
        prior.mean = 0,
        prior.scale = priorScale,
        prior.mean.for.intercept = 0,
        prior.scale.for.intercept = priorScale,
        family = binomial(link = "identity")))
      
      # G-computation
      # boot_results <- boot(data = dfInterim, statistic = gcomp_boot_fn, R = 2000)
      # if (any(is.na(boot_results$t))) {
      #   valid_indices <- which(!is.na(boot_results$t))
      #   
      #   if (length(valid_indices) < 10) {
      #     boot_ci <- list(percent = rep(NA, 5)) 
      #   } else {
      #     boot_results$t <- matrix(boot_results$t[valid_indices], ncol = 1)
      #     boot_results$R <- length(valid_indices) 
      #     boot_ci <- boot.ci(boot_results, type = "perc", conf = p_nomi[interim])
      #   }
      # } else {
      #   boot_ci <- boot.ci(boot_results, type = "perc", conf = p_nomi[interim])
      # }
      
      # Cefta-polymyxin
      dfCVP <- dfRCT %>% filter(cAC > 1, 
                                groupT %in% c('A', 'C'), 
                                slice <= interim)
      model4 <- myTryCatch(
        glm(data = dfCVP, death ~ groupT + as.factor(combination), 
            family = binomial('identity')))
      
      
      if (!(is.null(model1$value)|is.null(model2$value))){
        res1 <- get_safe_res(model1, 'goodPrior', interim, scenario[i,],
                             group_label = 'B', dfInterim = dfInterim)
        res2 <- get_safe_res(model2, 'badPrior', interim, scenario[i,],
                             group_label = 'B', dfInterim = dfInterim)
        res4 <- get_safe_res(model4, 'CVP', interim, scenario[i, ], 
                             group_label = 'C', dfInterim = dfCVP)
        res4 <- res4[str_detect(res4$term, 'groupT'), ]
        
        
        # g_rd_est <- boot_results$t0
        # g_rd_low <- boot_ci$percent[4]
        # g_rd_upp <- boot_ci$percent[5]
        # res3 <- data.frame(
        #   term = "groupTB",
        #   Estimate = g_rd_est,
        #   std = (g_rd_upp - g_rd_low) / (2 * 1.96), # Approximated SE from CI
        #   zVal = NA, 
        #   pVal = NA,
        #   model = 'gComp',
        #   interim = interim,
        #   scenarioID = scenario$scenarioID[i],
        #   interve = sum(dfInterim$groupT == 'B'),
        #   control = sum(dfInterim$groupT == 'A'),
        #   conf_low = g_rd_low,
        #   conf_high = g_rd_upp)
        # print(res1)
        # print(res2)
        # print(res4)
        res <- rbind(res1, res2, res4)
        # print(res)
        assign(paste0(scenario$scenarioID[i], '-', interim), 
               res)
      }
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

resultSum %<>% mutate(
  mortalityPoly = paste0('Poly - BAT = ', 100*(p1 - p2), '%'),
  mortalityCeft = paste0('Ceft - Poly = ', 100*(p3 - p1), '%'),
  priorScale = as.character(priorScale)) %>%
  arrange(mortalityPoly, mortalityCeft)
save(resultSum, file = 'results/scenarioPower.RDtata')  

for (modelIndex in unique(resultSum$model)){
  print(modelIndex)
  resultSum %>% filter(model == modelIndex) %>%
    ggplot(data = ., aes(x = sampleSize, y = power)) +
    geom_point(aes(shape = priorScale, col = priorScale)) +
    geom_hline(yintercept = c(5, 80), linetype = 2) +
    facet_wrap(~mortalityPoly+mortalityCeft, ncol = 4) +
    labs(x = 'Sample Size', y = 'Power (%)', col = 'Prior Scale', shape = 'Prior Scale')
  ggsave(paste0('results/Power-SampleSize', modelIndex, '.tiff'),
         dpi = 300, width = 15, height = 12)
}
