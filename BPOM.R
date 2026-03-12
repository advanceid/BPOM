# BPOM
# Preparation ----
rm(list = ls())

library(arm)
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

if (env == 1){
  setwd('D:/NUS Dropbox/Xiangyuan Huang/github/BPOM')
} else if(env == 2){
  setwd('BPOM')
}

treatmentCombination <- read_excel('Antibiotic list_Indo.xlsx', sheet = 1)

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
  group_by(Organism) %>%
  mutate(nInt = n())



# Scenario ----
if (env == 1){
  n_cores <- 1
  N_iteration <- 20
  sampleSize <- seq(200, 800, 40)
  samplePrior <- c(30, 100, 300)
  p1 <- 0.4
  p2 <- c(0.5, 0.45, 0.4, 0.38, 0.37, 0.36, 0.35)
  p3 <- 0.4
  maxit <- 1000
  nInterim <- 3
  } else if (env == 2){
    n_cores <- 100
    N_iteration <- 2000
    sampleSize <- seq(200, 1200, 20)
    samplePrior <- 500
    p1 <- c(0.5, 0.45, 0.4, 0.38, 0.37, 0.36, 0.35)
    p2 <- 0.4
    p3 <- 0.4
    maxit <- 1000
    nInterim <- 3
    }

scenario <- expand.grid(sampleSize = sampleSize,
                        samplePrior = samplePrior,
                        p1 = p1,
                        p2 = p2, 
                        p3 = p3) %>%
  mutate(scenarioID = paste0('scenario', row_number()))
save(scenario, file = 'scenario.RData')



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
simu <- function(scenario, comparator){
  for (i in 1:nrow(scenario)){
    sampleSize <- scenario$sampleSize[i]
    samplePrior <- scenario$samplePrior[i]
    pattern <- patternS
    pattern$allocation <- 
      apply(rmultinom(sampleSize, size = 1, patternS$prob/patternS$nInt), 1, sum)
    pattern$allocationPrior <- 
      apply(rmultinom(samplePrior, size = 1, patternS$prob/patternS$nInt), 1, sum)
    pattern %<>% group_by(Organism) %>%
      mutate(errO = rnorm(1, 0, 0.01),
             errO = first(errO)) %>% ungroup %>%
      mutate(
        mortality = if_else(
          groupT == 'A', scenario$p1[i], ifelse(
            groupT == 'B', scenario$p2[i], scenario$p3[i])),
        mortality = mortality + errO)
    
    dfRCT <- pattern[rep(1:nrow(pattern), pattern$allocation), ]
    dfRCT$death <- rbinom(n = nrow(dfRCT), size = 1, prob = dfRCT$mortality)
    dfRCT %<>% slice_sample(prop = 1) %>%
      mutate(slice = ntile(row_number(), nInterim))
    dfRCT %<>% filter(groupT %in% c('A', comparator))
    
    dfRCTPrior <- pattern[rep(1:nrow(pattern), pattern$allocationPrior), ]
    dfRCTPrior$death <- rbinom(n = nrow(dfRCTPrior), size = 1, prob = dfRCTPrior$mortality)
    dfRCTPrior %<>% filter(groupT %in% c('A', comparator))
    
    model <- myTryCatch(
      glm(death ~ groupT + Organism, family = binomial('identity'), 
          data = dfRCT,
          mustart = rep(mean(dfRCT$death), nrow(dfRCT))))
    model_coef <- coefficients(model[[1]])
    
    modelPrior <- myTryCatch(
      glm(data = dfRCTPrior, 
          death ~ groupT + Organism, family = binomial('identity')))
    modelPrior_coef <- coefficients(modelPrior[[1]])
    
    model_coef[names(model_coef)] <- 0
    modelPrior_coef <- modelPrior_coef[
      names(modelPrior_coef) %in% names(model_coef)]
    model_coef[names(modelPrior_coef)] <- modelPrior_coef
    model_coef[is.na(model_coef)] <- 0
    
    prior <- normal(location = model_coef,
                    scale = rep(1, length(model_coef)))
    
    for (interim in 1:nInterim){
      dfInterim <- dfRCT %>% filter(slice <= interim)
      
      model1 <- myTryCatch(bayesglm(
        death ~ groupT + Organism, 
        data = dfInterim, 
        maxit = maxit,
        prior.mean = model_coef[2:5],
        prior.scale = 1,
        prior.mean.for.intercept = model_coef[1],
        prior.scale.for.intercept = 1,
        family = binomial(link = "identity")))
      model2 <- myTryCatch(bayesglm(
        death ~ groupT + Organism, 
        data = dfInterim, 
        maxit = maxit,
        prior.mean = 0,
        prior.scale = 1,
        prior.mean.for.intercept = model_coef[1],
        prior.scale.for.intercept = 1,
        family = binomial(link = "identity")))
      
      if (!(is.null(model1$value)|is.null(model2$value))){
      res1 <- model1[[1]] %>% summary %>% coefficients %>% as.data.frame %>%
        rownames_to_column("term") %>% 
        rename(zVal = 'z value', std = 'Std. Error') %>%
        mutate(model = 'goodPrior',
               interim = interim,
               scenarioID = scenario$scenarioID[i],
               interve = sum(dfInterim$groupT == 'A'),
               control = sum(dfInterim$groupT == comparator))
      
      res2 <- model2[[1]] %>% summary %>% coefficients %>% as.data.frame %>%
        rownames_to_column("term") %>% 
        rename(zVal = 'z value', std = 'Std. Error') %>%
        mutate(model = 'badPrior',
               interim = interim,
               scenarioID = scenario$scenarioID[i],
               interve = sum(dfInterim$groupT == 'A'),
               control = sum(dfInterim$groupT == comparator))
      
      res <- rbind(res1, res2)
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
  print(i)
  out1 <- simu(scenario, comparator = 'B')
  out1$iteration_id <- i
  out1$comparator <- 'B'
  
  out2 <- simu(scenario, comparator = 'C')
  out2$iteration_id <- i
  out2$comparator <- 'C'
  return(list(out1 = out1, out2 = out2))
  }, 
  mc.cores = n_cores)
t1 <- Sys.time()

save(results, file = 'final.RData')

cat('Number of scenarios: ', nrow(scenario), '\n')
print(t1 - t0)

q()



# Load data ----
load('final.RData')
load('scenario.RData')



# Post-processing ----
# Interim analysis O-Brien Fleming
design <- gsDesign(k = 3, test.type = 1, alpha = 0.05, beta = 0.2,
                   timing = c(0.33, 0.66, 1.0), sfu = 'OF', sfupar = 0)
z_Val <- -design$upper$bound
cat('Z values:', z_Val)

final_data <- map(results, ~ bind_rows(.x$out1, .x$out2)) %>% 
  bind_rows()


resultSum <- final_data %>% filter(str_detect(term, 'groupT')) %>%
  mutate(z_Val = ifelse(interim == 1, z_Val[1], 
                        ifelse(interim == 2, z_Val[2],
                               ifelse(interim == 3, z_Val[3], NA)))) %>%
  filter(!is.na(z_Val)) %>%
  group_by(comparator, model, scenarioID, iteration_id) %>%
  mutate(count = n()) %>% filter(count == 3) %>%
  mutate(LL = Estimate + z_Val*std,
         pass = LL >= -0.1) %>%
  summarise(positive = sum(pass) > 0) %>%
  group_by(comparator, model, scenarioID) %>%
  summarise(nPOS = sum(positive),
            N_iteration = n(),
            power = 100*nPOS/N_iteration) %>%
  left_join(scenario, by = 'scenarioID') %>%
  mutate(p1 = as.character(100*p1))

ggplot(data = resultSum, aes(x = sampleSize, y = power)) +
  geom_point(aes(col = p1)) +
  geom_hline(yintercept = c(5, 80), linetype = 2) + 
  facet_wrap(comparator~model, ncol = 2) +
  labs(x = 'Sample Size', y = 'Power (%)', col = 'Colistin mortality')
ggsave('Power-SampleSize.tiff', dpi = 300, width = 15, height = 12)