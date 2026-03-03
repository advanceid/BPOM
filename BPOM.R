# BPOM
# Preparation ----
rm(list = ls())

library(arm)
library(lme4)
library(magrittr)
library(marginaleffects)
library(parallel)
library(parallelly)
library(readxl)
library(rstanarm)
library(tidyverse)


# Environment ----
env <- ifelse(str_detect(getwd(), 'C:/'), 1,
              2)
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
                         ifelse(str_detect(arm, 'Ceftazidime'),
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
  maxit = 1000
  } else if (env == 2){
    n_cores <- 100
    N_iteration <- 2000
    sampleSize <- seq(200, 1200, 20)
    samplePrior <- c(300)
    p1 <- 0.4
    p2 <- c(0.5, 0.45, 0.4, 0.38, 0.37, 0.36, 0.35)
    p3 <- 0.4
    maxit = 1000
    }


alphaV <- 0.05
zVal <- qnorm(1 - alphaV)

scenario <- expand.grid(sampleSize = sampleSize,
                        samplePrior = samplePrior,
                        p1 = p1,
                        p2 = p2, 
                        p3 = p3) %>%
  mutate(scenarioID = paste0('scenario', row_number()))



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
simu <- function(scenario){
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
    dfRCT %<>% filter(groupT %in% c('A', 'B'))
    dfRCTPrior <- pattern[rep(1:nrow(pattern), pattern$allocationPrior), ]
    dfRCTPrior$death <- rbinom(n = nrow(dfRCTPrior), size = 1, prob = dfRCTPrior$mortality)
    dfRCTPrior %<>% filter(groupT %in% c('A', 'B'))
    
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
    

    model1 <- myTryCatch(bayesglm(
      death ~ groupT + Organism, 
      data = dfRCT, 
      maxit = maxit,
      prior.mean = model_coef[2:5],
      prior.scale = 1,
      prior.mean.for.intercept = model_coef[1],
      prior.scale.for.intercept = 1,
      family = binomial(link = "identity")
    ))
    res1 <- model1[[1]] %>% summary %>% coefficients %>% as.data.frame %>%
      rownames_to_column("term") %>%
      mutate(UL = Estimate + zVal*`Std. Error`) %>% mutate(model = 'bayesglm') %>%
      dplyr::select(term, Estimate, UL, model)
    
    # model2 <- myTryCatch(stan_glm(death ~ groupT + Organism, 
    #                               family = binomial(link = "logit"), 
    #                               data = dfRCT))
    # rd_results <- avg_comparisons(
    #   model2[[1]], conf_level = (1 - 2*alphaV),
    #   variables = c("groupT", 'Organism'))
    # res2 <- rd_results %>% as.data.frame %>%
    #   rename(UL = conf.high, Estimate = estimate) %>%
    #   mutate(model = 'stanglm') %>% 
    #   dplyr::select(term, Estimate, UL, model)
    
    model <- 
      # rbind(res1, res2) %>% 
      res1 %>%
      mutate(scenarioID = scenario$scenarioID[i]) %>%
      mutate(A = sum(dfRCT$groupT == 'A'),
             B = sum(dfRCT$groupT == 'B'),
             C = sum(dfRCT$groupT == 'C'))
    assign(scenario$scenarioID[i], model)
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
  simu(scenario)}, 
  mc.cores = n_cores)
t1 <- Sys.time()
print(t1 - t0)

# Post processing ----
final_results <- do.call(rbind, results)
save(final_results, file = 'final.RData')
load('final.RData')
print(t1 - t0)


power <- final_results %>%
  mutate(term = ifelse(term == 'groupT', 'groupTB', term),
         Estimate = as.numeric(Estimate),
         UL = as.numeric(UL),
         A = as.numeric(A),
         B = as.numeric(B)) %>% 
  filter(term == 'groupTB', !is.na(Estimate)) %>%
  group_by(scenarioID, model) %>%
  summarise(nIteration = n(),
            pos = sum(UL < 0.1),
            ratePos = pos/nIteration,
            mEs = mean(Estimate),
            # merr = mean(sde),
            nA = mean(A),
            nB = mean(B),
            nEff = nA + nB) %>%
  left_join(scenario, by = 'scenarioID') %>%
  mutate(
    marginColiVBAT = p2 - p1,
    marginColiVBAT = as.character(marginColiVBAT),
    p2 = as.character(p2))

write.csv(power, 'power.csv')
power <- read.csv('power.csv')

power %<>% mutate(marginColiVBAT = as.character(100*marginColiVBAT))
ggplot(data = power, aes(x = sampleSize, y = ratePos, group = marginColiVBAT)) +
  geom_point(aes(col = marginColiVBAT)) +
  geom_smooth(aes(col = marginColiVBAT)) +
  geom_hline(yintercept = c(0.05, 0.80), linetype = 2) +
  scale_x_continuous(breaks = seq(200, 1200, 100)) +
  labs(col = 'Coli minus BAT(%)', y = 'Power', x = 'Sample Size')
ggsave('Power-samplesize.tiff', dpi = 300, width = 12, height = 8)

ggplot(data = power, aes(x = nEff, y = ratePos, group = marginColiVBAT)) +
  geom_point(aes(col = marginColiVBAT)) +
  geom_smooth(aes(col = marginColiVBAT)) +
  geom_hline(yintercept = c(0.05, 0.80), linetype = 2) +
  labs(x = 'True sample size') +
  facet_wrap(~model)
ggsave('Power-samplesizeTRUE.tiff', dpi = 300, width = 12, height = 8)



upperBound <- final_results %>% rename(
  estimate = Estimate,
  sde = 'Std. Error') %>% mutate(
    estimate = as.numeric(estimate),
    sde = as.numeric(sde),
    A = as.numeric(A),
    B = as.numeric(B),
    
    ll = estimate - zVal*sde,
    ul = estimate + zVal*sde)

ggplot(data = upperBound, aes(x = sampleSize, y = ul)) +
  geom_point()
