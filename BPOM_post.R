# BPOM post process
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

# Interim analysis O-Brien Fleming
design <- gsDesign(k = 3, test.type = 1, alpha = 0.05, beta = 0.2,
                   timing = c(0.33, 0.66, 1.0), sfu = 'OF', sfupar = 0)
z_Val <- -design$upper$bound
cat('Z values:', z_Val, '\n')
p_nomi <- 1 - 2*pnorm(z_Val)
cat('P values:', p_nomi, '\n')


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
save(resultSum, file = 'results/resultSum.RData')
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