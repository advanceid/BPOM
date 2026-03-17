

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
         groupT = ifelse(str_detect(arm, 'Colistin'), 'A', ##NOTE:
                         ifelse(str_detect(arm, 'Meropenem'),
                                'B', 'C')),
         prob = recode(Organism, 
                       'CRAB' = 0.5,
                       'CRE_A' = 0.15,
                       'CRE_B' = 0.15,
                       'CRPA' = 0.2)) %>%
  group_by(Organism) %>%
  mutate(nInt = n())

  n_cores <- 1
  N_iteration <- 20
  sampleSize <- seq(200, 800, 40)
  samplePrior <- c(30, 100, 300)
  p1 <- 0.4
  p2 <- 0.38 #c(0.5, 0.45, 0.4, 0.38, 0.37, 0.36, 0.35) ##dont' need to bother 
  p3 <- 0.4
  maxit <- 1000
  nInterim <- 3
  
  ##tbd: check that pattern allocation is correct 
  


scenario <- expand.grid(sampleSize = sampleSize,
                        samplePrior = samplePrior,
                        p1 = p1,
                        p2 = p2, 
                        p3 = p3) %>%
  mutate(scenarioID = paste0('scenario', row_number()))
#save(scenario, file = 'scenario.RData')

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

##modified to be frequentist 

simu <- function(scenario, comparator) {
 
  rm(list = ls(pattern = "^[0-9]+-[0-9]+$"), envir = .GlobalEnv)
  
  for (i in 1:nrow(scenario)) {
    sampleSize <- scenario$sampleSize[i]
    
    pattern <- patternS
    pattern$allocation <- apply(rmultinom(sampleSize, size = 1, patternS$prob/patternS$nInt), 1, sum)
    
    pattern %<>% 
      group_by(Organism) %>%
      mutate(errO = rnorm(1, 0, 0.01),
             errO = first(errO)) %>% 
      ungroup() %>%
      mutate(
        mortality = case_when(
          groupT == 'A' ~ scenario$p1[i],
          groupT == 'B' ~ scenario$p2[i],
          TRUE          ~ scenario$p3[i]
        ),
        mortality = mortality + errO
      )
    
    # Create the RCT dataframe
    dfRCT <- pattern[rep(1:nrow(pattern), pattern$allocation), ]
    dfRCT$death <- rbinom(n = nrow(dfRCT), size = 1, prob = dfRCT$mortality)
    
    # Shuffle and create interim slices
    dfRCT %<>% 
      slice_sample(prop = 1) %>%
      mutate(slice = ntile(row_number(), nInterim)) %>%
      filter(groupT %in% c('A', comparator))
    
    # 2. Frequentist Interim Analysis Loop
    for (interim in 1:nInterim) {
      dfInterim <- dfRCT %>% filter(slice <= interim)
      
      # Fitting a standard Frequentist GLM
      # Note: Link = 'identity' for Risk Difference; ensure data allows for convergence
      model_freq <- myTryCatch(
        glm(death ~ groupT + Organism, 
            family = binomial(link = "identity"), 
            data = dfInterim,
            maxit = maxit)
      )
      
      if (!is.null(model_freq$value)) {
        # Extracting coefficients, z-values, and SEs
        res <- model_freq[[1]] %>% 
          summary() %>% 
          coefficients() %>% 
          as.data.frame() %>%
          rownames_to_column("term") %>% 
          rename(zVal = `z value`, std = `Std. Error`, pVal = `Pr(>|z|)`) %>%
          mutate(
            model = 'frequentist_glm',
            interim = interim,
            scenarioID = scenario$scenarioID[i],
            interve = sum(dfInterim$groupT == 'A'),
            control = sum(dfInterim$groupT == comparator)
          )
        
        # Store result for this specific scenario and interim
        assign(paste0(scenario$scenarioID[i], "-", interim), res)
      }
    }
  }
  
 
  resNames <- ls(pattern = "^scenario[0-9]+-[0-9]+$")
  
  if(length(resNames) == 0) return(NULL)
  
  resCoef <- mget(resNames)
  resCoef <- do.call(rbind, resCoef)
  
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

save(results, file = 'final.RData')


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
  labs(x = 'Sample Size', y = 'Power (%)', col = 'Polymyxin mortality')
ggsave('Power-SampleSize.tiff', dpi = 300, width = 15, height = 12)

