# post process V5
# Preparation ----
rm(list = ls())

library(magrittr)
library(stringr)
library(tidyr)
library(tidyverse)


threshold <- c(0.99, 0.99, 0.95)   
# Defining thresholds of significance for interim analysis 1, 2 and 3


# Environment ----
env <- ifelse(str_detect(getwd(), 'C|D:/'), 1, 2)
t0 <- Sys.time()
set.seed(123)

if (env == 1){
  setwd('D:/NUS Dropbox/Xiangyuan Huang/github/BPOM')
} else if(env == 2){
  setwd('BPOM')
}



# files identification ----
posterior <- dir('results', pattern = 'run', full.names = T)
posterior <- posterior[order(posterior)]
posterior <- posterior[length(posterior)]
# posterior <- dir(, full.names = T)
posterior <- dir(posterior, full.names = T)
posterior <- posterior[str_detect(posterior, 'RData')]

load(posterior[str_detect(posterior, 'cenario.RData')])
load(posterior[str_detect(posterior, 'final.RData')])

posterior <- posterior[str_detect(posterior, 'pos-')]



# posterior loading ----
dfPass <- data.frame()
for (i in 1:length(posterior)){
  cat('Processing:', i, '\n')
  file <- posterior[i]
  scenarioID <- str_extract(file, 'scenario[0-9]+-')
  model <- str_extract(file, '-m[0-9]+-')
  interim <- str_extract(file, 'intrim[0-9]+-')
  iteration <- str_extract(file, 'iteration[0-9]+')
  load(posterior[i])
  if(exists('rd_vector1')){
    rd_vector <- rd_vector1
  } else {rd_vector <- rd_vector2}
  
  passRate <- sum(rd_vector > -0.1)/length(rd_vector)  
  # rate of posterior falling above -0.1
  # BAT mortality is not lower than polymyxin by 10 percent
  
  dfTemp <- data.frame(scenarioID = scenarioID,
                       model = model,
                       interim = interim,
                       iteration = iteration,
                       passRate = passRate)
  dfPass <- rbind(dfPass, dfTemp)
  rm(list = ls(pattern = 'rd_vector'))
}

dfPassP <- dfPass %>% mutate(
  across(c(interim, iteration), parse_number),
  model = str_replace_all(model, '-', ''),
  scenarioID = str_replace(scenarioID, '-', ''),
  threshold = ifelse(interim == 1, threshold[1],
                     ifelse(interim == 2, threshold[2],
                            ifelse(interim == 3, threshold[3], NA))),
  pass = passRate > threshold) %>%
  group_by(scenarioID, model, iteration) %>% 
  summarise(count = n(), passTotal = sum(pass)) %>% 
  filter(count == 3) %>% 
  group_by(scenarioID, model) %>%
  mutate(N_iteration = n(),
         power = 100*sum(passTotal > 0)/N_iteration) %>%
  left_join(scenario, by = 'scenarioID')

ggplot(data = dfPassP, aes(x = sampleSize, y = power)) +
  geom_point(aes(shape = model)) +
  geom_hline(yintercept = c(5, 80), linetype = 2) +
  facet_wrap(~p1+p2, ncol = 2)
  

