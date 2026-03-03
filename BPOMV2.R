#BPOMV2
# Load necessary libraries
rm(list = ls())
library(dplyr)
library(readxl)



# 1. Load and Clean Data
file_path <- "D:/NUS Dropbox/Xiangyuan Huang/github/BPOM"
setwd(file_path)
df <- read_excel('Antibiotic list_Indo.xlsx', sheet = 1)

# Categorize Intervention Arms and Filter AKI
df_clean <- df %>%
  rename(intervention = 'Intervention arms',
         intervention = 'Intervention arms',
         aki = 'Acute kidney injury') %>%
  mutate(arm = str_split(intervention, "\\r?\\n")) %>%
  unnest(arm) %>%
  mutate(arm = str_remove(arm, "^\\d+\\.\\s*"),
         arm = str_trim(arm),
         group = ifelse(str_detect(arm, 'Colistin'), 'A',
                        ifelse(str_detect(arm, 'Ceftazidime'),
                               'B', 'C'))) %>%
  filter(aki != "Yes")

# 2. Simulation Parameters
sample_sizes <- seq(400, 600, by = 50)
n_sim <- 1000
alpha <- 0.05

# Probabilities based on your requirements
# Pathogen Distribution
pathogen_probs <- c(CRAB = 0.50, CRE_DA = 0.15, CRE_B = 0.15, CRPA = 0.20)
# Mortality Rates (Outcome)
mortality_rates <- c(A = 0.38, B = 0.40, C = 0.40)

# 3. Simulation Function
run_simulation <- function(N) {
  significant_results <- 0
  
  for (i in 1:n_sim) {
    # Simulate groups based on your proportions (simplified for power calc)
    # Assigning Group A, B, and C randomly based on the clean data distribution
    sim_data <- data.frame(
      Group = sample(c("A", "B", "C"), N, replace = TRUE),
      Pathogen = sample(names(pathogen_probs), N, replace = TRUE, prob = pathogen_probs)
    )
    
    # Simulate Mortality Outcome based on Group
    sim_data$Mortality <- sapply(sim_data$Group, function(g) {
      rbinom(1, 1, mortality_rates[g])
    })
    
    # Generalized Logistic Model (Logit link)
    model <- glm(Mortality ~ Group + Pathogen, data = sim_data, family = binomial)
    
    # Check if Group A (Intervention) is statistically significant
    stats <- summary(model)$coefficients
    if ("GroupB" %in% rownames(stats) && stats["GroupB", "Pr(>|z|)"] < alpha) {
      significant_results <- significant_results + 1
    }
  }
  return(significant_results / n_sim)
}

# 4. Execute and Show Results
results <- data.frame(SampleSize = sample_sizes, Power = NA)

for (i in 1:nrow(results)) {
  results$Power[i] <- run_simulation(results$SampleSize[i])
}

print(results)
