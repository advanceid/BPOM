library(dplyr)
library(tidyr)

setwd( "/Users/lynnperez/Downloads/")

treatmentCombination <- read_excel('Antibiotic list_Indo.xlsx', sheet = 1)


treatmentCombination %<>% rename(intervention = 'Intervention arms',
                                 aki = 'Acute kidney injury')

eacharm <- treatmentCombination %>%
  separate_rows('Intervention arms', sep = "(?=\\d+\\.)") %>%
  filter('Intervention arms' != "") %>% 
  filter(trimws('Intervention arms') != "")

##remove blank rows 
eacharm$group <- NA 

eacharm <- eacharm %>%
  mutate(group = case_when(
    
    #Polymyxin 
    grepl("polymyxin", eacharm$`Intervention arms`, ignore.case = TRUE)  ~ "Polymyxin",
   
     # Ceftazidime
    grepl("ceftazidime", eacharm$`Intervention arms`, ignore.case = TRUE) ~ "Ceftazidime",
    
    # Meropenem
    grepl("meropenem", eacharm$`Intervention arms`, ignore.case = TRUE) ~ "Meropenem"

    ))


##remove blank rows 

eacharm <- eacharm[!is.na(eacharm$group),]

###note where AKI = TRUE and rule out, cannot compare bc ineligible for polymyxin?

eacharm$group <- ifelse(eacharm$`Acute kidney injury` == "Yes", "Ineligible", eacharm$group)

eacharm$silo <- paste0(eacharm$Organism, "_", eacharm$`Resistance gene`)

ratios <- table(eacharm$silo, eacharm$group)

ratios$group_prev <- c(0.5, 0.15, 0.15, 0.2)

##for now, ignore AKI? unsure what the prevalence would be

##For primary analysis: polymyxin versus meropenem:

silo_data <- data.frame(
  silo = c("CRAB", "CRE_Class B", "CRE_Class D_A", "CRPA"), 
  Ceftazidime = c(1, 1, 2, 2), 
  Meropenem = c(1, 1, 1, 1),  
  Polymyxin = c(3, 2, 2, 2)    
)

# CRAB (0.5), CRE_A/D (0.15), CRE_B (0.15), CRPA (0.2)
prevalences <- c(0.5, 0.15, 0.15, 0.2)

#weight probs 
prevs <- silo_data %>%
  mutate(
    prevalence = prevalences,
    # Total available arms per silo
    total_arms = Ceftazidime + Meropenem + Polymyxin,
    # Probability of being assigned to a specific group given the silo
    p_poly = (Polymyxin / total_arms) * prevalence,
    p_mero = (Meropenem / total_arms) * prevalence,
    p_ceft = (Ceftazidime / total_arms) * prevalence
  )

total_p_poly <- sum(prevs$p_poly)
total_p_mero <- sum(prevs$p_mero)
ratio_val <- total_p_poly / total_p_mero

  
  # Parameters for sample size calculation 
p_exp <- 0.38      # Experimental arm mortality (38%)
p_bat <- 0.40      # BAT (40%)
margin <- 0.10     # noninf margin (10%, where higher is worse)
ratio <- ratio_val      # allocation ratio 
alpha <- 0.025     #  one sided 
power <- 0.8 #0.95      # power? this is what was used in prev simulations 

z_alpha <- qnorm(1 - alpha) 
z_beta <- qnorm(power)    

numerator <- (z_alpha + z_beta)^2 * ((p_exp * (1 - p_exp) / ratio) + (p_bat * (1 - p_bat)))
denominator <- (p_exp - p_bat - margin)^2

n_c <- ceiling(numerator / denominator)
n_e <- ceiling(n_c * ratio)

total_ss <- sum(n_c, n_e)

final_ss <- ceiling(total_ss / 0.9)

####now check with ceftazidime arm 

ratio_val <- total_p_poly / total_p_cef

numerator <- (z_alpha + z_beta)^2 * ((p_exp * (1 - p_exp) / ratio) + (p_bat * (1 - p_bat)))
denominator <- (p_exp - p_bat - margin)^2

n_c <- ceiling(numerator / denominator)
n_e <- ceiling(n_c * ratio)

total_ss <- sum(n_c, n_e)

final_ss <- ceiling(total_ss / 0.9) ##assume 10% loss - to -followup 







