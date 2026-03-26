library(dplyr)
library(gsDesign)
library(tidyr)
library(tidyverse)
library(magrittr)
library(readxl)

setwd( "D:/NUS Dropbox/Xiangyuan Huang/github/BPOM")

treatmentCombination <- read_excel('Antibiotic list_Indo.xlsx', sheet = 1)


treatmentCombination %<>% rename(intervention = 'Intervention arms',
                                 aki = 'Acute kidney injury')

eacharm <- treatmentCombination %>%
  separate_rows('intervention', sep = "(?=\\d+\\.)") %>%
  filter('intervention' != "") %>% 
  filter(trimws('intervention') != "")

##remove blank rows 
eacharm$group <- NA 

eacharm <- eacharm %>%
  mutate(group = case_when(
    
    #Polymyxin 
    grepl("polymyxin", eacharm$`intervention`, ignore.case = TRUE)  ~ "Polymyxin",
   
     # Ceftazidime
    grepl("ceftazidime", eacharm$`intervention`, ignore.case = TRUE) ~ "Ceftazidime",
    
    # Meropenem
    grepl("meropenem", eacharm$`intervention`, ignore.case = TRUE) ~ "Meropenem"

    ))


##remove blank rows 

eacharm <- eacharm[!is.na(eacharm$group),]

###note where AKI = TRUE and rule out, cannot compare bc ineligible for polymyxin?

eacharm$group <- ifelse(eacharm$aki == "Yes", "Ineligible", eacharm$group)

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
total_p_ceft <- sum(prevs$p_ceft)


  
  # Parameters for sample size calculation 
p_exp <- 0.36      # Experimental arm mortality (38%)
p_bat <- 0.41      # BAT (40%)
p_ceft <- 0.36     # Cefta superior comparison
margin <- 0.10     # noninf margin (10%, where higher is worse)
alpha <- 0.05    #  one sided 
power <- 0.8 #0.95      # power? this is what was used in prev simulations 

z_alpha <- qnorm(1 - alpha) 
z_beta <- qnorm(power)    

ratio_val <- total_p_poly / total_p_mero # allocation ratio 
numerator <- (z_alpha + z_beta)^2 * ((p_exp * (1 - p_exp) / ratio_val) + (p_bat * (1 - p_bat)))
denominator <- (p_exp - p_bat - margin)^2

n_c <- ceiling(numerator / denominator)
n_e <- ceiling(n_c * ratio_val)

total_ss1 <- sum(n_c, n_e)


####now check with ceftazidime arm 
ratio_val <- total_p_poly / total_p_ceft

numerator <- (z_alpha + z_beta)^2 * ((p_exp * (1 - p_exp) / ratio_val) + (p_bat * (1 - p_bat)))
denominator <- (p_exp - p_bat - margin)^2

n_c <- ceiling(numerator / denominator)
n_e <- ceiling(n_c * ratio_val)

total_ss2 <- sum(n_c, n_e)



### interim analysis 
design <- gsDesign(k = 3, test.type = 1, alpha = 0.05, beta = 0.2,
                   timing = c(0.33, 0.66, 1.0), sfu = 'OF', sfupar = -4)
inflation <- max(design$n.I)


final_ss1 <- total_ss1*inflation/(0.9*(total_p_poly + total_p_mero))
final_ss1
final_ss2 <- total_ss2*inflation/(0.9*(total_p_poly + total_p_ceft))
final_ss2

final_ss1 <- total_ss1/(0.9*(total_p_poly + total_p_mero))
final_ss1

total_ss1/((total_p_poly + total_p_mero))


total_ss1/((total_p_poly + total_p_mero))
