
library(dplyr)
library(tidyr)
library(stringr)
library(readxl)

raw <- read_excel("unfiltered_Antibiotic list_24Feb2026.xlsx")
colnames(raw) <- c("country", "organism", "gene", "aki", "infection_syndrome",
                   "infection_source", "arms")
raw <- raw %>%
  filter(country %in% c("China", "Malaysia", "Singapore", "Thailand", "Indonesia"))

# Classify an individual arm string into a treatment group
classify_arm <- function(arm) {
  a <- tolower(trimws(arm))
  a <- str_remove(a, "^\\d+\\.\\s*")
  if (grepl("colistin|polymyxin", a))    return("Polymyxin B")
  if (grepl("ceftazidime-avibactam", a)) return("Ceftazidime-Avibactam")
  return("Other")
}

classify_organism <- function(organism, gene) {
  org <- toupper(trimws(organism))
  g   <- toupper(trimws(gene))
  if (org == "CRAB") return("CRAB")
  if (org == "CRPA") return("CRPA")
  if (org == "CRE") {
    has_B  <- grepl("\\bB\\b",    g)
    has_AD <- grepl("\\b[AD]\\b", g)
    if (has_B & !has_AD) return("CRE-B")
    if (has_AD)          return("CRE-A/D")
  }
  return(NA_character_)
}

raw <- raw %>%
  mutate(
    org_group  = mapply(classify_organism, organism, gene),
    stratum_id = row_number()
  ) %>%
  filter(!is.na(org_group))

arms_expanded <- raw %>%
  separate_rows(arms, sep = "\n") %>%
  mutate(
    arm      = str_remove(trimws(arms), "^\\d+\\.\\s*"),  # strip numbering
    tx_group = sapply(arm, classify_arm)
  ) %>%
  filter(arm != "")

stratum_fracs <- arms_expanded %>%
  group_by(stratum_id, country, org_group, aki) %>%
  summarise(
    n_total    = n(),
    n_poly     = sum(tx_group == "Polymyxin B"),
    n_caz      = sum(tx_group == "Ceftazidime-Avibactam"),
    n_other    = sum(tx_group == "Other"),
    has_poly   = n_poly  > 0,
    has_caz    = n_caz   > 0,
    .groups    = "drop"
  ) %>%
  mutate(
    frac_poly  = n_poly  / n_total,
    frac_caz   = n_caz   / n_total,
    frac_other = n_other / n_total
  )

avg_fracs <- stratum_fracs %>%
  group_by(country, org_group, aki) %>%
  summarise(
    frac_poly  = mean(frac_poly),
    frac_caz   = mean(frac_caz),
    frac_other = mean(frac_other),
    has_poly   = any(has_poly),     # TRUE if ANY sub-stratum offers a poly arm
    has_caz    = any(has_caz),      # TRUE if ANY sub-stratum offers a caz arm
    .groups    = "drop"
  )

######################### weights ###################################
country_fracs <- tibble(
  country      = c("China", "Malaysia", "Singapore", "Thailand", "Indonesia"),
  country_frac = rep(1/5, 5)
)

org_prev <- tibble(
  org_group = c("CRAB",  "CRE-A/D", "CRE-B", "CRPA"),
  org_prev  = c(0.50,     0.15,      0.15,    0.20)
)

aki_weights <- tibble(
  aki    = c("Yes", "No", "-"),
  aki_wt = c(0.5,   0.5,  1.0)
)

weighted <- avg_fracs %>%
  left_join(country_fracs, by = "country") %>%
  left_join(org_prev,      by = "org_group") %>%
  left_join(aki_weights,   by = "aki") %>%
  mutate(joint_wt = country_frac * org_prev * aki_wt)

N <- 1000

allocation <- weighted %>%
  pivot_longer(
    cols      = c(frac_poly, frac_caz, frac_other),
    names_to  = "frac_type",
    values_to = "arm_frac"
  ) %>%
  mutate(
    tx_group = recode(frac_type,
                      "frac_poly"  = "Polymyxin B",
                      "frac_caz"   = "Ceftazidime-Avibactam",
                      "frac_other" = "Other")
  ) %>%
  filter(arm_frac > 0) %>%
  mutate(raw_wt = joint_wt * arm_frac) %>%
  group_by(country, org_group, tx_group) %>%
  summarise(raw_wt = sum(raw_wt), .groups = "drop") %>%
  mutate(
    proportion = raw_wt / sum(raw_wt),
    n_expected = proportion * N
  ) %>%
  arrange(country, org_group, tx_group)

allocation %>%
  transmute(Country    = country,
            Organism   = org_group,
            `Tx Group` = tx_group,
            `%`        = sprintf("%.1f%%", proportion * 100),
            `Expected N` = round(n_expected, 1)) %>%
  print(n = Inf)

allocation %>%
  group_by(country) %>%
  summarise(`%` = sprintf("%.1f%%", sum(proportion)*100),
            `Expected N` = round(sum(n_expected), 1), .groups="drop") %>%
  arrange(desc(`Expected N`)) %>% rename(Country = country) %>% print()

allocation %>%
  group_by(org_group) %>%
  summarise(`%` = sprintf("%.1f%%", sum(proportion)*100),
            `Expected N` = round(sum(n_expected), 1), .groups="drop") %>%
  arrange(desc(`Expected N`)) %>% rename(Organism = org_group) %>% print()

allocation %>%
  group_by(tx_group) %>%
  summarise(`%` = sprintf("%.1f%%", sum(proportion)*100),
            `Expected N` = round(sum(n_expected), 1), .groups="drop") %>%
  arrange(desc(`Expected N`)) %>% rename(`Treatment Group` = tx_group) %>% print()

# Non-AKI strata: 
non_aki_strata <- weighted %>%
  filter(aki %in% c("-", "No"))

# Flag strata as comparable (has BOTH poly and CAZ arms)
comparable_strata <- non_aki_strata %>%
  filter(has_poly & has_caz)

# how many patients of total asian trial eligible 
comparable_alloc <- comparable_strata %>%
  pivot_longer(
    cols      = c(frac_poly, frac_caz, frac_other),
    names_to  = "frac_type",
    values_to = "arm_frac"
  ) %>%
  mutate(
    tx_group = recode(frac_type,
                      "frac_poly"  = "Polymyxin B",
                      "frac_caz"   = "Ceftazidime-Avibactam",
                      "frac_other" = "Other")
  ) %>%
  filter(arm_frac > 0) %>%
  mutate(raw_wt = joint_wt * arm_frac) %>%
  group_by(country, org_group, tx_group) %>%
  summarise(raw_wt = sum(raw_wt), .groups = "drop")

# Total weight of non-AKI comparable pool (normalised against the full trial)
# Non-AKI patients = 50% of N; comparable strata are a subset of non-AKI strata.
# We first find what fraction of the FULL trial weight the comparable non-AKI
# strata represent, then scale to N/2 non-AKI patients.

# Total raw weight across all strata (from the full allocation)
total_raw_all <- sum(allocation$raw_wt)   # denominator for full trial proportions

# Raw weight of ALL non-AKI strata (comparable + excluded)
total_raw_non_aki <- sum(non_aki_strata$joint_wt *
                           (non_aki_strata$frac_poly +
                              non_aki_strata$frac_caz  +
                              non_aki_strata$frac_other))

# Raw weight of comparable non-AKI strata only
total_raw_comparable <- sum(comparable_alloc$raw_wt)

# Scale: comparable non-AKI patients out of N/2 non-AKI patients
N_non_aki    <- N * 0.5   # = 500
N_comparable <- N_non_aki * (total_raw_comparable / total_raw_non_aki)

comparable_alloc <- comparable_alloc %>%
  mutate(
    prop_within_comparable = raw_wt / total_raw_comparable,
    n_comparable           = prop_within_comparable * N_comparable
  ) %>%
  arrange(country, org_group, tx_group)

cat("\nComparable non-AKI allocation table (country × organism × tx group):\n")
comparable_alloc %>%
  transmute(Country    = country,
            Organism   = org_group,
            `Tx Group` = tx_group,
            `% within comparable` = sprintf("%.1f%%", prop_within_comparable * 100),
            `Expected N` = round(n_comparable, 1)) %>%
  print(n = Inf)

cat("\n── Summary: Polymyxin B vs Ceftazidime-Avibactam (non-AKI, comparable strata) ──\n")
comparable_alloc %>%
  #filter(tx_group != "Other") %>%
  group_by(tx_group) %>%
  summarise(
    `% within comparable` = sprintf("%.1f%%", sum(prop_within_comparable) * 100),
    `Expected N`          = round(sum(n_comparable), 1),
    .groups = "drop"
  ) %>%
  rename(`Treatment Group` = tx_group) %>%
  print()


