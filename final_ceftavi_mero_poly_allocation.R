library(dplyr)
library(tidyr)
library(stringr)
library(readxl)

raw <- read_excel("unfiltered_Antibiotic list_24Feb2026.xlsx")

colnames(raw) <- c("country", "organism", "gene", "aki", "infection_syndrome",
                   "infection_source", "arms")

raw <- raw %>%
  filter(country %in% c("China", "Malaysia", "Singapore", "Thailand", "Indonesia"))

classify_arm <- function(arm) {
  a <- tolower(trimws(arm))
  # Strip leading "1. / 2. / ..." numbering if present
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

# ── FULL ALLOCATION RESULTS ───────────────────────────────────────────────────

cat("════════════════════════════════════════════════════════════════\n")
cat(sprintf("  FULL ALLOCATION TABLE  (N = %d)\n", N))
cat("════════════════════════════════════════════════════════════════\n")
allocation %>%
  transmute(Country    = country,
            Organism   = org_group,
            `Tx Group` = tx_group,
            `%`        = sprintf("%.1f%%", proportion * 100),
            `Expected N` = round(n_expected, 1)) %>%
  print(n = Inf)

cat("\n── By Country ──────────────────────────────────────────────────\n")
allocation %>%
  group_by(country) %>%
  summarise(`%` = sprintf("%.1f%%", sum(proportion)*100),
            `Expected N` = round(sum(n_expected), 1), .groups="drop") %>%
  arrange(desc(`Expected N`)) %>% rename(Country = country) %>% print()

cat("\n── By Organism ─────────────────────────────────────────────────\n")
allocation %>%
  group_by(org_group) %>%
  summarise(`%` = sprintf("%.1f%%", sum(proportion)*100),
            `Expected N` = round(sum(n_expected), 1), .groups="drop") %>%
  arrange(desc(`Expected N`)) %>% rename(Organism = org_group) %>% print()

cat("\n── By Treatment Group ──────────────────────────────────────────\n")
allocation %>%
  group_by(tx_group) %>%
  summarise(`%` = sprintf("%.1f%%", sum(proportion)*100),
            `Expected N` = round(sum(n_expected), 1), .groups="drop") %>%
  arrange(desc(`Expected N`)) %>% rename(`Treatment Group` = tx_group) %>% print()

cat("\n── AKI Split (global, 50/50) ───────────────────────────────────\n")
cat("  AKI patients     :", round(N * 0.5), "\n")
cat("  Non-AKI patients :", round(N * 0.5), "\n")

cat("\n── Verification ────────────────────────────────────────────────\n")
cat("  Sum proportions:", round(sum(allocation$proportion), 6), "(should be 1)\n")
cat("  Sum expected N :", round(sum(allocation$n_expected), 1), "\n")

# Non-AKI strata: AKI="-" or AKI="No"
non_aki_strata <- weighted %>%
  filter(aki %in% c("-", "No"))

# Flag strata as comparable (has BOTH poly and CAZ arms)
comparable_strata <- non_aki_strata %>%
  filter(has_poly & has_caz)

excluded_strata <- non_aki_strata %>%
  filter(!(has_poly & has_caz))

cat("\nExcluded non-AKI strata (do NOT have both Poly + CAZ arms):\n")
excluded_strata %>%
  select(country, org_group, aki, has_poly, has_caz,
         frac_poly, frac_caz, frac_other) %>%
  print(n = Inf)

# Compute allocation within comparable non-AKI strata only
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
  filter(tx_group != "Other") %>%
  group_by(tx_group) %>%
  summarise(
    `% within comparable` = sprintf("%.1f%%", sum(prop_within_comparable) * 100),
    `Expected N`          = round(sum(n_comparable), 1),
    .groups = "drop"
  ) %>%
  rename(`Treatment Group` = tx_group) %>%
  print()

cat("\n── Other arms within comparable strata (not used in comparison) ──\n")
comparable_alloc %>%
  filter(tx_group == "Other") %>%
  group_by(tx_group) %>%
  summarise(
    `% within comparable` = sprintf("%.1f%%", sum(prop_within_comparable) * 100),
    `Expected N`          = round(sum(n_comparable), 1),
    .groups = "drop"
  ) %>%
  rename(`Treatment Group` = tx_group) %>%
  print()

cat(sprintf(
  "\nTotal non-AKI patients              : %d (= N × 50%%)\n", round(N_non_aki)))
cat(sprintf(
  "  In comparable (Poly+CAZ) strata   : %.1f\n", N_comparable))
cat(sprintf(
  "  In excluded strata (Poly-only)    : %.1f\n", N_non_aki - N_comparable))
cat(sprintf(
  "    → of which: Polymyxin B        : %.1f\n",
  N_non_aki - N_comparable))   # Thailand CRAB: all arms are Polymyxin

cat("\n── By Organism within comparable non-AKI strata ───────────────\n")
comparable_alloc %>%
  group_by(org_group) %>%
  summarise(`Expected N` = round(sum(n_comparable), 1), .groups = "drop") %>%
  rename(Organism = org_group) %>%
  arrange(desc(`Expected N`)) %>%
  print()

cat("\n── By Country within comparable non-AKI strata ────────────────\n")
comparable_alloc %>%
  group_by(country) %>%
  summarise(`Expected N` = round(sum(n_comparable), 1), .groups = "drop") %>%
  rename(Country = country) %>%
  arrange(desc(`Expected N`)) %>%
  print()

# ════════════════════════════════════════════════════════════════════════════════
# PART 3: DATAFRAME 1 — % OF TOTAL TRIAL IN EACH TX GROUP
#         (eligible = non-AKI AND in a comparable stratum with both Poly + CAZ)
#
# Denominator = all N=1000 trial participants (full trial weight).
# Numerator   = weight attributable to comparable non-AKI strata, split by tx group.
# "Eligible"  is defined as: AKI="-" or AKI="No", AND stratum offers BOTH a
#             Polymyxin arm AND a Ceftazidime-Avibactam arm.
# ════════════════════════════════════════════════════════════════════════════════

cat("\n")
cat("════════════════════════════════════════════════════════════════\n")
cat("  PART 3: DATAFRAME 1\n")
cat("  % of TOTAL trial (N=1000) by treatment group\n")
cat("  (eligible = non-AKI + comparable stratum)\n")
cat("════════════════════════════════════════════════════════════════\n")

# Full-trial total raw weight (denominator)
full_trial_raw_wt <- weighted %>%
  mutate(row_raw = joint_wt * (frac_poly + frac_caz + frac_other)) %>%
  summarise(total = sum(row_raw)) %>%
  pull(total)

# Raw weight in comparable non-AKI strata, by tx group
comp_raw_by_group <- weighted %>%
  filter(aki %in% c("-", "No"), has_poly, has_caz) %>%
  summarise(
    raw_Polymyxin_B            = sum(joint_wt * frac_poly),
    raw_Ceftazidime_Avibactam  = sum(joint_wt * frac_caz),
    raw_Other                  = sum(joint_wt * frac_other)
  )

# Dataframe 1: one row per tx group, showing % of total trial
df1_pct_of_total <- tibble(
  treatment_group = c("Polymyxin B", "Ceftazidime-Avibactam", "Other"),
  raw_wt = c(
    comp_raw_by_group$raw_Polymyxin_B,
    comp_raw_by_group$raw_Ceftazidime_Avibactam,
    comp_raw_by_group$raw_Other
  )
) %>%
  mutate(
    pct_of_total_trial = raw_wt / full_trial_raw_wt,
    n_of_total_1000    = round(pct_of_total_trial * N, 1)
  ) %>%
  select(
    `Treatment Group`       = treatment_group,
    `% of Total Trial`      = pct_of_total_trial,
    `Expected N (of 1000)`  = n_of_total_1000
  )

# Add a summary row for all eligible (comparable non-AKI) combined
df1_pct_of_total <- bind_rows(
  df1_pct_of_total,
  tibble(
    `Treatment Group`      = "── TOTAL eligible ──",
    `% of Total Trial`     = sum(df1_pct_of_total$`% of Total Trial`),
    `Expected N (of 1000)` = sum(df1_pct_of_total$`Expected N (of 1000)`)
  )
) %>%
  mutate(`% of Total Trial` = sprintf("%.1f%%", `% of Total Trial` * 100))

cat("\ndf1_pct_of_total:\n")
print(df1_pct_of_total, n = Inf)

cat("\nInterpretation: of all 1000 trial participants, the percentages above\n")
cat("are expected to be non-AKI patients in a stratum that offers BOTH a\n")
cat("Polymyxin B arm AND a Ceftazidime-Avibactam arm (comparable strata).\n")
cat("The remaining", sprintf("%.1f%%",
                             (1 - (comp_raw_by_group$raw_Polymyxin_B +
                                     comp_raw_by_group$raw_Ceftazidime_Avibactam +
                                     comp_raw_by_group$raw_Other) / full_trial_raw_wt) * 100),
    "of participants are either AKI patients or in non-comparable strata.\n")

# ════════════════════════════════════════════════════════════════════════════════
# PART 4: DATAFRAME 2 — % OF TOTAL TRIAL BY COUNTRY × ORGANISM × TX GROUP
#         Same eligibility filter as Part 3 (non-AKI, comparable strata).
# ════════════════════════════════════════════════════════════════════════════════

cat("\n")
cat("════════════════════════════════════════════════════════════════\n")
cat("  PART 4: DATAFRAME 2\n")
cat("  % of TOTAL trial (N=1000) by country × organism × tx group\n")
cat("  (eligible = non-AKI + comparable stratum)\n")
cat("════════════════════════════════════════════════════════════════\n")

# Long-form: one row per country × organism × tx group
df2_country_organism <- weighted %>%
  filter(aki %in% c("-", "No"), has_poly, has_caz) %>%
  pivot_longer(
    cols      = c(frac_poly, frac_caz, frac_other),
    names_to  = "frac_type",
    values_to = "arm_frac"
  ) %>%
  mutate(
    treatment_group = recode(frac_type,
                             "frac_poly"  = "Polymyxin B",
                             "frac_caz"   = "Ceftazidime-Avibactam",
                             "frac_other" = "Other")
  ) %>%
  filter(arm_frac > 0) %>%
  group_by(country, org_group, treatment_group) %>%
  summarise(raw_wt = sum(joint_wt * arm_frac), .groups = "drop") %>%
  mutate(
    pct_of_total_trial    = raw_wt / full_trial_raw_wt,
    n_of_total_1000       = round(pct_of_total_trial * N, 1)
  ) %>%
  arrange(country, org_group, treatment_group) %>%
  select(
    Country                = country,
    Organism               = org_group,
    `Treatment Group`      = treatment_group,
    `% of Total Trial`     = pct_of_total_trial,
    `Expected N (of 1000)` = n_of_total_1000
  )

cat("\ndf2_country_organism (long form):\n")
df2_country_organism %>%
  mutate(`% of Total Trial` = sprintf("%.1f%%", `% of Total Trial` * 100)) %>%
  print(n = Inf)

# Wide form: pivot tx groups into columns for easier reading
cat("\ndf2_country_organism (wide form — Expected N per tx group):\n")
df2_country_organism %>%
  select(Country, Organism, `Treatment Group`, `Expected N (of 1000)`) %>%
  pivot_wider(
    names_from  = `Treatment Group`,
    values_from = `Expected N (of 1000)`,
    values_fill = 0
  ) %>%
  mutate(
    `Row Total` = rowSums(across(where(is.numeric)))
  ) %>%
  print(n = Inf)

cat("\ndf2_country_organism (wide form — % of total trial per tx group):\n")
df2_country_organism %>%
  mutate(`% of Total Trial` = sprintf("%.1f%%", `% of Total Trial` * 100)) %>%
  select(Country, Organism, `Treatment Group`, `% of Total Trial`) %>%
  pivot_wider(
    names_from  = `Treatment Group`,
    values_from = `% of Total Trial`,
    values_fill = "0.0%"
  ) %>%
  print(n = Inf)

in_indo <- df2_country_organism[df2_country_organism$Country == "Indonesia",]
