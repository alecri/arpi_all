#' Main data management for analyses of ARSi-all combination
#' see additional codes in /additional data management
#' written by AC, 230126

library(tidyverse)
library(labelled)
library(lubridate)
library(readxl)

# load data ----

## locked data ----
derived_data <- list.files(path = "data/locked/derived_data/2023-03-26/",
                           pattern = "_cr.RData$", full.names = T)
invisible(lapply(derived_data, load, .GlobalEnv))
load("data/locked/derived_data/2023-03-26/subjects_hs.RData")

# keep only data for patients' id included in the DSMB analysis for the meeting in 221208
load("data/derived/id_20221125.Rdata")
subjects_cr <- filter(subjects_cr, id %in% id_20221125$id) %>% 
  left_join(select(id_20221125, id, status_probio_20221125 = status_probio), by = "id")


## overall survival ----
os <- read.csv("data/locked/raw_data/2023-03-26/MEB Overall Survival Scheduled Events (Regular)_230621.csv",
               header = TRUE, sep = ";", skip = 2, na.string = c(".b", ".c", ".d", ".e", ".n", ".h"),
               stringsAsFactors = FALSE) %>% 
  mutate(subjectId = toupper(subjectId))


## adverse events ----
sae_selected <- read_excel("data/locked/raw_data/2023-03-26/classified_sae_230326.xlsx")


## additional useful objects ----
load("data/derived/additional_data_cr.RData")


## randomization lists ----
randlist_csv_cr <- list.files(path = "randomization/change_randomization_prob/rand_prob_list/mCRPC",
                              pattern = ".csv$", full.names = T)




# data management ----


## correcting errors/protocol violations ----

subjects_cr <- subjects_cr %>% 
  mutate(
    status_probio = replace(status_probio, patientId == "OLVA4003", "excluded"),
    status_probio_20221125 = replace(status_probio_20221125, patientId == "OLVA4003", "excluded"),
    reason_exclusion = replace(reason_exclusion, patientId == "OLVA4003", "not eligible based on in/exclusion criteria"),
    trt_line = replace(trt_line, patientId == "KS1119", 
                       as.character(as.double(trt_line[patientId == "KS1119"]) - 1)),
    reason_exclusion = replace(reason_exclusion, 
                               patientId %in% c("OLVA4006", "SU2021", "CHUL3701", "KAL1602", "OLVA4016", 
                                                "RY1723", "SG1029", "SG1042", "SG1060", "STAV5101", 
                                                "STAV5102", "UHBA7002", "UZG3015", "UZG3019") & rand_round == 1,
                               "unsuccesful liquid biopsy"),
    last_dte = replace(last_dte, patientId == "AK1214" & rand_round == 2, as.Date("2022-02-25")),
    E1_F19_dateDiscont = replace(E1_F19_dateDiscont, patientId == "AK1214" & rand_round == 2, as.Date("2022-02-25")),
    last_dte = replace(last_dte, patientId == "AK1222" & rand_round == 2, as.Date("2022-02-10")),
    E1_F19_dateDiscont = replace(E1_F19_dateDiscont, patientId == "AK1222" & rand_round == 2, as.Date("2022-02-10")),
    last_dte = replace(last_dte, patientId == "UZG3012" & rand_round == 2, as.Date("2022-01-10")),
    E1_F19_dateDiscont = replace(E1_F19_dateDiscont, patientId == "UZG3012" & rand_round == 2, as.Date("2022-01-10")),
    prog = replace(prog, patientId %in% c("AK1214", "AK1222", "UZG3012") & rand_round == 2, 1),
    time_obs = time_length(difftime(last_dte, E1_F16_dateRand), "months")
  )


## baseline variables ----

baseline_cr <- baseline_cr %>% 
  mutate(
    E1_F29_trt_line = replace(E1_F29_trt_line, patientId == "KS1119", 
                       as.character(as.double(E1_F29_trt_line[patientId == "KS1119"]) - 1)),
    # fixed in smart trial as well
    E1_F4_ecogstatus = replace(E1_F4_ecogstatus, patientId == "OLVA4015" & rand_round == 1, 0),
    E1_F4_ecog_2bin = E1_F4_ecogstatus <= 1,
    # treatment history variables (at study entry)
    E1_F29_trt_line = droplevels(E1_F29_trt_line),
    trt_line_entry = droplevels(E1_F29_trt_line),
    adtmono = as.numeric(
      (!is.na(E1_F29_typeSystemicTherapynmCRPC == "ADT monotherapy ") &
         E1_F29_typeSystemicTherapynmCRPC == "ADT monotherapy") |
        (!is.na(E1_F29_mHSPCsystemicTherapyType) & 
           E1_F29_mHSPCsystemicTherapyType == "ADT monotherapy") |
        (is.na(E1_F29_typeSystemicTherapynmCRPC) & is.na(E1_F29_typeSystemicTherapynmCRPC) &
           !is.na(E1_F29_therapyLocDisease) & E1_F29_therapyLocDisease == "Primary ADT-only") |
        (is.na(E1_F29_typeSystemicTherapynmCRPC) & is.na(E1_F29_typeSystemicTherapynmCRPC) &
           !is.na(E1_F29_typeSystemicTherapyBCR_0) & E1_F29_typeSystemicTherapyBCR_0 == 1)
    ),
    type_systemic_before_mcrp = case_when(
      E1_F29_typeSystemicTherapynmCRPCmHSPC == "ADT monotherapy" ~ "ADT monotherapy",
      E1_F29_typeSystemicTherapynmCRPCmHSPC == "ADT + docetaxel" ~ "ADT + docetaxel",
      E1_F29_typeSystemicTherapynmCRPCmHSPC == "ADT + abiraterone" |
        grepl("amide", E1_F29_typeSystemicTherapynmCRPCmHSPC) ~ "ADT + ARPi",
      E1_F29_typeSystemicTherapynmCRPCmHSPC == "Other" ~ "Other"
    ),
    type_systemic_before_mcrp = replace(type_systemic_before_mcrp,
                                        is.na(type_systemic_before_mcrp) &
                                          ((!is.na(E1_F29_therapyLocDisease) & E1_F29_therapyLocDisease == "Primary ADT-only") |
                                          (!is.na(E1_F29_typeSystemicTherapyBCR_0) & E1_F29_typeSystemicTherapyBCR_0 == 1)),
                                        "ADT monotherapy"),
    type_systemic_before_mcrp = replace(type_systemic_before_mcrp, is.na(type_systemic_before_mcrp), "ADT monotherapy"),
    type_systemic_before_mcrp = factor(type_systemic_before_mcrp, 
                                       levels = c("ADT monotherapy", "ADT + docetaxel",
                                                  "ADT + ARPi", "Other")),
    # variable for treatment history
    received_ARSi = (
      grepl("(abiraterone)|(tamide)", E1_F29_typeSystemicTherapynmCRPCmHSPC) +
        (!is.na(E1_F29_mCRPCsystemicTherapyType_1) & E1_F29_mCRPCsystemicTherapyType_1 == 1) +
        (!is.na(E1_F29_mCRPCsystemicTherapyType_2) & E1_F29_mCRPCsystemicTherapyType_2 == 1)
    ),
    received_ARSi_entry = received_ARSi,
    received_taxane = (
      grepl("taxel", E1_F29_typeSystemicTherapynmCRPCmHSPC) +
        (!is.na(E1_F29_mCRPCsystemicTherapyType_0) & E1_F29_mCRPCsystemicTherapyType_0 == 1) +
        (!is.na(E1_F29_mCRPCsystemicTherapyType_4) & E1_F29_mCRPCsystemicTherapyType_4 == 1)
    ),
    received_taxane_entry = received_taxane
  ) %>% 
  left_join(select(subjects_cr, id, trt_received, therapy_class, all_of(signatures_cr), subtype), by = "id") %>% 
  left_join(select(bio_cr, id, E1_F27_ctDNAlevel), by = "id") %>% 
  mutate(E1_F27_ctDNAlevel = factor(E1_F27_ctDNAlevel, levels = 0:2,
                                    labels = c("Low (<5%)", "Medium (5% - 40%)", "High (≥40%)")))


# variable labels for table 1
## baseline labels ----

# variables (and labels) for table 1
lab_table1_cr <- c(
  # Clinical
  "E1_F1_age_consent_entry" = "Median age at study entry (range) - yr",
  "E1_F4_ecog_2bin" = "ECOG performance‐status score of 0 or 1 at study entry — no. (%)",
  "E1_F3_progr_disease" = "Type of CRPC progression at study entry (PCWG3)",
  "E2_F26_location_metastases" = "Location metastases at study entry — no. (%)",
  "E1_F29_LocOrM1" = "Metastatic disease (M1) at diagnosis — no. (%)",
  "type_systemic_before_mcrp" = "Previous systemic therapy for mHSPC/nmCRPC",
  # "trt_line_entry" = "Treatment line for mCRPC at study entry",
  # "received_ARSi_entry" = "Previous androgen‐receptor–pathway inhibitor at study entry — no. (%)",
  # "received_taxane_entry" = "Previous taxane therapy at study entry — no. (%)",
  "E1_F29_trt_line" = "Treatment line for mCRPC at randomization",
  "received_ARSi" = "Previous androgen‐receptor–pathway inhibitor at randomization — no. (%)",
  "received_taxane" = "Previous taxane therapy at randomization — no. (%)",
  "E1_F7_Hb" = "Median hemoglobin at study entry (range) - g/L",
  "E1_F7_PSA" = "Median prostate-specific antigen level (range) at study entry - ng/ml",
  "E1_F7_LDH" = "Median lactate gehydrogenase at study entry (range) — ukat/L",
  "E1_F7_ALP" = "Median alkaline phosphatase level at study entry (range) - ukat/L",
  "E1_F27_ctDNAlevel_entry" = "ctDNA fraction at study entry — no. (%)",
  "TP53- & AR-_entry" = "AR wildtype and TP53 wildtype",
  "DRD+_entry" = "DRD mutated",
  "TP53+_entry" = "TP53 mutated",
  "TEfus+_entry" = "TMPRSS2:ERG fusion"
)

baseline_cr <- baseline_cr %>% 
  group_by(patientId) %>% 
  mutate(
    E1_F1_age_consent_entry = E1_F1_age_consent[rand_round == 1],
    trt_received_1st = trt_received[rand_round == 1],
    E1_F27_ctDNAlevel_entry = E1_F27_ctDNAlevel[rand_round == 1],
    # keep all the variables for table 1 at study entry for rand_round 2
    across(any_of(signatures_cr), ~ replace(., rand_round == 2, .[rand_round == 1]), .names = "{.col}_entry"),
    across(c(any_of(names(lab_table1_cr)), -E1_F29_trt_line, -received_ARSi, -received_taxane), 
           ~ replace(., rand_round == 2, .[rand_round == 1])),
    # updating received_arsi & received_taxane with info from treatment received in 1st rand_round
    received_taxane = replace(received_taxane, rand_round == 2, received_taxane[rand_round == 1] +
                                trt_received[rand_round == 1] %in% c("Docetaxel", "Cabazitaxel")),
    received_ARSi = replace(received_ARSi, rand_round == 2, received_ARSi[rand_round == 1]  +
                              trt_received[rand_round == 1] %in% c("Abiraterone", "Enzalutamide")),
  ) %>% 
  ungroup()




## overall survival ----

os <- os %>% 
  mutate(
    across(contains("_Dateof"), as.Date),
    E1_F46_Alive = replace(E1_F46_Alive, E1_F46_Withdrawn == 1, 1),
    dead = 1 - E1_F46_Alive,
    E1_F46_DateofCheck = replace(E1_F46_DateofCheck, E1_F46_Withdrawn == 1 & !is.na(E1_F46_Dateofwithdraw), 
                                 E1_F46_Dateofwithdraw[E1_F46_Withdrawn == 1 & !is.na(E1_F46_Dateofwithdraw)]),
    date_check_death = as.Date(if_else(dead == 1, E1_F46_DateofDeath, E1_F46_DateofCheck))
  )


subjects_cr <- left_join(subjects_cr, select(os, subjectId, starts_with("E1_F"),
                                             dead, date_check_death), 
                         by = c("patientId" = "subjectId")) %>% 
  mutate(time_os = time_length(difftime(date_check_death, E1_F16_dateRand), "months"))

if (FALSE) {
  
  # some checks on missing info  
  table(subjects_cr$dead, useNA = "ifany")
  check_os <- filter(subjects_cr, (is.na(dead) | is.na(time_os)) & 
                       (is.na(E1_F46_Withdrawn) | E1_F46_Withdrawn == 0))
  table(check_os$status_probio_20221125)
  
  filter(check_os, status_probio_20221125 == "randomized") %>% 
    select(siteName, patientId, rand_round, E1_F19_discontinue, prog, E1_F19_reasonDiscont, therapy_class) %>% 
    knitr::kable()

}


## adverse events ----

sysorgclass_tab <- table(sae_selected$`System organ class`)
#selected_sysorgclass <- names(sort(sysorgclass_tab[sysorgclass_tab >= 5], decreasing = TRUE))
selected_sysorgclass <- names(sort(sysorgclass_tab[sysorgclass_tab >= 0], decreasing = TRUE))

sae_byid <- split(sae_selected, sae_selected$id) %>% 
  lapply(function(x) {
    data.frame(n_sae = nrow(x), any_sae = 1) %>% 
      bind_cols(
        lapply(selected_sysorgclass, function(ae) {
          any(!is.na(x$`System organ class`) & x$`System organ class` %in% ae)
        }) %>% 
          set_names(nm = selected_sysorgclass) %>% 
          bind_cols()
      )
  }) %>% 
  bind_rows(.id = "id")

arsiall_sae <- subjects_cr %>% 
  filter(id %in% id_20221125$id[id_20221125$status_probio == "randomized"]) %>% 
  filter(therapy_class %in% c("ARSi", "Taxane") | (therapy_received %in% c("ARSi", "Taxane"))) %>% 
  left_join(sae_byid, by = "id") %>% 
  mutate(
    across(c(any_sae, n_sae, all_of(selected_sysorgclass)), ~ replace(., is.na(.), FALSE)),
    n_sae_cat = cut(n_sae, breaks = c(0, 1, 2, 3, 10), labels = c("0", "1", "2", "3+"),
                    include.lowest = TRUE, right = FALSE)
  )




## define therapy class re-challenge variable for sensitivity analysis ----

subjects_cr <- subjects_cr %>% 
  left_join(select(baseline_cr, id, E1_F29_LocOrM1, adtmono, starts_with("E1_F29_mCRPCsystemicTherapyType"), 
                   starts_with("received_")), by = "id") %>% 
  mutate(
    therapyclass_rechallange = as.numeric(
      (E1_F29_mCRPCsystemicTherapyType_1 == 1 & E1_F29_mCRPCsystemicTherapyType_2 == 1) |
        (E1_F29_mCRPCsystemicTherapyType_0 == 1 & E1_F29_mCRPCsystemicTherapyType_4 == 1) |
        (received_taxane > 1) | (received_ARSi > 1) |
        (trt_received %in% c("Abiraterone", "Enzalutamide") & received_ARSi > 0) |
        (trt_received %in% c("Docetaxel", "Cabazitaxel") & received_taxane > 0)
    ),
    # manual correction based on Bram's checks
    therapyclass_rechallange = replace(
      therapyclass_rechallange,
      id %in% c("626a72ce0de50c58582fdee5", "5e6f4fb8e816a10a90874035", 
                "5e4c01f28f515912a4980538", "60422ce253afeb3750190cf1",
                "5e0dfe763b75ff282c41f351", "5f27fa82ec758815acb2494c"), 1)
  )



## competing type of progression ----


# data on discontinuation for NU1304 (rand_round = 2) were wrongly excluded
subjects_cr[subjects_cr$id == "5f324e53d308364834c5afe5", 
            c("psa_prog", "radio_prog", "clinic_prog")] <-
  visits_cr[visits_cr$id == "5f324e53d308364834c5afe5" & !is.na(visits_cr$F23_dateOfVisit) &
              visits_cr$F23_dateOfVisit == "2021-08-09", 
            c("psa_prog", "radio_prog", "clinic_prog")]

# computing dimension which drives the progression
subjects_cr <- subjects_cr %>% 
  mutate(
    # manual correction based on correspondence with Bram on slack (231007)
    psa_prog = replace(psa_prog, patientId == "AK1214" & rand_round == 2, 1),
    radio_prog = replace(radio_prog, patientId == "AK1214" & rand_round == 2, 1),
    clinic_prog = replace(clinic_prog, patientId == "AK1214" & rand_round == 2, 0),
    psa_prog = replace(psa_prog, patientId == "AK1222" & rand_round == 2, 1),
    radio_prog = replace(radio_prog, patientId == "AK1222" & rand_round == 2, 1),
    clinic_prog = replace(clinic_prog, patientId == "AK1222" & rand_round == 2, 0),
    clinic_prog = replace(clinic_prog, patientId == "UZG3023" & rand_round == 2, 1),
    
    type_prog_driven = gsub("NA","0", paste0(psa_prog, radio_prog, clinic_prog)),
    type_prog_driven = replace(type_prog_driven, prog == 0, NA),
    type_prog_driven = factor(type_prog_driven,
                              levels = c("000", "100", "010", "110", "001", "101", "011", "111"),
                              labels = c(
                                "Other", "Only PSA", "Only radiological", "PSA and radiological", "Only clinical",
                                "PSA and clinical", "Radiological and clinical", "All dimensions"
                              )),
    type_prog_driven = fct_expand(type_prog_driven, "Death"),
    type_prog_driven = replace(type_prog_driven, (patientId == "RY1730" & rand_round == 2) |
                                 (patientId == "UZG3010" & rand_round == 1), "Death"),
    across(ends_with("_prog"), ~ factor(., levels = c(0, 1), labels = c("No", "Yes"))),
    across(ends_with("_prog"), ~ fct_na_value_to_level(., level = c("Not evaluated"))),
    num_dim_prog = as.factor((psa_prog == "Yes") + (radio_prog == "Yes") + (clinic_prog == "Yes"))
  )


type_prog <- c("PSA" = "psa_prog", "Radiological" = "radio_prog", "Clinical" = "clinic_prog")
# see more data management in code_backup/00_data_management_230511.R

# change in randomization probabilities ----

dat_rand_prob_cr <- lapply(randlist_csv_cr, function(csv)
  read.table(csv, header = T) %>% 
    rownames_to_column(var = "subtype") %>% 
    mutate(date = as.Date(substr(
      gsub("randomization/change_randomization_prob/rand_prob_list/mCRPC/", "", csv), 1, 10
    )))
) %>% 
  bind_rows() %>% 
  mutate(subtype = fct_inorder(subtype)) %>% 
  pivot_longer(cols = -c("subtype", "date"), values_to = "rand_prob", names_to = "treatment") %>% 
  mutate(
    subtype = factor(subtype, levels = subtype_scheme_cr$subtype),
    treatment = fct_inorder(treatment),
    therapy_class = fct_collapse(treatment, ARSi = c("Abiraterone", "Enzalutamide"),
                                 Taxane = c("Docetaxel", "Cabazitaxel"),
                                 Platinum = "Carboplatin", PARPi = "Niraparib")
  ) %>% 
  left_join(subtype_scheme_cr, by = "subtype")

dat_randprob_therapy_cr <- dat_rand_prob_cr %>% 
  filter(therapy_class %in% c("Control", "ARSi", "Taxane"), date <= as.Date("2022-11-25")) %>% 
  group_by(date, subtype, therapy_class) %>% 
  summarise(rand_prob = sum(rand_prob), .groups = "drop") %>% 
  group_by(date, subtype) %>% 
  mutate(rand_prob = rand_prob/sum(rand_prob)) %>% 
  left_join(subtype_scheme_cr, by = "subtype")


# averaging over observed biomarkers prevalences
prev_subtype_cr <- filter(subjects_cr, status_probio_20221125 == "randomized") %>% 
  count(subtype) %>%
  complete(subtype, fill = list(n = 0)) %>% 
  left_join(subtype_scheme_cr, by = "subtype")

dat_randprob_signtherapy_cr <- lapply(signatures_cr, function(s) {
  filter(dat_randprob_therapy_cr, !!sym(s) == 1) %>% 
    left_join(select(prev_subtype_cr, subtype, n), by = "subtype") %>% 
    group_by(date, subtype) %>% 
    mutate(rand_prob_prev = rand_prob*n) %>% 
    group_by(date, therapy_class) %>% 
    summarise(rand_prob_prev = mean(rand_prob_prev), .groups = "drop") %>%
    group_by(date) %>% 
    mutate(rand_prob_prev = rand_prob_prev/sum(rand_prob_prev))
}) %>% 
  set_names(nm = signatures_cr) %>% 
  bind_rows(.id = "signature") %>% 
  mutate(signature = factor(signature, levels = signatures_cr))



# remove unnecessary objects ----
suppressWarnings(
  rm("i", "randlist_csv_cr", "dat_rand_prob_cr", "id_20221125", "prev_subtype_cr", "check_os",
     "sae_selected", "sysorgclass_tab")
)


# save derived data ----

save(list = ls(), file = "analyses/ARSi-all/derived_data/derived_data.RData")
