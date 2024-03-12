library(tidyverse)
library(labelled)
library(scales)
library(rstanarm)
library(survival)
library(flexsurv)
library(Epi)
library(survminer)
library(ggsurvfit)
library(ggpp)
library(ggtext)

save_results <- TRUE


# load data ----

load("analyses/ARSi-all/derived_data/derived_data.RData")

## useful objects ----
date <- "2023-03-26"
therapies <- c("Control", "ARSi", "Taxane")
signatures_cr <- set_names(signatures_cr, nm = signatures_cr)

type_analyses <- c("primary" = "primary", "interaction" = "interaction", "os" = "os")
type_sensitivity <- c("main" = "main", "1line" = "1line", "norechal" = "norechal", 
                      "adtmono" = "adtmono")
# if the same as type_analyses, remove
type_sensitivity_os <- c("main" = "main", "1line" = "1line", "norechal" = "norechal", 
                         "adtmono" = "adtmono")

## naming for figures & plots ----
therapies_legend <- c("Control" = "Physician's choice", "ARSi" = "AR Pathway Inhibitors", "Taxane" = "Taxanes")
therapies_legend_long <- c("Control" = "Physician's choice",
                           "ARSi" = "Androgen Receptor \n Pathway Inhibitors", "Taxane" = "Taxanes")
cmp_legend <- c("ARSi-Control" = "AR Pathway Inhibitors vs Physician's choice", 
                "ARSi-Taxane" = "AR Pathway Inhibitors vs Taxanes", 
                "Taxane-Control" = "Taxanes vs Physician's choice")
signatures_legend <- c("All patients", "*AR* wild-type and *TP53* wild-type", 
                       "Homologous Recombination Deficiency", 
                       "*TP53*-altered", "*TMPRSS2-ERG* fusion positive") %>% 
  set_names(nm = signatures_cr)
signatures_legend_complement <- c("*AR*- or *TP53*-altered", 
                                  "Homologous Recombination Proficient", 
                                  "*TP53* wild-type", 
                                  "*TMPRSS2-ERG* fusion negative") %>% 
  set_names(nm = signatures_cr[-1])
sensitivity_legend <- c("main" = "Main", "1line" = "Only first liners", 
                        "norechal" = "No re-challenge")
sensitivity_legend_os <- c("main" = "Main", "1line" = "Only first lines", 
                           "norechal" = "No re-challenge",
                           "adtmono" = "Received ADT monotherapy")


## colors ----
coltherapies <- set_names(RColorBrewer::brewer.pal(n = 3, name = 'Dark2'), nm = therapies)
coltherapies_cmp <- set_names(RColorBrewer::brewer.pal(n = 3, name = 'Set1'), nm = cmp_legend)
col_withinsign <- set_names(RColorBrewer::brewer.pal(n = 3, name = 'Set1')[c(2, 1)], 
                            nm = c("Signature negative", "Signature positive"))

# data management -----

dat_baseline <- baseline_cr  %>% 
  mutate(
    received_ARSi = factor(received_ARSi, levels = 0:3, labels = 
                             c("None", "One regimen", "Two regimens", "More than two regimens")),
    received_ARSi = fct_collapse(received_ARSi, "More than one regimen" = c("Two regimens", "More than two regimens")),
    received_ARSi_entry = factor(received_ARSi_entry, levels = 0:3, labels = 
                                   c("None", "One regimen", "Two regimens", "More than two regimens")),
    received_ARSi_entry = fct_collapse(received_ARSi_entry, "More than one regimen" = c("Two regimens", "More than two regimens")),
    received_taxane = factor(received_taxane, levels = 0:3, labels = 
                               c("None", "One regimen", "Two regimens", "More than two regimens")),
    received_taxane = fct_collapse(received_taxane, "More than one regimen" = c("Two regimens", "More than two regimens")),
    received_taxane_entry = factor(received_taxane_entry, levels = 0:3, labels = 
                                     c("None", "One regimen", "Two regimens", "More than two regimens")),
    received_taxane_entry = fct_collapse(received_taxane_entry, "More than one regimen" = c("Two regimens", "More than two regimens"))
  ) %>% 
  set_variable_labels(.labels = as.list(lab_table1_cr), .strict = FALSE)


## data for analysis ----

# data for analysis of arsi-all graduation
arsi_all <- subjects_cr %>% 
  filter(status_probio_20221125 == "randomized", therapy_class %in% therapies) %>% 
  mutate(therapy_class = droplevels(therapy_class)) 

# table(arsi_all$therapy_received[arsi_all$therapy_class == "Control"], useNA = "ifany")
# prop.table(table(arsi_all$therapy_received[arsi_all$therapy_class == "Control"], useNA = "ifany"))
# table(arsi_all$therapy_class)

# subjects id for different analyses
id_sensitivity <- list(
  "main" = arsi_all$id,
  "1line" = arsi_all$id[arsi_all$trt_line == 1],
  "norechal" = arsi_all$id[arsi_all$therapyclass_rechallange == 0],
  "adtmono" = arsi_all$id[arsi_all$adtmono == 1 & arsi_all$trt_line == 1]
)

id_sensitivity_os <- list(
  "main" = arsi_all$id[arsi_all$rand_round == 1],
  "1line" = arsi_all$id[arsi_all$trt_line == 1 & arsi_all$rand_round == 1 & arsi_all$trt_line == 1],
  "norechal" = arsi_all$id[arsi_all$therapyclass_rechallange == 0 & arsi_all$rand_round == 1],
  "adtmono" = arsi_all$id[arsi_all$adtmono == 1 & arsi_all$rand_round == 1 & arsi_all$trt_line == 1]
  # "metachronous" = arsi_all$id[arsi_all$trt_line == 1 & arsi_all$E1_F29_LocOrM1 == 0],
  # "synchronous" = arsi_all$id[arsi_all$trt_line == 1 & arsi_all$E1_F29_LocOrM1 == 1]
)

# list of data sets for primary analysis
dat_primary <- lapply(id_sensitivity, function(i) {
  lapply(signatures_cr, function(s) {
    filter(arsi_all, id %in% i, !!sym(s) == 1)
  })
})
# list of data sets for interaction analysis
dat_interaction <- lapply(id_sensitivity, function(i) {
  lapply(signatures_cr[-1], function(s) {
    filter(arsi_all, id %in% i, therapy_class != "Control") %>% 
      rename(signature = all_of(s)) %>% 
      mutate(therapy_class = factor(therapy_class, levels = therapies[c(3, 2)]))
  })
})
# list of data sets for overall survival
dat_os <- lapply(id_sensitivity_os, function(i) {
  lapply(signatures_cr, function(s) {
    filter(arsi_all, id %in% i, !!sym(s) == 1)
  })
})

## negative controls ----

negativectrl <- lapply(therapies[-1], function(tx) {
  lapply(signatures_cr, function(s) {
    filter(dat_primary$main[[s]], therapy_received == tx | therapy_class == tx) %>% 
      mutate(therapy_class = droplevels(therapy_class))
  })
}) %>% setNames(nm = therapies[-1])

negativectrl_os <- lapply(therapies[-1], function(tx) {
  lapply(signatures_cr, function(s) {
    filter(dat_os$main[[s]], therapy_received == tx | therapy_class == tx) %>% 
      mutate(therapy_class = droplevels(therapy_class))
  })
}) %>% setNames(nm = therapies[-1])


## data flowchart ----

# info for consort flowchart 
info_flowchart <- lst(
  nsigned_hs = nrow(filter(subjects_hs, E1_F1_consentDate <= as.Date("2022-11-24"))),
  n_signed_cr = nrow(filter(subjects_cr, rand_round == 1)),
  reason_exclusion = table(subjects_cr$reason_exclusion[subjects_cr$status_probio_20221125 == "excluded" &
                                                          subjects_cr$rand_round == 1], useNA = "always"),
  n_randomized = nrow(filter(subjects_cr, rand_round == 1, status_probio_20221125 == "randomized")),
  n_1strand = table(subjects_cr$therapy_class[subjects_cr$status_probio_20221125 == "randomized" & subjects_cr$rand_round == 1], useNA = "always"),
  n_1strand_tx = c(n_1strand[c("Control", "ARSi", "Taxane")], "Other_tx" = sum(n_1strand[c("Platinum", "PARPi")])),
  n_2strand = table(subjects_cr$therapy_class[subjects_cr$status_probio_20221125 == "randomized" & subjects_cr$rand_round == 2], useNA = "always"),
  n_2strand_tx = c(n_2strand[c("Control", "ARSi", "Taxane")], "Other_tx" = sum(n_2strand[c("Platinum", "PARPi")]))
)
# filter(subjects_cr, status_probio_20221125 == "excluded", rand_round == 1,
#        reason_exclusion == "other") %>% select(patientId, reason_exclusion_other)


add_info_flowchart <- list(
  
  "Discontinuation (rand_round = 1)" = subjects_cr %>% 
    mutate(
      group = fct_collapse(therapy_class, "Other" = c("Platinum", "PARPi")),
      E1_F19_discontinue = replace(E1_F19_discontinue, prog == 1, 1)
    ) %>%
    filter(status_probio_20221125 == "randomized", rand_round == 1) %>% 
    group_by(group) %>% 
    summarise(E1_F19_discontinue = sum(E1_F19_discontinue, na.rm = T),
              prog = sum(prog, na.rm = T)),

  "Reason disc not progression (rand_round = 1)" = subjects_cr %>% 
    mutate(group = fct_collapse(therapy_class, "Other" = c("Platinum", "PARPi"))) %>%
    filter(status_probio_20221125 == "randomized", rand_round == 1, 
           E1_F19_discontinue == 1, prog == 0) %>% 
    select(patientId, group, E1_F19_reasonDiscont, E1_F19_reasonOther, status_probio, therapy_class) %>% 
    arrange(group),
  
  "Re-randomization" = subjects_cr %>% 
    mutate(group = fct_collapse(therapy_class, "Other" = c("Platinum", "PARPi"))) %>%
    group_by(patientId) %>% 
    mutate(first_group = therapy_class[rand_round == 1],
           first_group = fct_collapse(first_group, "Other" = c("Platinum", "PARPi"))) %>% 
    filter(status_probio_20221125 == "randomized", rand_round == 2) %>% 
    group_by(first_group) %>% 
    count(),
  
  "Discontinuation (rand_round = 2)" = subjects_cr %>% 
    mutate(group = fct_collapse(therapy_class, "Other" = c("Platinum", "PARPi")),
           E1_F19_discontinue = replace(E1_F19_discontinue, prog == 1, 1)) %>%
    filter(status_probio_20221125 == "randomized", rand_round == 2) %>% 
    group_by(group) %>% 
    summarise(E1_F19_discontinue = sum(E1_F19_discontinue, na.rm = T),
              prog = sum(prog, na.rm = T)),
  
  "Reason disc not progression (rand_round = 2)" = subjects_cr %>% 
    mutate(
      group = fct_collapse(therapy_class, "Other" = c("Platinum", "PARPi"))
    ) %>%
    filter(status_probio_20221125 == "randomized", rand_round == 2, 
           E1_F19_discontinue == 1, prog == 0) %>% 
    select(patientId, group, E1_F19_reasonDiscont, E1_F19_reasonOther, status_probio, therapy_class) %>% 
    arrange(group)
  
)


## data adverse events ----

arsiall_sae <- arsiall_sae %>% 
  mutate(
    therapy_class = factor(therapy_class, levels = names(therapies_legend),
                           labels = therapies_legend),
    therapy_received = factor(therapy_received, levels = names(therapies_legend)[-1],
                              labels = therapies_legend[-1]),
    therapy_group = as.character(therapy_class),
    therapy_group = replace(therapy_group, therapy_group == therapies_legend[1],
                            paste0(therapies_legend[1], " (", 
                                   therapy_received[therapy_group == therapies_legend[1]], ")")),
    therapy_group = factor(therapy_group, 
                           levels = c("Physician's choice (AR Pathway Inhibitors)",
                                      "Physician's choice (Taxanes)",
                                      "AR Pathway Inhibitors", "Taxanes"),
                           labels = c("Physician's choice \n(AR Pathway Inhibitors)",
                                      "Physician's choice \n(Taxanes)",
                                      "AR Pathway Inhibitors", "Taxanes"))
  )




# analysis ----

## primary analysis ----

### tabular results ----

tabres_main <- lapply(type_sensitivity, function(sens) {
  lapply(signatures_cr, function(s) { 
    dat_primary[[sens]][[s]] %>%
      group_by(therapy_class) %>%
      summarize(
        n = n(), n_1strand = sum(rand_round == 1, na.rm = T), 
        n_rerand = sum(rand_round == 2, na.rm = T),
        progr = sum(prog, na.rm = T), progr_1strand = sum(prog[rand_round == 1], na.rm = T),
        progr_rerand = sum(prog[rand_round == 2], na.rm = T),
        PT = sum(time_obs, na.rm = T), .groups = "drop"
      ) %>%
      complete(therapy_class, fill = list(n = 0, n_1strand = 0, n_rerand = 0, PT = 0,
                                         progr = 0, progr_1strand = 0, progr_rerand = 0))
  }) %>% 
    bind_rows(.id = "signature")
}) %>% 
  bind_rows(.id = "type_sensitivity") %>% 
  bind_cols(
    lapply(seq(nrow(.)), function(i){
      tidy(poisson.test(.$progr[i], T = .$PT[i], conf.level = 0.9))
    }) %>% 
      bind_rows() %>% 
      select(rate = estimate, `rate 5%` = conf.low, `rate 95%` = conf.high)
  ) %>% 
  mutate(
    signature = factor(signature, levels = signatures_cr),
    type_sensitivity = fct_inorder(type_sensitivity)
  )


### bayesian models ----

# individual seed for reproducibility
seed_primary <- lapply(subtypes_sign_cr, function(s) {
    as.numeric(as.Date("2023-03-26")) +
      sum(as.integer(factor(s, levels = subtype_scheme_cr$subtype)))
})

fit_main <- lapply(type_sensitivity, function(sens) {
  lapply(signatures_cr, function(s) {
    if (all(table(dat_primary[[sens]][[s]]$therapy_class) > 0)) {
      fit <- stan_surv(Surv(time_obs, prog) ~ therapy_class, 
                       data = dat_primary[[sens]][[s]],
                       basehaz = "weibull-aft", algorithm = "sampling",
                       prior = normal(0, .5), prior_aux = exponential(1),
                       prior_intercept = normal(3, 1),
                       chains = 4, iter = 2000, warmup = 500, 
                       seed = seed_primary[[s]] + 
                         (seq(type_sensitivity) - 1)[type_sensitivity == sens])
      fit
    }
  })
})
# warnings() # visualize warnings


### treatment effects ----

te_main <- lapply(type_sensitivity, function(sens) { 
  lapply(signatures_cr, function(s) {
    fit <- fit_main[[sens]][[s]]
    if (!is.null(fit)) {
      data.frame(
        "logSTR_ARSivsControl" = as.matrix(fit)[, "therapy_classARSi"],
        "logSTR_TaxanevsControl" = as.matrix(fit)[, "therapy_classTaxane"],
        "logSTR_ARSivsTaxane" = as.matrix(fit)[, "therapy_classARSi"] - as.matrix(fit)[, "therapy_classTaxane"],
        "shape" = as.matrix(fit)[, "weibull-shape"]
      )
    }
  }) %>% 
    bind_rows(.id = "signature")
}) %>% 
  bind_rows(.id = "type_sensitivity") %>% 
  pivot_longer(cols = starts_with("logSTR"), names_to = "therapy_class", 
               values_to = "logSTR") %>% 
  separate(therapy_class, into = c("therapy_class", "comparison"), sep = "vs") %>% 
  mutate(
    therapy_class = factor(gsub("logSTR_", "", therapy_class), levels = therapies),
    comparison = factor(comparison, levels = therapies[-2]),
    signature = factor(signature, levels = signatures_cr),
    type_sensitivity = fct_inorder(type_sensitivity),
    logHR = -logSTR*shape
  )

# summary table
tesum_main <- te_main %>% 
  group_by(type_sensitivity, signature, therapy_class, comparison) %>% 
  reframe(
    prob_sup = mean(logSTR >= 0),
    STR = quantile(exp(logSTR), c(0.5, 0.05, 0.95)),
    HR = quantile(exp(logHR), c(0.5, 0.05, 0.95)),
    q = c("median", "low", "high")
  ) %>% 
  pivot_wider(names_from = "q", values_from = c("STR", "HR"))



### median PFS ----
medianPFS_primary <- lapply(type_sensitivity, function(sens) {
  lapply(signatures_cr, function(s) {
    
    fit <- fit_main[[sens]][[s]]
    if (!is.null(fit)) {
      log_scale_dat <- data.frame(
        "Control" = as.matrix(fit)[, c("(Intercept)")],
        "ARSi" = rowSums(as.matrix(fit)[, c("(Intercept)", "therapy_classARSi")]),
        "Taxane" = rowSums(as.matrix(fit)[, c("(Intercept)", "therapy_classTaxane")])
      )
      median_dat <- (log(2)^(1/as.matrix(fit)[, "weibull-shape"])) * exp(log_scale_dat)
      
      median_dat %>% 
        pivot_longer(cols = everything(), names_to = "therapy_class", values_to = "median")
    }
  }) %>%
    bind_rows(.id = "signature")
}) %>%
  bind_rows(.id = "type_sensitivity")

medianPFS_primary_table <- lapply(type_sensitivity, function(sens) {
  lapply(signatures_cr, function(s) {
    filter(medianPFS_primary, type_sensitivity == sens, signature == s) %>% 
      mutate(Group = factor(therapy_class, levels = names(therapies_legend), labels = therapies_legend)) %>% 
      group_by(Group) %>% 
      summarise(
        "Median Time \n (90% Credible Intervals)" = paste0(
          format(round(quantile(median, .5), 1), nsmall = 1), " (", 
          format(round(quantile(median, .1), 1), nsmall = 1), ", ",
          format(round(quantile(median, .95), 1), nsmall = 1), ")")
      )
  })
})

### posterior survival curves ----

km_primary <- list()
postsurv_primary <- lapply(type_sensitivity, function(sens) { 
  lapply(signatures_cr, function(s) {
    
    if (!is.null(fit_main[[sens]][[s]])) {
      # KM curves
      dat <- dat_primary[[sens]][[s]] %>% 
        mutate(therapy_class = factor(therapy_class, levels = names(therapies_legend),
                                      labels = therapies_legend))
      km_primary[[sens]][[s]] <<- survfit2(Surv(time_obs, prog) ~ therapy_class, data = dat, conf.int = .9)
      pkm <- km_primary[[sens]][[s]] %>% 
        ggsurvfit(linewidth = 1) +
        add_risktable(theme = list(theme_risktable_default(axis.text.y.size = 11, plot.title.size = 11),
                                   theme(plot.title = element_text(face = "bold")))) +
        add_censor_mark() +
        #coord_cartesian(xlim = c(0, 27)) + 
        scale_y_continuous(labels = percent, breaks = seq(0, 1, by = .25)) +
        labs(x = "Time from randomization, months", y = "Progression Free Survival (%)", col = "Group") +
        theme_survminer() + theme(legend.position = "top") +
        ggpp::annotate(geom = "table", x = 40 + 5*(s == "all"), y = 1.1 - 0.1*(s == "all"), 
                       label = list(medianPFS_primary_table[[sens]][[s]]), 
                       table.theme = ttheme_gtlight, size = 2.7 + 0.3*(s == "all"))
        #expand_limits(x = 30 - 3*(s == "all"))
      
      # posterior survival curves
      newdat <- data.frame(therapy_class = therapies)
      surv_post <- posterior_survfit(fit_main[[sens]][[s]],
                                     newdata = newdat, extrapolate = TRUE,
                                     draws = 2000, return_matrix = T,
                                     times = 0, prob = .9
      )
      time_grid <- sapply(surv_post, function(x) attr(x, "times")[1])
      surv_dat_post <- lapply(seq_along(surv_post), function(i) {
        data.frame(ndraw = seq(nrow(surv_post[[i]])),
                   time = time_grid[i],
                   surv = unname(surv_post[[i]]))
      }) %>%
        bind_rows() %>% 
        pivot_longer(cols = starts_with("surv"), names_to = "therapy_class", values_to = "surv") %>% 
        mutate(
          therapy_class = factor(therapy_class, levels = paste0("surv.", seq(nrow(newdat))), 
                                 labels = therapies_legend[newdat$therapy_class])
        ) %>% 
        group_by(time, therapy_class) %>% 
        reframe(surv = quantile(surv, c(0.5, 0.05, 0.95)),
                q = c("median", "low", "high")) %>% 
        pivot_wider(names_from = "q", values_from = "surv", names_prefix = "surv_") 
      
      pkm <- pkm +
        geom_line(data = surv_dat_post, aes(time, surv_median, color = therapy_class)) +
        geom_ribbon(data = surv_dat_post, aes(x = time, y = surv_median, ymin = surv_low, 
                                              ymax = surv_high, fill = therapy_class), alpha = .15) +
        scale_color_manual(labels = therapies_legend, values = set_names(coltherapies, therapies_legend)) +
        scale_fill_manual(labels = therapies_legend, values = set_names(coltherapies, therapies_legend)) +
        labs(col = "Group", fill = "Group")
      pkm
      
    }
  })
})



## interaction analysis -----

### tabular results ----

tabres_inter <- lapply(type_sensitivity, function(sens) {
  lapply(signatures_cr[-1], function(s) { 
    dat_interaction[[sens]][[s]] %>%
      group_by(therapy_class, signature) %>%
      summarize(
        n = n(), n_1strand = sum(rand_round == 1, na.rm = T), 
        n_rerand = sum(rand_round == 2, na.rm = T),
        progr = sum(prog, na.rm = T), progr_1strand = sum(prog[rand_round == 1], na.rm = T),
        progr_rerand = sum(prog[rand_round == 2], na.rm = T),
        PT = sum(time_obs, na.rm = T), .groups = "drop"
      ) %>%
      complete(therapy_class, signature, 
               fill = list(n = 0, n_1strand = 0, n_rerand = 0, PT = 0,
                           progr = 0, progr_1strand = 0, progr_rerand = 0)) %>% 
      mutate(sign_bin = as.character(signature))
  }) %>% 
    bind_rows(.id = "signature")
}) %>% 
  bind_rows(.id = "type_sensitivity") %>% 
  bind_cols(
    lapply(seq(nrow(.)), function(i){
      tidy(poisson.test(.$progr[i], T = .$PT[i], conf.level = 0.9))
    }) %>% 
      bind_rows() %>% 
      select(rate = estimate, `rate 5%` = conf.low, `rate 95%` = conf.high)
  ) %>% 
  mutate(
    signature = factor(signature, levels = signatures_cr[-1]),
    type_sensitivity = fct_inorder(type_sensitivity)
  )


### bayesian models ----

fit_inter <- lapply(type_sensitivity, function(sens) {
  lapply(signatures_cr[-1], function(s) {
    
    if (all(table(dat_interaction[[sens]][[s]]$therapy_class, dat_interaction[[sens]][[s]]$signature) > 0)) {
      fit <- stan_surv(Surv(time_obs, prog) ~ therapy_class*signature, 
                       data = dat_interaction[[sens]][[s]],
                       basehaz = "weibull-aft", algorithm = "sampling",
                       prior = normal(0, .5), prior_aux = exponential(1),
                       prior_intercept = normal(3, 1),
                       chains = 4, iter = 2000, warmup = 500, 
                       seed = seed_primary[[s]] + 123 +
                         (seq(type_sensitivity) - 1)[type_sensitivity == sens])
      fit
    }
  })
})
# warnings() # visualize warnings


### treatment effects ----

te_inter <- lapply(type_sensitivity, function(sens) { 
  lapply(signatures_cr[-1], function(s) {
    fit <- fit_inter[[sens]][[s]]
    if (!is.null(fit)) {
      data.frame(
        "logSTR_ARSivsTaxane_0" = as.matrix(fit)[, "therapy_classARSi"],
        "logSTR_ARSivsTaxane_1" =  as.matrix(fit)[, "therapy_classARSi"] + as.matrix(fit)[, "therapy_classARSi:signature"],
        "shape" = as.matrix(fit)[, "weibull-shape"]
      )
    }
  }) %>% 
    bind_rows(.id = "signature")
}) %>% 
  bind_rows(.id = "type_sensitivity") %>% 
  pivot_longer(cols = starts_with("logSTR"), names_to = "therapy_class", 
               values_to = "logSTR") %>% 
  separate(therapy_class, into = c("te", "therapy_class", "sign_bin"), sep = "_") %>% 
  separate(therapy_class, into = c("therapy_class", "comparison"), sep = "vs") %>% 
  mutate(
    signature = factor(signature, levels = signatures_cr[-1]),
    type_sensitivity = fct_inorder(type_sensitivity),
    logHR = -logSTR*shape
  ) %>% 
  select(-te)

# interaction coefficient
coef_inter <- lapply(type_sensitivity, function(sens) { 
  lapply(signatures_cr[-1], function(s) {
    fit <- fit_inter[[sens]][[s]]
    if (!is.null(fit)) data.frame("b3" = as.matrix(fit)[, "therapy_classARSi:signature"])
  }) %>% 
    bind_rows(.id = "signature")
}) %>% 
  bind_rows(.id = "type_sensitivity")



# summary table
tesum_inter <- te_inter %>% 
  group_by(type_sensitivity, signature, therapy_class, comparison, sign_bin) %>% 
  reframe(
    prob_sup = mean(logSTR >= 0),
    STR = quantile(exp(logSTR), c(0.5, 0.05, 0.95)),
    HR = quantile(exp(logHR), c(0.5, 0.05, 0.95)),
    q = c("median", "low", "high")
  ) %>% 
  pivot_wider(names_from = "q", values_from = c("STR", "HR"))


### median PFS ----
medianPFS_inter <- lapply(type_sensitivity, function(sens) {
  lapply(signatures_cr[-1], function(s) {
    
    fit <- fit_inter[[sens]][[s]]
    if (!is.null(fit)) {
      log_scale_dat <- data.frame(
        "ARSi_0" = rowSums(as.matrix(fit)[, c("(Intercept)", "therapy_classARSi")]),
        "ARSi_1" = rowSums(as.matrix(fit)[, c("(Intercept)", "therapy_classARSi", "signature", "therapy_classARSi:signature")]),
        "Taxane_0" = rowSums(as.matrix(fit)[, c("(Intercept)"), drop = FALSE]),
        "Taxane_1" = rowSums(as.matrix(fit)[, c("(Intercept)", "signature")])
      )
      median_dat <- (log(2)^(1/as.matrix(fit)[, "weibull-shape"])) * exp(log_scale_dat)
      
      median_dat %>% 
        pivot_longer(cols = everything()) %>% 
        separate(name, into = c("therapy_class", "sign_bin"), sep = "_")  
    }
  }) %>%
    bind_rows(.id = "signature")
}) %>%
  bind_rows(.id = "type_sensitivity")

medianPFS_inter_table <- lapply(type_sensitivity, function(sens) {
  lapply(signatures_cr[-1], function(s) {
    lapply(c("0" = 0, "1" = 1), function(subs) {
      filter(medianPFS_inter, type_sensitivity == sens, signature == s, sign_bin == subs) %>% 
        mutate(
          Group = factor(therapy_class, levels = names(therapies_legend), labels = therapies_legend)
        ) %>% 
        group_by(Group) %>% 
        summarise(
          "Median Time \n (90% Credible Intervals)" = paste0(
            format(round(quantile(value, .5), 1), nsmall = 1), " (", 
            format(round(quantile(value, .1), 1), nsmall = 1), ", ",
            format(round(quantile(value, .95), 1), nsmall = 1), ")"),
          .groups = "drop"
        )
    })
  })
})


### posterior survival curves ----

km_inter <- list()
postsurv_interaction <- lapply(type_sensitivity, function(sens) { 
  lapply(signatures_cr[-1], function(s) {
    
    if (!is.null(fit_inter[[sens]][[s]])) {
      lapply(c("0" = 0, "1" = 1), function(subs) {
        
        # KM curves
        dat <- filter(dat_interaction[[sens]][[s]], signature == subs) %>% 
          mutate(therapy_class = factor(therapy_class, levels = names(therapies_legend[-1]),
                                        labels = therapies_legend[-1]))
        km_inter[[sens]][[s]][[as.character(subs)]] <<- survfit2(Surv(time_obs, prog) ~ therapy_class, 
                                                                 data = dat, conf.int = .9)
        pkm <- km_inter[[sens]][[s]][[as.character(subs)]] %>% 
          ggsurvfit(linewidth = 1) +
          add_risktable(theme = list(theme_risktable_default(axis.text.y.size = 11, plot.title.size = 11),
                                     theme(plot.title = element_text(face = "bold")))) +
          add_censor_mark() +
          #coord_cartesian(xlim = c(0, 25)) + 
          scale_y_continuous(labels = percent, breaks = seq(0, 1, by = .25)) +
          labs(x = "Time from randomization, months", y = "Progression Free Survival (%)", col = "Group") +
          theme_survminer() + theme(legend.position = "top")  +
          ggpp::annotate(geom = "table", x = 35, y = 1.1, 
                         label = list(medianPFS_inter_table[[sens]][[s]][[as.character((subs))]]), 
                         table.theme = ttheme_gtlight, size = 3)
          #+ expand_limits(x = 27)
        
        
        # posterior survival curves
        newdat <- expand_grid(therapy_class = therapies[-1],
                              signature = subs)
        surv_post <- posterior_survfit(fit_inter[[sens]][[s]],
                                       newdata = newdat, extrapolate = TRUE,
                                       draws = 2000, return_matrix = T,
                                       times = 0, prob = .9
        )
        time_grid <- sapply(surv_post, function(x) attr(x, "times")[1])
        surv_dat_post <- lapply(seq_along(surv_post), function(i) {
          data.frame(ndraw = seq(nrow(surv_post[[i]])),
                     time = time_grid[i],
                     surv = unname(surv_post[[i]]))
        }) %>%
          bind_rows() %>% 
          pivot_longer(cols = starts_with("surv"), names_to = "group", values_to = "surv") %>% 
          mutate(id_raw = gsub("surv.", "", group)) %>% 
          left_join(rownames_to_column(newdat, var = "id_raw"), by = "id_raw") %>% 
          mutate(
            therapy_class = factor(therapy_class, levels = names(therapies_legend[-1]), labels = therapies_legend[-1])
          ) %>% 
          group_by(time, therapy_class) %>% 
          reframe(surv = quantile(surv, c(0.5, 0.05, 0.95)),
                  q = c("median", "low", "high")) %>% 
          pivot_wider(names_from = "q", values_from = "surv", names_prefix = "surv_") 
        
        
        pkm <- pkm +
          geom_line(data = surv_dat_post, aes(time, surv_median, color = therapy_class)) +
          geom_ribbon(data = surv_dat_post, aes(x = time, y = surv_median, ymin = surv_low,
                                                ymax = surv_high, fill = therapy_class), alpha = .15) +
          scale_color_manual(labels = therapies_legend[-1], values = set_names(coltherapies[-1], therapies_legend[-1])) +
          scale_fill_manual(labels = therapies_legend[-1], values = set_names(coltherapies[-1], therapies_legend[-1])) +
          labs(col = "Group", fill = "Group", title = ifelse(subs == 1, signatures_legend[s], signatures_legend_complement[s])) +
          theme(plot.title = element_markdown())
        pkm
        
      })
    }
  })
})




## subgroup combinations ----

#### tabular results ----

tabres_subgroup <- dat_primary$main$all %>% 
  group_by(subtype, therapy_class) %>%
  summarize(
    n = n(), 
    n_1strand = sum(rand_round == 1, na.rm = T),
    n_rerand = sum(rand_round == 2, na.rm = T),
    progr = sum(prog, na.rm = T),
    progr_1strand = sum(prog[rand_round == 1], na.rm = T),
    progr_rerand = sum(prog[rand_round == 2], na.rm = T),
    PT = sum(time_obs, na.rm = T), .groups = "drop"
  ) %>%
  complete(subtype, therapy_class, fill = list(n = 0, n_1strand = 0, n_rerand = 0, PT = 0,
                                               progr = 0, progr_1strand = 0, progr_rerand = 0)) %>% 
  bind_cols(
    lapply(seq(nrow(.)), function(i){
      tidy(poisson.test(.$progr[i], T = .$PT[i], conf.level = 0.9))
    }) %>% 
      bind_rows() %>% 
      select(rate = estimate, `rate 5%` = conf.low, `rate 95%` = conf.high)
  ) %>% 
  rename(subgroup = subtype)


#### bayesian models ----

fit_subgroup <- lapply(subtype_scheme_cr$subtype, function(s) {
  dat <- filter(dat_primary[["main"]][["all"]], subtype == s)
  if (all(table(dat$therapy_class) >= 1)) {
    fit <- stan_surv(Surv(time_obs, prog) ~ therapy_class, data = dat,
                     basehaz = "weibull-aft", algorithm = "sampling",
                     prior = normal(0, .5), prior_aux = exponential(1),
                     prior_intercept = normal(3, 1), # control = list(adapt_delta = 0.99),
                     chains = 4, iter = 2000, warmup = 500, 
                     seed = 695742 + which(subtype_scheme_cr$subtype == s))
    fit
  }
}) %>% 
  set_names(nm = subtype_scheme_cr$subtype)
# warnings() # visualize warnings


### treatment effects ----

te_subgroup <- lapply(subtype_scheme_cr$subtype, function(s) {
  fit <- fit_subgroup[[s]]
  if (!is.null(fit)) {
    data.frame(
      "logSTR_ARSivsControl" = as.matrix(fit)[, "therapy_classARSi"],
      "logSTR_TaxanevsControl" = as.matrix(fit)[, "therapy_classTaxane"],
      "logSTR_ARSivsTaxane" = as.matrix(fit)[, "therapy_classARSi"] - as.matrix(fit)[, "therapy_classTaxane"],
      "shape" = as.matrix(fit)[, "weibull-shape"]
    )
  }
}) %>% 
  set_names(nm = subtype_scheme_cr$subtype) %>% 
  bind_rows(.id = "subgroup") %>% 
  pivot_longer(cols = starts_with("logSTR"), names_to = "therapy_class", 
               values_to = "logSTR") %>% 
  separate(therapy_class, into = c("therapy_class", "comparison"), sep = "vs") %>% 
  mutate(
    therapy_class = factor(gsub("logSTR_", "", therapy_class), levels = therapies),
    comparison = factor(comparison, levels = therapies[-2]),
    subgroup = factor(subgroup, levels = subtype_scheme_cr$subtype),
    logHR = -logSTR*shape
  )


# summary table
tesum_subgroup <- te_subgroup %>% 
  group_by(subgroup, therapy_class, comparison) %>% 
  reframe(
    prob_sup = mean(logSTR >= 0),
    STR = quantile(exp(logSTR), c(0.5, 0.05, 0.95)),
    HR = quantile(exp(logHR), c(0.5, 0.05, 0.95)),
    q = c("median", "low", "high")
  ) %>% 
  pivot_wider(names_from = "q", values_from = c("STR", "HR"))



## overall survival -----

### tabular results ----

tabres_os <- lapply(type_sensitivity_os, function(sens) {
  lapply(signatures_cr, function(s) { 
    dat_os[[sens]][[s]] %>%
      group_by(therapy_class) %>%
      summarize(
        n = n(), deaths = sum(dead, na.rm = T),
        PT = sum(time_os, na.rm = T), .groups = "drop"
      ) %>%
      complete(therapy_class, fill = list(n = 0, deaths = 0, PT = 0))
  }) %>% 
    bind_rows(.id = "signature")
}) %>% 
  bind_rows(.id = "type_sensitivity") %>% 
  bind_cols(
    lapply(seq(nrow(.)), function(i){
      tidy(poisson.test(.$deaths[i], T = .$PT[i], conf.level = 0.9))
    }) %>% 
      bind_rows() %>% 
      select(rate = estimate, `rate 5%` = conf.low, `rate 95%` = conf.high)
  ) %>% 
  mutate(
    signature = factor(signature, levels = signatures_cr),
    type_sensitivity = fct_inorder(type_sensitivity)
  )


### bayesian models ----

# individual seed for reproducibility
seed_os <- lapply(subtypes_sign_cr, function(s) {
  as.numeric(as.Date("2023-03-26")) + 267 +
    sum(as.integer(factor(s, levels = subtype_scheme_cr$subtype)))
})

fit_os <- lapply(type_sensitivity_os, function(sens) {
  lapply(signatures_cr, function(s) {
    if (all(table(dat_os[[sens]][[s]]$therapy_class) > 0)) {
      fit <- stan_surv(Surv(time_os, dead) ~ therapy_class, 
                       data = dat_os[[sens]][[s]],
                       basehaz = "weibull-aft", algorithm = "sampling",
                       prior = normal(0, .5), prior_aux = exponential(1),
                       prior_intercept = normal(4, 1),
                       chains = 4, iter = 2000, warmup = 500, 
                       seed = seed_os[[s]] + 
                         (seq(type_sensitivity_os) - 1)[type_sensitivity_os == sens])
      fit
    }
  })
})
# warnings() # visualize warnings


### treatment effects ----

te_os <- lapply(type_sensitivity_os, function(sens) { 
  lapply(signatures_cr, function(s) {
    fit <- fit_os[[sens]][[s]]
    if (!is.null(fit)) {
      data.frame(
        "logSTR_ARSivsControl" = as.matrix(fit)[, "therapy_classARSi"],
        "logSTR_TaxanevsControl" = as.matrix(fit)[, "therapy_classTaxane"],
        "logSTR_ARSivsTaxane" = as.matrix(fit)[, "therapy_classARSi"] - as.matrix(fit)[, "therapy_classTaxane"],
        "shape" = as.matrix(fit)[, "weibull-shape"]
      )
    }
  }) %>% 
    bind_rows(.id = "signature")
}) %>% 
  bind_rows(.id = "type_sensitivity") %>% 
  pivot_longer(cols = starts_with("logSTR"), names_to = "therapy_class", 
               values_to = "logSTR") %>% 
  separate(therapy_class, into = c("therapy_class", "comparison"), sep = "vs") %>% 
  mutate(
    therapy_class = factor(gsub("logSTR_", "", therapy_class), levels = therapies),
    comparison = factor(comparison, levels = therapies[-2]),
    signature = factor(signature, levels = signatures_cr),
    type_sensitivity = fct_inorder(type_sensitivity),
    logHR = -logSTR*shape
  )

# summary table
tesum_os <- te_os %>% 
  group_by(type_sensitivity, signature, therapy_class, comparison) %>% 
  reframe(
    prob_sup = mean(logSTR >= 0),
    STR = quantile(exp(logSTR), c(0.5, 0.05, 0.95)),
    HR = quantile(exp(logHR), c(0.5, 0.05, 0.95)),
    q = c("median", "low", "high")
  ) %>% 
  pivot_wider(names_from = "q", values_from = c("STR", "HR"))



### median OS ----
medianOS_primary <- lapply(type_sensitivity, function(sens) {
  lapply(signatures_cr, function(s) {
    
    fit <- fit_os[[sens]][[s]]
    if (!is.null(fit)) {
      log_scale_dat <- data.frame(
        "Control" = as.matrix(fit)[, c("(Intercept)")],
        "ARSi" = rowSums(as.matrix(fit)[, c("(Intercept)", "therapy_classARSi")]),
        "Taxane" = rowSums(as.matrix(fit)[, c("(Intercept)", "therapy_classTaxane")])
      )
      median_dat <- (log(2)^(1/as.matrix(fit)[, "weibull-shape"])) * exp(log_scale_dat)
      
      median_dat %>% 
        pivot_longer(cols = everything(), names_to = "therapy_class", values_to = "median")
    }
  }) %>%
    bind_rows(.id = "signature")
}) %>%
  bind_rows(.id = "type_sensitivity")

medianOS_primary_table <- lapply(type_sensitivity, function(sens) {
  lapply(signatures_cr, function(s) {
    filter(medianOS_primary, type_sensitivity == sens, signature == s) %>% 
      mutate(Group = factor(therapy_class, levels = names(therapies_legend), labels = therapies_legend)) %>% 
      group_by(Group) %>% 
      summarise(
        "Median Time \n (90% Credible Intervals)" = paste0(
          format(round(quantile(median, .5), 1), nsmall = 1), " (", 
          format(round(quantile(median, .1), 1), nsmall = 1), ", ",
          format(round(quantile(median, .95), 1), nsmall = 1), ")")
      )
  })
})


### posterior survival curves ----

km_os <- list()
postsurv_os <- lapply(type_sensitivity_os, function(sens) { 
  lapply(signatures_cr, function(s) {
    
    if (!is.null(fit_os[[sens]][[s]])) {
      # KM curves
      dat <- dat_os[[sens]][[s]] %>% 
        mutate(therapy_class = factor(therapy_class, levels = names(therapies_legend),
                                      labels = therapies_legend))
      km_os[[sens]][[s]] <<- survfit2(Surv(time_os, dead) ~ therapy_class, data = dat, conf.int = .9)
      
      pkm <- km_os[[sens]][[s]] %>% 
        ggsurvfit(linewidth = 1) +
        add_risktable(theme = list(theme_risktable_default(axis.text.y.size = 11, plot.title.size = 11),
                                   theme(plot.title = element_text(face = "bold")))) +
        add_censor_mark() +
        #coord_cartesian(xlim = c(0, 36)) + 
        scale_y_continuous(labels = percent) +
        labs(x = "Time from randomization, months", y = "Survival probability (%)", col = "Group") +
        theme_survminer() + theme(legend.position = "top") +
        ggpp::annotate(geom = "table", x = -1 + .8*(s == "all"), y = -0.19 + 0.12*(s == "all"), 
                       label = list(medianOS_primary_table[[sens]][[s]]), 
                       table.theme = ttheme_gtlight, size = 2.7 + 0.3*(s == "all")) +
        expand_limits(x = -1 + 1*(s == "all"))
      
      
      # posterior survival curves
      newdat <- data.frame(therapy_class = therapies)
      surv_post <- posterior_survfit(fit_os[[sens]][[s]],
                                     newdata = newdat, extrapolate = TRUE,
                                     draws = 2000, return_matrix = T,
                                     times = 0, prob = .9
      )
      time_grid <- sapply(surv_post, function(x) attr(x, "times")[1])
      surv_dat_post <- lapply(seq_along(surv_post), function(i) {
        data.frame(ndraw = seq(nrow(surv_post[[i]])),
                   time = time_grid[i],
                   surv = unname(surv_post[[i]]))
      }) %>%
        bind_rows() %>% 
        pivot_longer(cols = starts_with("surv"), names_to = "therapy_class", values_to = "surv") %>% 
        mutate(
          therapy_class = factor(therapy_class, levels = paste0("surv.", seq(nrow(newdat))), 
                                 labels = therapies_legend[newdat$therapy_class])
        ) %>% 
        group_by(time, therapy_class) %>% 
        reframe(surv = quantile(surv, c(0.5, 0.05, 0.95)),
                q = c("median", "low", "high")) %>% 
        pivot_wider(names_from = "q", values_from = "surv", names_prefix = "surv_") 
      
      pkm <- pkm +
        geom_line(data = surv_dat_post, aes(time, surv_median, color = therapy_class)) +
        geom_ribbon(data = surv_dat_post, aes(x = time, y = surv_median, ymin = surv_low, 
                                              ymax = surv_high, fill = therapy_class), alpha = .15) +
        scale_color_manual(labels = therapies_legend, values = set_names(coltherapies, therapies_legend)) +
        scale_fill_manual(labels = therapies_legend, values = set_names(coltherapies, therapies_legend)) +
        labs(col = "Group", fill = "Group")
      pkm
      
      
    }
  })
})



## negative controls ----

### PFS time ----

# individual seed for reproducibility
seed_negctrl <- lapply(therapies[-1], function(tx) {
  lapply(subtypes_sign_cr, function(s) {
    as.numeric(as.Date("2023-03-26")) +
      sum(as.integer(factor(s, levels = subtype_scheme_cr$subtype))) +
      nchar(tx)
  })
}) %>% setNames(nm = therapies[-1])

fit_negctrl <- lapply(therapies[-1], function(tx) {
  lapply(signatures_cr, function(s) {
    if (all(table(negativectrl[[tx]][[s]]$therapy_class) > 0)) {
      fit <- stan_surv(Surv(time_obs, prog) ~ therapy_class, 
                       data = negativectrl[[tx]][[s]],
                       basehaz = "weibull-aft", algorithm = "sampling",
                       prior = normal(0, .5), prior_aux = exponential(1),
                       prior_intercept = normal(3, 1),
                       chains = 4, iter = 2000, warmup = 500, 
                       seed = seed_negctrl[[tx]][[s]])
      fit
    }
  })
}) %>% setNames(nm = therapies[-1])
# warnings() # visualize warnings

te_negctrl <- lapply(therapies[-1], function(tx) { 
  lapply(signatures_cr, function(s) {
    fit <- fit_negctrl[[tx]][[s]]
    if (!is.null(fit)) {
      data.frame(
        "logSTR" = as.matrix(fit)[, paste0("therapy_class", tx)],
        "shape" = as.matrix(fit)[, "weibull-shape"]
      )
    }
  }) %>% 
    bind_rows(.id = "signature")
}) %>% 
  set_names(nm = therapies[-1]) %>% 
  bind_rows(.id = "therapy_class") %>% 
  mutate(
    therapy_class = factor(therapy_class, levels = therapies),
    signature = factor(signature, levels = signatures_cr),
    logHR = -logSTR*shape
  )


### OS time ----

# individual seed for reproducibility
seed_negctrl_os <- lapply(therapies[-1], function(tx) {
  lapply(subtypes_sign_cr, function(s) {
    as.numeric(as.Date("2023-03-26")) + 63 +
      sum(as.integer(factor(s, levels = subtype_scheme_cr$subtype))) +
      nchar(tx)
  })
}) %>% setNames(nm = therapies[-1])

fit_negctrl_os <- lapply(therapies[-1], function(tx) {
  lapply(signatures_cr, function(s) {
    if (all(table(negativectrl[[tx]][[s]]$therapy_class) > 0)) {
      fit <- stan_surv(Surv(time_os, dead) ~ therapy_class, 
                       data = negativectrl_os[[tx]][[s]],
                       basehaz = "weibull-aft", algorithm = "sampling",
                       prior = normal(0, .5), prior_aux = exponential(1),
                       prior_intercept = normal(3, 1),
                       chains = 4, iter = 2000, warmup = 500, 
                       seed = seed_negctrl[[tx]][[s]])
      fit
    }
  })
}) %>% setNames(nm = therapies[-1])
# warnings() # visualize warnings

te_negctrl_os <- lapply(therapies[-1], function(tx) { 
  lapply(signatures_cr, function(s) {
    fit <- fit_negctrl_os[[tx]][[s]]
    if (!is.null(fit)) {
      data.frame(
        "logSTR" = as.matrix(fit)[, paste0("therapy_class", tx)],
        "shape" = as.matrix(fit)[, "weibull-shape"]
      )
    }
  }) %>% 
    bind_rows(.id = "signature")
}) %>% 
  set_names(nm = therapies[-1]) %>% 
  bind_rows(.id = "therapy_class") %>% 
  mutate(
    therapy_class = factor(therapy_class, levels = therapies),
    signature = factor(signature, levels = signatures_cr),
    logHR = -logSTR*shape
  )


# save results ----

if (save_results) {
  
  common_names <- c("signatures_cr", "therapies", "lab", "analyses",
                    "col", "legend", "sysorgclass", "type_sensitivity", 
                    "dat_baseline", "arsi_all", "_sae", "flowchart", "id_sensitivity",
                    "tabres", "te_", "tesum_", "km_", "median", "postsurv",
                    "dat_randprob", "coef_", "te_negctrl")
  list_names <- c(unlist(sapply(common_names, function(chr) 
    grep(chr, ls(envir = .GlobalEnv), value = TRUE))), "date")
  #save(list = list_names, file = "analyses/ARSi-all/derived_data/results.RData")
  save(list = list_names, file = "analyses/ARSi-all/derived_data/NM_results.RData")
  # to include all sae
  # save(arsiall_sae, file = "analyses/ARSi-all/derived_data/arsiall_sae.RData")
  
}
