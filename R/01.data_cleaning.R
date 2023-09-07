# Import packages
library(tidyverse)


############################################################# Load data
fct_name_repair <- function(colnms) {
  tolower(
    gsub("μ", "u", 
         (gsub("\\-", "minus", 
               (gsub("\\+", "plus", colnms))
               ))
    ))
}
##########
path <- fs::path("", "Volumes", "Peres_Research","Myeloma", "Grant applications", 
                 "InNovation proposal_2022")

mif_data <- 
  readxl::read_xlsx(paste0(path, "/data/raw data",
                           "/MCC21130_P2_Myeloma_August 2023 AADL_mIF data/",
                           "Peres_P2_Myeloma_Aug2023_Results.xlsx"),
                    .name_repair = fct_name_repair) %>% 
  janitor::clean_names() %>% 
  `colnames<-`(str_remove_all(colnames(.), "opal_..._|_cells") %>% 
                 str_replace("positive", "plus") %>% 
                 str_replace("plus_plus", "plus")) %>% 
  rename(total_cells = total)

mrn_map <- 
  readxl::read_xlsx(paste0(path, "/data/raw data/AADL_MRN_map.xlsx")) %>% 
  janitor::clean_names()

clinical <- 
  read_csv(paste0(path, "/data/raw data/Abecma_Moffitt_06282023.csv")) %>% 
  janitor::clean_names()

mrn_tokeep_from_extra <- 
  read_csv(
    paste0(path, 
           "/data/processed data/clinical_mrn_to_keep_from_phi_cytokine_pt_data_07062023.csv"))

extra_clinical <-  
  read_csv(paste0(path, "/data/raw data/phi_cytokine_pt_data_07-06-2023.csv")) %>% 
  janitor::clean_names() %>% 
  filter(str_detect(mrn, mrn_tokeep_from_extra)) %>% 
  rename(date_of_car_t_infusion_number_of_days_from_apharesis_to_infusion = date_of_car_t_infusion,
         day_of_max_crs_relative_to_infusion = date_of_max_crs,
         max_icans_relative_to_infusion = max_icans,
         day_of_max_icans_relative_to_infusion = date_of_max_icans,
         crsonset_dt = date_of_crs_onset,
         crslatest_dt = latest_date_of_crs,
         crsduration = crs_duration,
         icansonset_dt = date_of_icans_onset,
         icanslatest_dt = latest_date_of_icans,
         icansduration = icans_duration)

rm(mrn_tokeep_from_extra)


############################################################# Clean data
extra_clinical <- extra_clinical %>% 
  mutate_at(c("amyloid_yes_no", "plasma_cell_leukemia_yes_no",
              "poems"), ~ case_when(
    . == "no"            ~ 0,
    . == "yes"           ~ 1
  )) %>% 
  mutate_at(c("race",
              "race_white_black_asian_pacific_islander_american_indian_alaskan_native_other_unknown"),
            ~ case_when(
    . == "caucasian"      ~ "White"
  )) %>% 
  mutate(ethnicity_yes_no_of_hispanic_latin_x_or_spanish_origin_unknown = case_when(
    ethnicity_yes_no_of_hispanic_latin_x_or_spanish_origin_unknown == 
      "nonhispanic"       ~ "no"
  ))


clinical <- clinical %>% 
  bind_rows(., extra_clinical)
write_csv(clinical, "cleaned clinical 110 patients.csv")

immune_data <- mrn_map %>% 
  full_join(., mif_data,
            by = c("image_tag", "analysis_region")) %>% 
  right_join(clinical, .,
          by = "mrn") %>% 
  mutate(mrn = as.character(mrn)) %>% 
  # Homogenize string case
  mutate(across(c(where(is.character), 
                  -c(starts_with("x")),
                  c(contains("mrd"))
                  ), ~ str_to_sentence(.) 
                )) %>% 
  mutate(across(c(contains("response")), ~ str_to_upper(.) 
  )) %>% 
  mutate(across(c(contains("response")), ~ str_replace(., "SCR", "sCR") 
  ))

immune_data1 <- immune_data %>% 
  # Fix class and Unknown
  mutate_at(c("extramedullary_disease_yes_no", 
              "ethnicity_yes_no_of_hispanic_latin_x_or_spanish_origin_unknown",
              "race_white_black_asian_pacific_islander_american_indian_alaskan_native_other_unknown"),
            ~ na_if(., "Unknown")) %>% 
  mutate_at(c("stem_cell_boost", "need_for_gcsf",
              "need_for_promacta_or_need_for_tpo_agonist"),
            ~ na_if(., "Tbd")) %>% 
  mutate(across(c(starts_with("day_of_max_"), starts_with("date"),
                  ends_with("date"), contains("_dt")
                  ),
                ~ as.Date(., format = "%m/%d/%y")
  )) %>% 
  
  # Create new var
  # infusion
  mutate(Notinfused = case_when(
    !is.na(date_of_car_t_infusion_number_of_days_from_apharesis_to_infusion) ~ "No",
    TRUE                                                                     ~ "Yes"
  )) %>% 
  mutate(time_aph_infusion = 
           interval(
             start = date_of_apheresis_or_number_of_days_from_apheresis_to_ld, 
             end = date_of_car_t_infusion_number_of_days_from_apharesis_to_infusion)/
           duration(n=1, units = "days")) %>%
  # pfs
  mutate(pfs_event = case_when(
    !is.na(date_of_pd) | !is.na(date_of_death)         ~ 1,
    TRUE                                               ~ 0
  )) %>% 
  mutate(pfs_time = coalesce(date_of_pd, date_of_death, date_of_last_contact)) %>% 
  mutate(mo_pfs_from_infusion = 
           interval(
             start = date_of_car_t_infusion_number_of_days_from_apharesis_to_infusion, 
             end = pfs_time)/
           duration(n=1, units = "months")) %>%
  # os 
  mutate(os_event = case_when(
    !is.na(date_of_death)                              ~ 1,
    TRUE                                               ~ 0
  )) %>% 
  mutate(os_time = coalesce(date_of_death, date_of_last_contact)) %>% 
  mutate(mo_os_from_infusion = 
           interval(
             start = date_of_car_t_infusion_number_of_days_from_apharesis_to_infusion, 
             end = os_time)/
           duration(n=1, units = "months")) %>% 
  mutate(agecat = case_when(
    age < 70                                           ~ "< 70 years",
    age >= 70                                          ~ "≥ 70 years"
  ),
  agecat = factor(agecat, levels = c("< 70 years", "≥ 70 years"))
  ) %>% 
  mutate(male_sex = case_when(
    sex == "Male"                                      ~ 1,
    TRUE                                               ~ 0
  )) %>% 
  # ECOG
  mutate(ecog_at_ld = case_when(
    ecog_at_ld == 0 |
      ecog_at_ld == 1                                  ~ "0-1",
    ecog_at_ld %in% c(2:4)                             ~ "2-4",
    TRUE                                               ~ NA_character_
  )) %>% 
  mutate(r_iss_at_car_t_infusion = case_when(
    r_iss_at_car_t_infusion == "Ii"                    ~ "II",
    r_iss_at_car_t_infusion == "Iii"                   ~ "III",
    TRUE                                               ~ r_iss_at_car_t_infusion
  )) %>% 
  # mutate(high_marrow_burden_50_percent= str_to_upper(high_marrow_burden_50_percent))  %>%
  # cytogenetics
  mutate(cytogenetics = case_when(
    deletion_17p_prior_to_car_t_infusion_yes_no == "Yes" | 
      t_4_14_prior_to_car_t_infusion_yes_no == "Yes" | 
      t_14_16_prior_to_car_t_infusion_yes_no == "Yes"             ~ "Yes",
    deletion_17p_prior_to_car_t_infusion_yes_no == "No" &
      t_4_14_prior_to_car_t_infusion_yes_no == "No" &
      t_14_16_prior_to_car_t_infusion_yes_no == "No"               ~ "No",
    TRUE                                                           ~ NA_character_
  )) %>% 
  # Refractory status
  mutate(immuno = case_when(
    lenalidomide_exposed_vs_refractory == "Refractory" |
      pomalidomide_exposed_vs_refractory == "Refractory"           ~ "Yes",
    TRUE                                                           ~ "No"
  )) %>% 
  mutate(proteasome = case_when(
    bortezomib_exposed_vs_refractory == "Refractory" |
      carfilzomib_exposed_vs_refractory == "Refractory"            ~ "Yes",
    TRUE                                                           ~ "No"
  )) %>% 
  mutate(dara = case_when(
    monoclonal_ab_exposed_vs_refractory == "Refractory"            ~ "Yes",
    TRUE                                                           ~ "No"
  )) %>% 
  mutate(doubleref = case_when(
    proteasome == "Yes" &
      immuno == "Yes"                                              ~ "Yes",
    TRUE                                                           ~ "No"
  )) %>% 
  mutate(tripleref = case_when(
    proteasome == "Yes" &
      immuno == "Yes" &
      monoclonal_ab_exposed_vs_refractory == "Refractory"          ~ "Yes",
    TRUE                                                           ~ "No"
  )) %>% 
  mutate(pentaref = case_when(
    bortezomib_exposed_vs_refractory == "Refractory" & 
      carfilzomib_exposed_vs_refractory == "Refractory" &
      lenalidomide_exposed_vs_refractory == "Refractory" & 
      pomalidomide_exposed_vs_refractory == "Refractory" & 
      monoclonal_ab_exposed_vs_refractory == "Refractory"          ~ "Yes",
    TRUE                                                           ~ "No"
  )) %>% 
  mutate(number_prior_lines_therapy = case_when(
    prior_lines_of_therapy <= 4                        ~ "≤ 4",
    prior_lines_of_therapy > 4                         ~ "> 4"
  )) %>% 
  
  # CRS ----
  mutate(CRS_any = ifelse(max_crs_grade_0_5 == 0, "No", "Yes"),
         CRS_any2 = ifelse(CRS_any == "No", 0, 1),
         CRS_gradecat = case_when(
           max_crs_grade_0_5 == 0                             ~ "No CRS",
           max_crs_grade_0_5 == 1                             ~ "Grade 1 or 2",
           max_crs_grade_0_5 == 2                             ~ "Grade 1 or 2",
           max_crs_grade_0_5 >= 3                             ~ "Grade ≥3"),
         CRS_gradecat = factor(CRS_gradecat, levels = c("No CRS", "Grade 1 or 2", "Grade ≥3")),
         gr3_CRS = case_when(
           max_crs_grade_0_5 %in% c(0:2)                      ~ "Grade <3",
           max_crs_grade_0_5 %in% c(3:5)                      ~ "Grade ≥3",
           TRUE                                                   ~ NA_character_
         )) %>% 
  mutate(gr2_CRS = case_when(
    max_crs_grade_0_5 %in% c(0:1)                      ~ "Grade <2",
    max_crs_grade_0_5 %in% c(2:5)                      ~ "Grade ≥2",
    TRUE                                                   ~ NA_character_
  )) %>%

  # ICANS ----
  mutate(ICANS_any = ifelse(max_icans_relative_to_infusion == 0, "No", "Yes"),
         ICANS_any2 = ifelse(ICANS_any == "No", 0, 1),
         ICANS_gradecat = case_when(
           max_icans_relative_to_infusion == 0 ~              "No ICANS",
           max_icans_relative_to_infusion == 1 ~              "Grade 1 or 2",
           max_icans_relative_to_infusion == 2 ~              "Grade 1 or 2",
           max_icans_relative_to_infusion >= 3 ~              "Grade ≥3"),
         ICANS_gradecat = factor(ICANS_gradecat, levels = c("No ICANS", "Grade 1 or 2", "Grade ≥3")),
         gr3_ICANS = case_when(
           max_icans_relative_to_infusion %in% c(0:2)         ~ "Grade <3",
           max_icans_relative_to_infusion %in% c(3:5)         ~ "Grade ≥3",
           TRUE                                                   ~ NA_character_
         )) %>% 
  mutate(gr2_ICANS = case_when(
    max_icans_relative_to_infusion %in% c(0:1)         ~ "Grade <2",
    max_icans_relative_to_infusion %in% c(2:5)         ~ "Grade ≥2",
    TRUE                                                   ~ NA_character_
  )) %>% 
  # Cytopenia ----
  mutate(Neutropenia_7d = case_when(
    anc_day_7 < 1.8                                              ~ "Yes",
    is.na(anc_day_7)                                          ~ NA_character_,
    TRUE                                                          ~ "No"
  )) %>% 
  mutate(Neutropenia_30d = case_when(
    anc_day_30 < 1.8                                             ~ "Yes",
    is.na(anc_day_30)                                          ~ NA_character_,
    TRUE                                                          ~ "No"
  )) %>% 
  mutate(Neutropenia_60d = case_when(
    anc_day_60 < 1.8                                             ~ "Yes",
    is.na(anc_day_60)                                          ~ NA_character_,
    TRUE                                                          ~ "No"
  )) %>% 
  mutate(Neutropenia_90d = case_when(
    anc_day_90 < 1.8                                             ~ "Yes",
    is.na(anc_day_90)                                          ~ NA_character_,
    TRUE                                                          ~ "No"
  )) %>% 
  mutate(gr3_Neutropenia_7d = case_when(
    anc_day_7 < 1                                              ~ "Yes",
    is.na(anc_day_7)                                          ~ NA_character_,
    TRUE                                                          ~ "No"
  )) %>% 
  mutate(gr3_Neutropenia_30d = case_when(
    anc_day_30 < 1                                             ~ "Yes",
    is.na(anc_day_30)                                          ~ NA_character_,
    TRUE                                                          ~ "No"
  )) %>% 
  mutate(gr3_Neutropenia_60d = case_when(
    anc_day_60 < 1                                             ~ "Yes",
    is.na(anc_day_60)                                          ~ NA_character_,
    TRUE                                                          ~ "No"
  )) %>% 
  mutate(gr3_Neutropenia_90d = case_when(
    anc_day_90 < 1                                             ~ "Yes",
    is.na(anc_day_90)                                          ~ NA_character_,
    TRUE                                                          ~ "No"
  )) %>% 
  mutate(Anemia_7d = case_when(
    hg_b_day_7 < 11.4                                              ~ "Yes",
    is.na(hg_b_day_7)                                          ~ NA_character_,
    TRUE                                                          ~ "No"
  )) %>%
  mutate(Anemia_30d = case_when(
    hg_b_day_30 < 11.4                                             ~ "Yes",
    is.na(hg_b_day_30)                                          ~ NA_character_,
    TRUE                                                          ~ "No"
  )) %>%
  mutate(Anemia_60d = case_when(
    hg_b_day_60 < 11.4                                             ~ "Yes",
    is.na(hg_b_day_60)                                          ~ NA_character_,
    TRUE                                                          ~ "No"
  )) %>%
  mutate(Anemia_90d = case_when(
    hg_b_day_90 < 11.4                                             ~ "Yes",
    is.na(hg_b_day_90)                                          ~ NA_character_,
    TRUE                                                          ~ "No"
  )) %>%                                                        
  mutate(gr3_Anemia_7d = case_when(
    hg_b_day_7 < 8                                              ~ "Yes",
    is.na(hg_b_day_7)                                          ~ NA_character_,
    TRUE                                                          ~ "No"
  )) %>%
  mutate(gr3_Anemia_30d = case_when(
    hg_b_day_30 < 8                                             ~ "Yes",
    is.na(hg_b_day_30)                                          ~ NA_character_,
    TRUE                                                          ~ "No"
  )) %>%
  mutate(gr3_Anemia_60d = case_when(
    hg_b_day_60 < 8                                             ~ "Yes",
    is.na(hg_b_day_60)                                          ~ NA_character_,
    TRUE                                                          ~ "No"
  )) %>%
  mutate(gr3_Anemia_90d = case_when(
    hg_b_day_90 < 8                                             ~ "Yes",
    is.na(hg_b_day_90)                                          ~ NA_character_,
    TRUE                                                          ~ "No"
  )) %>%
  mutate(Thrombocytopenia_7d = case_when(
    platelets_day_7 < 143                                              ~ "Yes",
    is.na(platelets_day_7)                                          ~ NA_character_,
    TRUE                                                          ~ "No"
  )) %>%
  mutate(Thrombocytopenia_30d = case_when(
    plateletes_day_30 < 143                                             ~ "Yes",
    is.na(plateletes_day_30)                                          ~ NA_character_,
    TRUE                                                          ~ "No"
  )) %>%
  mutate(Thrombocytopenia_60d = case_when(
    platelets_day_60 < 143                                             ~ "Yes",
    is.na(platelets_day_60)                                          ~ NA_character_,
    TRUE                                                          ~ "No"
  )) %>%
  mutate(Thrombocytopenia_90d = case_when(
    platelets_day_90 < 143                                             ~ "Yes",
    is.na(platelets_day_90)                                          ~ NA_character_,
    TRUE                                                          ~ "No"
  )) %>%
  mutate(gr3_Thrombocytopenia_7d = case_when(
    platelets_day_7 < 50                                              ~ "Yes",
    is.na(platelets_day_7)                                          ~ NA_character_,
    TRUE                                                          ~ "No"
  )) %>%
  mutate(gr3_Thrombocytopenia_30d = case_when(
    plateletes_day_30 < 50                                             ~ "Yes",
    is.na(plateletes_day_30)                                          ~ NA_character_,
    TRUE                                                          ~ "No"
  )) %>%
  mutate(gr3_Thrombocytopenia_60d = case_when(
    platelets_day_60 < 50                                             ~ "Yes",
    is.na(platelets_day_60)                                          ~ NA_character_,
    TRUE                                                          ~ "No"
  )) %>%
  mutate(gr3_Thrombocytopenia_90d = case_when(
    platelets_day_90 < 50                                             ~ "Yes",
    is.na(platelets_day_90)                                          ~ NA_character_,
    TRUE                                                          ~ "No"
  )) %>% 
  mutate(Neutropenia_gr1 = case_when(
    (anc_day_7 < 1.8 &
       anc_day_7 > 1.5) |
      (anc_day_30 < 1.8 &
         anc_day_30 > 1.5) |
      (anc_day_60 < 1.8 &
         anc_day_60 > 1.5) |
      (anc_day_90 < 1.8 &
         anc_day_90 > 1.5)                                          ~ "Yes",
    is.na(anc_day_7) |
      is.na(anc_day_30) |
      is.na(anc_day_60) |
      is.na(anc_day_90)                                          ~ NA_character_,
    TRUE                                          ~ "No"
  )) %>% 
  mutate(Neutropenia_gr2 = case_when(
    (anc_day_7 < 1.5 &
       anc_day_7 > 1) |
      (anc_day_30 < 1.5 &
         anc_day_30 > 1) |
      (anc_day_60 < 1.5 &
         anc_day_60 > 1) |
      (anc_day_90 < 1.5 &
         anc_day_90 > 1)                                          ~ "Yes",
    is.na(anc_day_7) |
      is.na(anc_day_30) |
      is.na(anc_day_60) |
      is.na(anc_day_90)                                          ~ NA_character_,
    TRUE                                          ~ "No"
  )) %>% 
  mutate(Neutropenia_gr3 = case_when(
    anc_day_7 < 1 |
      anc_day_30 < 1 |
      anc_day_60 < 1 |
      anc_day_90 < 1                                          ~ "Yes",
    is.na(anc_day_7) |
      is.na(anc_day_30) |
      is.na(anc_day_60) |
      is.na(anc_day_90)                                          ~ NA_character_,
    TRUE                                          ~ "No"
  )) %>% 
  mutate(any_neutropenia = case_when(
    Neutropenia_gr1 == "Yes" |
      Neutropenia_gr2 == "Yes" |
      Neutropenia_gr3 == "Yes"                                          ~ "Yes",
    is.na(Neutropenia_gr1) |
      is.na(Neutropenia_gr2) |
      is.na(Neutropenia_gr3)                                          ~ NA_character_,
    TRUE                                          ~ "No"
  )) %>% 
  mutate(severe_neutro_30beyond = case_when(
    anc_day_30 < 1 |
      anc_day_60 < 1 |
      anc_day_90 < 1                                          ~ "Yes",
    is.na(anc_day_30) |
      is.na(anc_day_60) |
      is.na(anc_day_90)                                          ~ NA_character_,
    TRUE                                          ~ "No"
  )) %>% 
  
  mutate(Anemia_gr1 = case_when(
    (hg_b_day_7 < 11.4 &
       hg_b_day_7 > 10) |
      (hg_b_day_30 < 11.4 &
         hg_b_day_30 > 10) |
      (hg_b_day_60 < 11.4 &
         hg_b_day_60 > 10) |
      (hg_b_day_90 < 11.4 &
         hg_b_day_90 > 10)                                        ~ "Yes",
    is.na(hg_b_day_7) |
      is.na(hg_b_day_30) |
      is.na(hg_b_day_60) |
      is.na(hg_b_day_90)                                          ~ NA_character_,
    TRUE                                          ~ "No"
  )) %>% 
  mutate(Anemia_gr2 = case_when(
    (hg_b_day_7 < 10 &
       hg_b_day_7 > 8) |
      (hg_b_day_30 < 10 &
         hg_b_day_30 > 8) |
      (hg_b_day_60 < 10 &
         hg_b_day_60 > 8) |
      (hg_b_day_90 < 10 &
         hg_b_day_90 > 8)                                        ~ "Yes",
    is.na(hg_b_day_7) |
      is.na(hg_b_day_30) |
      is.na(hg_b_day_60) |
      is.na(hg_b_day_90)                                          ~ NA_character_,
    TRUE                                          ~ "No"
  )) %>% 
  mutate(Anemia_gr3 = case_when(
    hg_b_day_7 < 8 |
      hg_b_day_30 < 8 |
      hg_b_day_60 < 8 |
      hg_b_day_90 < 8                                        ~ "Yes",
    is.na(hg_b_day_7) |
      is.na(hg_b_day_30) |
      is.na(hg_b_day_60) |
      is.na(hg_b_day_90)                                          ~ NA_character_,
    TRUE                                          ~ "No"
  )) %>% 
  mutate(severe_anemia_30beyond = case_when(
    hg_b_day_30 < 8 |
      hg_b_day_60 < 8 |
      hg_b_day_90 < 8                                        ~ "Yes",
    is.na(hg_b_day_30) |
      is.na(hg_b_day_60) |
      is.na(hg_b_day_90)                                          ~ NA_character_,
    TRUE                                          ~ "No"
  )) %>% 
  
  mutate(Thrombocytopenia_gr1 = case_when(
    (platelets_day_7 < 143 &
       platelets_day_7 > 75) |
      (plateletes_day_30 < 143 &
         plateletes_day_30 > 75) |
      (platelets_day_60 < 143 &
         platelets_day_60 > 75) |
      (platelets_day_90 < 143 &
         platelets_day_90 > 75)                                        ~ "Yes",
    is.na(platelets_day_7) |
      is.na(plateletes_day_30) |
      is.na(platelets_day_60) |
      is.na(platelets_day_90)                                          ~ NA_character_,
    TRUE                                          ~ "No"
  )) %>% 
  mutate(Thrombocytopenia_gr2 = case_when(
    (platelets_day_7 < 75 &
       platelets_day_7 > 50) |
      (plateletes_day_30 < 75 &
         plateletes_day_30 > 50) |
      (platelets_day_60 < 75 &
         platelets_day_60 > 50) |
      (platelets_day_90 < 75 &
         platelets_day_90 > 50)                                        ~ "Yes",
    is.na(platelets_day_7) |
      is.na(plateletes_day_30) |
      is.na(platelets_day_60) |
      is.na(platelets_day_90)                                          ~ NA_character_,
    TRUE                                          ~ "No"
  )) %>% 
  mutate(Thrombocytopenia_gr3 = case_when(
    platelets_day_7 < 50 |
      plateletes_day_30 < 50 |
      platelets_day_60 < 50 |
      platelets_day_90 < 50                                        ~ "Yes",
    is.na(platelets_day_7) |
      is.na(plateletes_day_30) |
      is.na(platelets_day_60) |
      is.na(platelets_day_90)                                          ~ NA_character_,
    TRUE                                          ~ "No"
  )) %>% 
  mutate(severe_throm_30beyond = case_when(
    plateletes_day_30 < 50 |
      platelets_day_60 < 50 |
      platelets_day_90 < 50                                        ~ "Yes",
    is.na(plateletes_day_30) |
      is.na(platelets_day_60) |
      is.na(platelets_day_90)                                          ~ NA_character_,
    TRUE                                          ~ "No"
  )) %>% 
  mutate(any_anemia = case_when(
    Anemia_gr1 == "Yes" |
      Anemia_gr2 == "Yes" |
      Anemia_gr3 == "Yes"                                          ~ "Yes",
    is.na(Anemia_gr1) |
      is.na(Anemia_gr2) |
      is.na(Anemia_gr3)                                          ~ NA_character_,
    TRUE                                          ~ "No"
  )) %>% 
  mutate(any_thrombo = case_when(
    Thrombocytopenia_gr1 == "Yes" |
      Thrombocytopenia_gr2 == "Yes" |
      Thrombocytopenia_gr3 == "Yes"                                          ~ "Yes",
    is.na(Thrombocytopenia_gr1) |
      is.na(Thrombocytopenia_gr2) |
      is.na(Thrombocytopenia_gr3)                                          ~ NA_character_,
    TRUE                                          ~ "No"
  )) %>% 
  mutate(gr3_cytopenia_after_30d = case_when(
    anc_day_30 < 1 |
      hg_b_day_30 < 8 |
      plateletes_day_30 < 50 |
      anc_day_60 < 1 |
      hg_b_day_60 < 8 |
      platelets_day_60 < 50 |
      anc_day_90 < 1 |
      hg_b_day_90 < 8 |
      platelets_day_90 < 50                                        ~ "Yes",
    is.na(anc_day_30) |
      is.na(hg_b_day_30) |
      is.na(plateletes_day_30) |
      is.na(anc_day_60) |
      is.na(hg_b_day_60) |
      is.na(platelets_day_60) |
      is.na(anc_day_90) |
      is.na(hg_b_day_90) |
      is.na(platelets_day_90)                                       ~ NA_character_,
    TRUE                                                              ~ "No"
  )) %>% 
  mutate(cytopenia_d7 = case_when(
    Anemia_7d == "Yes" | Thrombocytopenia_7d == "Yes" | Neutropenia_7d == "Yes" ~ "Yes",
    is.na(Anemia_7d) | is.na(Thrombocytopenia_7d) | is.na(Neutropenia_7d) ~ NA_character_,
    TRUE ~ "No"
  )) %>%
  mutate(cytopenia_d30 = case_when(
    Anemia_30d == "Yes" | Thrombocytopenia_30d == "Yes" | Neutropenia_30d == "Yes" ~ "Yes",
    is.na(Anemia_30d) | is.na(Thrombocytopenia_30d) | is.na(Neutropenia_30d) ~ NA_character_,
    TRUE ~ "No"
  )) %>%
  mutate(cytopenia_d60 = case_when(
    Anemia_60d == "Yes" | Thrombocytopenia_60d == "Yes" | Neutropenia_60d == "Yes" ~ "Yes",
    is.na(Anemia_60d) | is.na(Thrombocytopenia_60d) | is.na(Neutropenia_60d) ~ NA_character_,
    TRUE ~ "No"
  )) %>%
  mutate(cytopenia_d90 = case_when(
    Anemia_90d == "Yes" | Thrombocytopenia_90d == "Yes" | Neutropenia_90d == "Yes" ~ "Yes",
    is.na(Anemia_90d) | is.na(Thrombocytopenia_90d) | is.na(Neutropenia_90d) ~ NA_character_,
    TRUE ~ "No"
  )) %>%
  mutate(gr3_cytopenia_d7 = case_when(
    gr3_Anemia_7d == "Yes" | gr3_Thrombocytopenia_7d == "Yes" | gr3_Neutropenia_7d == "Yes" ~ "Yes",
    is.na(gr3_Anemia_7d) | is.na(gr3_Thrombocytopenia_7d) | is.na(gr3_Neutropenia_7d) ~ NA_character_,
    TRUE ~ "No"
  )) %>%
  mutate(gr3_cytopenia_d30 = case_when(
    gr3_Anemia_30d == "Yes" | gr3_Thrombocytopenia_30d == "Yes" | gr3_Neutropenia_30d == "Yes" ~ "Yes",
    is.na(gr3_Anemia_30d) | is.na(gr3_Thrombocytopenia_30d) | is.na(gr3_Neutropenia_30d) ~ NA_character_,
    TRUE ~ "No"
  )) %>%
  mutate(gr3_cytopenia_d60 = case_when(
    gr3_Anemia_60d == "Yes" | gr3_Thrombocytopenia_60d == "Yes" | gr3_Neutropenia_60d == "Yes" ~ "Yes",
    is.na(gr3_Anemia_60d) | is.na(gr3_Thrombocytopenia_60d) | is.na(gr3_Neutropenia_60d) ~ NA_character_,
    TRUE ~ "No"
  )) %>%
  mutate(gr3_cytopenia_d90 = case_when(
    gr3_Anemia_90d == "Yes" | gr3_Thrombocytopenia_90d == "Yes" | gr3_Neutropenia_90d == "Yes" ~ "Yes",
    is.na(gr3_Anemia_90d) | is.na(gr3_Thrombocytopenia_90d) | is.na(gr3_Neutropenia_90d) ~ NA_character_,
    TRUE ~ "No"
  )) %>% 
  # Cell dose ----
  mutate(cell_dose_million_cells = case_when(
    cell_dose_million_cells == 999999                  ~ NA_real_, # Not necessary for this data but keep
    TRUE                                               ~ cell_dose_million_cells
  )) %>%
  mutate(cart_cell_dose = case_when(
    cell_dose_million_cells >= 400                     ~ "≥ 400",
    cell_dose_million_cells < 400                      ~ "< 400",
    TRUE                                               ~ NA_character_
  )) %>% 
  
  # Therapy ----
  # mutate(bridging_therapy_yes_no = str_to_upper(bridging_therapy_yes_no)) %>%
  # mutate(extramedullary_disease_yes_no = case_when(
  #   extramedullary_disease_yes_no == "Unknown" ~ NA_character_,
  #   TRUE ~ extramedullary_disease_yes_no
  # )) %>%
  # mutate(extramedullary_disease_yes_no = str_to_upper(extramedullary_disease_yes_no)) %>%
  # mutate(prior_treatment_with_any_gene_therapy_based_therapeutic_or_investigational_cellular_therapy_or_bcma_targeted_therapy_yes_no = str_to_upper(prior_treatment_with_any_gene_therapy_based_therapeutic_or_investigational_cellular_therapy_or_bcma_targeted_therapy_yes_no)) %>%
  # mutate(anakinra_yes_no = str_to_upper(anakinra_yes_no)) %>%
  # mutate(non_secretory_mm_yes_no = str_to_upper(non_secretory_mm_yes_no)) %>%
  mutate(plasma_cell_leukemia_yes_no = case_when(
    plasma_cell_leukemia_yes_no == 0 ~ "No",
    plasma_cell_leukemia_yes_no == 1 ~ "Yes",
    is.na(plasma_cell_leukemia_yes_no) ~ NA_character_
  )) %>%
  mutate(amyloid_yes_no = case_when(
    amyloid_yes_no == 0 ~ "No",
    amyloid_yes_no == 1 ~ "Yes",
    is.na(amyloid_yes_no) ~ NA_character_
  )) %>%
  mutate(poems = case_when(
    poems == 0 ~ "No",
    poems == 1 ~ "Yes",
    is.na(poems) ~ NA_character_
  )) %>%
  # mutate(prior_auto_sct_yes_no = str_to_upper(prior_auto_sct_yes_no)) %>%
  # mutate(asct_within_12_weeks_prior_to_apheresis_yes_no = str_to_upper(asct_within_12_weeks_prior_to_apheresis_yes_no)) %>%
  # mutate(prior_allo_sct_yes_no = str_to_upper(prior_allo_sct_yes_no)) %>%
  mutate(inpatientbed = case_when(
    icu_admission_yes_no == "No" | is.na(icu_admission_yes_no) ~ length_of_hospital_stay_total_days_including_readmission,
    icu_admission_yes_no == "Yes" ~ (length_of_hospital_stay_total_days_including_readmission - length_of_icu_stay)
  )) %>%
  mutate(inpatientbed_costs = inpatientbed * 3001.80) %>% 
  mutate(icubed_costs = length_of_icu_stay * 9531.52) %>%
  mutate(length_icu_stay = case_when(
    icu_admission_yes_no == "No" | is.na(icu_admission_yes_no) ~ NA_real_,
    icu_admission_yes_no == "Yes" ~ length_of_icu_stay
  )) %>%
  mutate(icubed_costs_amongicu = case_when(
    icu_admission_yes_no == "Yes" ~ icubed_costs,
    TRUE ~ NA_real_
  )) %>%
  mutate(total_costs = case_when(
    icu_admission_yes_no == "Yes" ~ (inpatientbed_costs + icubed_costs),
    TRUE ~ inpatientbed_costs)) %>%
  # mutate(icu_admission_yes_no = str_to_upper(icu_admission_yes_no)) %>% 
  # mutate(tocilizumab_yes_no = str_to_upper(tocilizumab_yes_no)) %>% 
  # mutate(steroid_use_yes_no = str_to_upper(steroid_use_yes_no)) %>%
  mutate(bridging_therapy_yes_no = case_when( # Not needed for this data but keep
    bridging_therapy_yes_no == "Xrt" ~ "Yes",
    TRUE ~ bridging_therapy_yes_no
  )) %>% 
  mutate(readmit = case_when(
    readmission_duration > 0 ~ "Yes",
    TRUE ~ "No"
  )) %>%
  mutate(typeinfection = case_when(
    is.na(infection_yes_no)                            ~ NA_character_,
    infection_yes_no == "No" ~ "No infection",
    infection_yes_no == "Yes" & 
      is.na(infection_type) ~ "Unknown type of infection",
    TRUE ~ infection_type
  )) %>%
  mutate(crsonset = crsonset_dt - date_of_car_t_infusion_number_of_days_from_apharesis_to_infusion) %>%
  mutate(icansonset = icansonset_dt - date_of_car_t_infusion_number_of_days_from_apharesis_to_infusion) %>%
  mutate(infectiononset = infection_dt - date_of_car_t_infusion_number_of_days_from_apharesis_to_infusion)

# Response
immune_data2 <- immune_data1 %>% 
  mutate(day30response = case_when(
    (os_event == 1 | pfs_event ==1) & mo_pfs_from_infusion < 1 ~ "Died or progressed before day 30",
    #pfs_event == 1 & mo_pfs_from_infusion < 1 ~ "Died or progressed before day 30",
    str_detect(day_30_response_s_cr_cr_vgpr_pr, "CR") ~ "sCR or CR",
    str_detect(day_30_response_s_cr_cr_vgpr_pr, "VGPR") ~ "VGPR",
    str_detect(day_30_response_s_cr_cr_vgpr_pr, "PR") ~ "PR",
    str_detect(day_30_response_s_cr_cr_vgpr_pr, "MR") | str_detect(day_30_response_s_cr_cr_vgpr_pr, "SD") ~ "SD",
    str_detect(day_30_response_s_cr_cr_vgpr_pr, "PD") ~ "PD",
    is.na(day_30_response_s_cr_cr_vgpr_pr) & os_event == 0 & mo_pfs_from_infusion < 1 ~ "Day 30 not reached",
    TRUE ~ "Not provided"
  )) %>%
  mutate(mon3response = case_when(
    (os_event == 1 | pfs_event == 1) & mo_pfs_from_infusion < 3 ~ "Died or progressed before 3 months",
    str_detect(x3_month_response_s_cr_cr_vgpr_pr, "CR") ~ "sCR or CR",
    str_detect(x3_month_response_s_cr_cr_vgpr_pr, "VGPR") ~ "VGPR",
    str_detect(x3_month_response_s_cr_cr_vgpr_pr, "PR") ~ "PR",
    str_detect(x3_month_response_s_cr_cr_vgpr_pr, "MR") | str_detect(x3_month_response_s_cr_cr_vgpr_pr, "SD") ~ "SD",
    str_detect(x3_month_response_s_cr_cr_vgpr_pr, "PD") ~ "PD",
    is.na(x3_month_response_s_cr_cr_vgpr_pr) & os_event == 0 & mo_pfs_from_infusion < 3 ~ "3 months not reached",
    TRUE ~ "Not provided"
  )) %>%
  mutate(mon6response = case_when(
    (os_event == 1 | pfs_event == 1) & mo_pfs_from_infusion < 6 ~ "Died or progressed before 6 months",
    str_detect(x6_month_response_s_cr_cr_vgpr_pr, "CR") ~ "sCR or CR",
    str_detect(x6_month_response_s_cr_cr_vgpr_pr, "VGPR") ~ "VGPR",
    str_detect(x6_month_response_s_cr_cr_vgpr_pr, "PR") ~ "PR",
    str_detect(x6_month_response_s_cr_cr_vgpr_pr, "MR") | str_detect(x6_month_response_s_cr_cr_vgpr_pr, "SD") ~ "SD",
    str_detect(x6_month_response_s_cr_cr_vgpr_pr, "PD") ~ "PD",
    is.na(x6_month_response_s_cr_cr_vgpr_pr) & os_event == 0 & mo_pfs_from_infusion < 6 ~ "6 months not reached",
    TRUE ~ "Not provided"
  )) %>%
  mutate(mon9response = case_when(
    (os_event == 1 | pfs_event) & mo_pfs_from_infusion < 9 ~ "Died or progressed before 9 months",
    str_detect(x9_month_response_s_cr_cr_vgpr_pr, "CR") ~ "sCR or CR",
    str_detect(x9_month_response_s_cr_cr_vgpr_pr, "VGPR") ~ "VGPR",
    str_detect(x9_month_response_s_cr_cr_vgpr_pr, "PR") ~ "PR",
    str_detect(x9_month_response_s_cr_cr_vgpr_pr, "MR") | str_detect(x9_month_response_s_cr_cr_vgpr_pr, "SD") ~ "SD",
    str_detect(x9_month_response_s_cr_cr_vgpr_pr, "PD") ~ "PD",
    is.na(x9_month_response_s_cr_cr_vgpr_pr) & os_event == 0 & mo_pfs_from_infusion < 9 ~ "9 months not reached",
    TRUE ~ "Not provided"
  )) %>%
  mutate(mon12response = case_when(
    (os_event == 1 | pfs_event == 1) & mo_pfs_from_infusion < 12 ~ "Died or progressed before 12 months",
    str_detect(x12_month_response_s_cr_cr_vgpr_pr, "CR") ~ "sCR or CR",
    str_detect(x12_month_response_s_cr_cr_vgpr_pr, "VGPR") ~ "VGPR",
    str_detect(x12_month_response_s_cr_cr_vgpr_pr, "PR") ~ "PR",
    str_detect(x12_month_response_s_cr_cr_vgpr_pr, "MR") | str_detect(x12_month_response_s_cr_cr_vgpr_pr, "SD") ~ "SD",
    str_detect(x12_month_response_s_cr_cr_vgpr_pr, "PD") ~ "PD",
    is.na(x12_month_response_s_cr_cr_vgpr_pr) & os_event == 0 & mo_pfs_from_infusion < 12 ~ "12 months not reached",
    TRUE ~ "Not provided"
  )) %>%
  mutate(mon15response = case_when(
    (os_event == 1 | pfs_event == 1) & mo_pfs_from_infusion < 15 ~ "Died or progressed before 15 months",
    str_detect(x15_month_response_s_cr_cr_vgpr_pr, "CR") ~ "sCR or CR",
    str_detect(x15_month_response_s_cr_cr_vgpr_pr, "VGPR") ~ "VGPR",
    str_detect(x15_month_response_s_cr_cr_vgpr_pr, "PR") ~ "PR",
    str_detect(x15_month_response_s_cr_cr_vgpr_pr, "MR") | str_detect(x15_month_response_s_cr_cr_vgpr_pr, "SD") ~ "SD",
    str_detect(x15_month_response_s_cr_cr_vgpr_pr, "PD") ~ "PD",
    is.na(x15_month_response_s_cr_cr_vgpr_pr) & os_event == 0 & mo_pfs_from_infusion < 15 ~ "15 months not reached",
    TRUE ~ "Not provided"
  )) %>%
  mutate(mon18response = case_when(
    (os_event == 1 | pfs_event == 1) & mo_pfs_from_infusion < 18 ~ "Died or progressed before 18 months",
    str_detect(x18_month_response_s_cr_cr_vgpr_pr, "CR") ~ "sCR or CR",
    str_detect(x18_month_response_s_cr_cr_vgpr_pr, "VGPR") ~ "VGPR",
    str_detect(x18_month_response_s_cr_cr_vgpr_pr, "PR") ~ "PR",
    str_detect(x18_month_response_s_cr_cr_vgpr_pr, "MR") | str_detect(x15_month_response_s_cr_cr_vgpr_pr, "SD") ~ "SD",
    str_detect(x18_month_response_s_cr_cr_vgpr_pr, "PD") ~ "PD",
    is.na(x18_month_response_s_cr_cr_vgpr_pr) & os_event == 0 & mo_pfs_from_infusion < 18 ~ "18 months not reached",
    TRUE ~ "Not provided"
  )) %>%
  # if a PD then always a PD
  mutate(mon3response = case_when(
    (day30response == "PD" | day30response == "Died or progressed before day 30") ~ "PD",
    TRUE ~ mon3response
  )) %>%
  mutate(mon6response = case_when(
    (mon3response == "PD" | mon3response == "Died or progressed before 3 months") ~ "PD",
    TRUE ~ mon6response
  )) %>%
  mutate(mon9response = case_when(
    (mon6response == "PD" | mon6response == "Died or progressed before 6 months") ~ "PD",
    TRUE ~ mon9response
  )) %>%
  mutate(mon12response = case_when(
    (mon9response == "PD" | mon9response == "Died or progressed before 9 months") ~ "PD",
    TRUE ~ mon12response
  )) %>%
  mutate(mon15response = case_when(
    (mon12response == "PD" | mon12response == "Died or progressed before 12 months") ~ "PD",
    TRUE ~ mon15response
  )) %>%
  mutate(mon18response = case_when(
    (mon15response == "PD" | mon15response == "Died or progressed before 15 months") ~ "PD",
    TRUE ~ mon18response
  )) %>%
  # include patients who died or progressed before day 30 as a PD response
  # do not include patients who did not reach Day 30, who reached day 30 but did not provide a response, and those that are not evaluable in the denominator for response calculations
  mutate(day30response_v2 = case_when(
    day30response == "Died or progressed before day 30"  ~ "PD",
    day30response == "Day 30 not reached" |  day30response == "Not provided" | day30response == "Not Evaluable" ~ NA_character_, 
    TRUE ~ day30response
  )) %>%
  mutate(ORR_30days = case_when(
    str_detect(day30response_v2, "CR")   ~ "ORR",
    str_detect(day30response_v2, "PR")   ~ "ORR",
    is.na(day30response_v2) ~ NA_character_,
    TRUE                                                      ~ "SD/PD"
  )) %>%
  mutate(CRorbetter_30days = case_when(
    str_detect(day30response_v2, "CR")   ~ "CR or better",
    is.na(day30response_v2)               ~ NA_character_,
    TRUE                                                        ~ "< CR"
  )) %>%
  # MRD this replace what is commented out below
  mutate(across(c(contains("mrd")), ~ case_when(
    . == "Mrd_neg" ~ "Negative",
    . == "Mrd_pos" ~ "Positive",
    . == "Mrd_unk" ~ NA_character_,
    TRUE ~ .
  )
  )) %>% 
  # mutate(day_30_mrd_clonoseq_positive_vs_negative = case_when(
  #   day_30_mrd_clonoseq_positive_vs_negative == "mrd_neg" ~ "Negative",
  #   day_30_mrd_clonoseq_positive_vs_negative == "mrd_pos" ~ "Positive",
  #   day_30_mrd_clonoseq_positive_vs_negative == "mrd_unk" ~ NA_character_,
  #   TRUE ~ day_30_mrd_clonoseq_positive_vs_negative
  # )) %>%
  # mutate(x3_month_mrd_clonoseq_positive_vs_negative = case_when(
  #   x3_month_mrd_clonoseq_positive_vs_negative == "mrd_neg" ~ "Negative",
  #   x3_month_mrd_clonoseq_positive_vs_negative == "mrd_pos" ~ "Positive",
  #   x3_month_mrd_clonoseq_positive_vs_negative == "mrd_unk" ~ NA_character_,
  #   TRUE ~ x3_month_mrd_clonoseq_positive_vs_negative
  # )) %>%
  # mutate(x6_month_mrd_clonoseq_positive_vs_negative = case_when(
  #   x6_month_mrd_clonoseq_positive_vs_negative == "mrd_neg" ~ "Negative",
  #   x6_month_mrd_clonoseq_positive_vs_negative == "mrd_pos" ~ "Positive",
  #   x6_month_mrd_clonoseq_positive_vs_negative == "mrd_unk" ~ NA_character_,
  #   TRUE ~ x6_month_mrd_clonoseq_positive_vs_negative
  # )) %>%
  # mutate(x9_month_mrd_clonoseq_positive_vs_negative = case_when(
  #   x9_month_mrd_clonoseq_positive_vs_negative == "mrd_neg" ~ "Negative",
  #   x9_month_mrd_clonoseq_positive_vs_negative == "mrd_pos" ~ "Positive",
  #   x9_month_mrd_clonoseq_positive_vs_negative == "mrd_unk" ~ NA_character_,
  #   TRUE ~ x9_month_mrd_clonoseq_positive_vs_negative
  # )) %>%
  # mutate(x12_month_mrd_clonoseq_positive_vs_negative = case_when(
  #   x12_month_mrd_clonoseq_positive_vs_negative == "mrd_neg" ~ "Negative",
  #   x12_month_mrd_clonoseq_positive_vs_negative == "mrd_pos" ~ "Positive",
  #   x12_month_mrd_clonoseq_positive_vs_negative == "mrd_unk" ~ NA_character_,
  #   TRUE ~ x12_month_mrd_clonoseq_positive_vs_negative
  # )) %>%
  
  
  mutate(day30response_MRD = case_when(
    day30response_v2 == "sCR or CR" & day_30_mrd_clonoseq_positive_vs_negative == "Negative" ~ "sCR or CR, MRD-",
    day30response_v2 == "sCR or CR" & day_30_mrd_clonoseq_positive_vs_negative == "Positive" ~ "sCR or CR, MRD+",
    day30response_v2 == "sCR or CR" & is.na(day_30_mrd_clonoseq_positive_vs_negative) ~ "sCR or CR, MRD unknown",
    TRUE ~ day30response_v2
  )) %>%
  mutate(mon3response_v2 = case_when(
    mon3response == "Died or progressed before 3 months" ~ "PD" ,
    mon3response == "3 months not reached" |  mon3response == "Not provided" | mon3response == "Not Evaluable" ~ NA_character_, 
    TRUE ~ mon3response
  )) %>%
  mutate(ORR_3mon = case_when(
    str_detect(mon3response_v2, "CR")   ~ "ORR",
    str_detect(mon3response_v2, "PR")   ~ "ORR",
    is.na(mon3response_v2) ~ NA_character_,
    TRUE                                                      ~ "SD/PD"
  )) %>%
  mutate(CRorbetter_3mon = case_when(
    str_detect(mon3response_v2, "CR")   ~ "CR or better",
    is.na(mon3response_v2)               ~ NA_character_,
    TRUE                                                        ~ "< CR"
  )) %>%
  mutate(mon3response_MRD = case_when(
    mon3response_v2 == "sCR or CR" & x3_month_mrd_clonoseq_positive_vs_negative == "Negative" ~ "sCR or CR, MRD-",
    mon3response_v2 == "sCR or CR" & x3_month_mrd_clonoseq_positive_vs_negative == "Positive" ~ "sCR or CR, MRD+",
    mon3response_v2 == "sCR or CR" & is.na(x3_month_mrd_clonoseq_positive_vs_negative) ~ "sCR or CR, MRD unknown",
    TRUE ~ mon3response_v2
  )) %>%
  mutate(mon6response_v2 = case_when(
    mon6response == "Died or progressed before 6 months"  ~ "PD",
    mon6response == "6 months not reached" |  mon6response == "Not provided" | mon6response == "Not Evaluable" ~ NA_character_, 
    TRUE ~ mon6response
  )) %>%
  mutate(ORR_6mon = case_when(
    str_detect(mon6response_v2, "CR")   ~ "ORR",
    str_detect(mon6response_v2, "PR")   ~ "ORR",
    is.na(mon6response_v2) ~ NA_character_,
    TRUE                                                      ~ "SD/PD"
  )) %>%
  mutate(CRorbetter_6mon = case_when(
    str_detect(mon6response_v2, "CR")   ~ "CR or better",
    is.na(mon6response_v2)               ~ NA_character_,
    TRUE                                                        ~ "< CR"
  )) %>%
  mutate(mon6response_MRD = case_when(
    mon6response_v2 == "sCR or CR" & x6_month_mrd_clonoseq_positive_vs_negative == "Negative" ~ "sCR or CR, MRD-",
    mon6response_v2 == "sCR or CR" & x6_month_mrd_clonoseq_positive_vs_negative == "Positive" ~ "sCR or CR, MRD+",
    mon6response_v2 == "sCR or CR" & is.na(x6_month_mrd_clonoseq_positive_vs_negative) ~ "sCR or CR, MRD unknown",
    TRUE ~ mon6response_v2
  )) %>%
  mutate(mon9response_v2 = case_when(
    mon9response == "Died or progressed before 9 months"  ~ "PD",
    mon9response == "9 months not reached" |  mon9response == "Not provided" | mon9response == "Not Evaluable" ~ NA_character_, 
    TRUE ~ mon9response
  )) %>%
  mutate(ORR_9mon = case_when(
    str_detect(mon9response_v2, "CR")   ~ "ORR",
    str_detect(mon9response_v2, "PR")   ~ "ORR",
    is.na(mon9response_v2) ~ NA_character_,
    TRUE                                                      ~ "SD/PD"
  )) %>%
  mutate(CRorbetter_9mon = case_when(
    str_detect(mon9response_v2, "CR")   ~ "CR or better",
    is.na(mon9response_v2)               ~ NA_character_,
    TRUE                                                        ~ "< CR"
  )) %>%
  mutate(mon9response_MRD = case_when(
    mon9response_v2 == "sCR or CR" & x9_month_mrd_clonoseq_positive_vs_negative == "Negative" ~ "sCR or CR, MRD-",
    mon9response_v2 == "sCR or CR" & x9_month_mrd_clonoseq_positive_vs_negative == "Positive" ~ "sCR or CR, MRD+",
    mon9response_v2 == "sCR or CR" & is.na(x9_month_mrd_clonoseq_positive_vs_negative) ~ "sCR or CR, MRD unknown",
    TRUE ~ mon9response_v2
  )) %>%
  mutate(mon12response_v2 = case_when(
    mon12response == "Died or progressed before 12 months"  ~ "PD",
    mon12response == "12 months not reached" | mon12response == "Not provided" | mon12response == "Not Evaluable" ~ NA_character_, 
    TRUE ~ mon12response
  )) %>%
  mutate(ORR_12mon = case_when(
    str_detect(mon12response_v2, "CR")   ~ "ORR",
    str_detect(mon12response_v2, "PR")   ~ "ORR",
    is.na(mon12response_v2) ~ NA_character_,
    TRUE                                                      ~ "SD/PD"
  )) %>%
  mutate(CRorbetter_12mon = case_when(
    str_detect(mon12response_v2, "CR")   ~ "CR or better",
    is.na(mon12response_v2)               ~ NA_character_,
    TRUE                                                        ~ "< CR"
  )) %>%
  mutate(mon12response_MRD = case_when(
    mon12response_v2 == "sCR or CR" & x12_month_mrd_clonoseq_positive_vs_negative == "Negative" ~ "sCR or CR, MRD-",
    mon12response_v2 == "sCR or CR" & x12_month_mrd_clonoseq_positive_vs_negative == "Positive" ~ "sCR or CR, MRD+",
    mon12response_v2 == "sCR or CR" & is.na(x12_month_mrd_clonoseq_positive_vs_negative) ~ "sCR or CR, MRD unknown",
    TRUE ~ mon12response_v2
  )) %>%
  mutate(mon15response_v2 = case_when(
    mon15response == "Died or progressed before 15 months"  ~ "PD",
    mon15response == "15 months not reached" | mon15response == "Not provided" | mon15response == "Not Evaluable" ~ NA_character_, 
    TRUE ~ mon15response
  )) %>%
  mutate(ORR_15mon = case_when(
    str_detect(mon15response_v2, "CR")   ~ "ORR",
    str_detect(mon15response_v2, "PR")   ~ "ORR",
    is.na(mon15response_v2) ~ NA_character_,
    TRUE                                                      ~ "SD/PD"
  )) %>%
  mutate(CRorbetter_15mon = case_when(
    str_detect(mon15response_v2, "CR")   ~ "CR or better",
    is.na(mon15response_v2)               ~ NA_character_,
    TRUE                                                        ~ "< CR"
  )) %>%
  mutate(mon15response_MRD = case_when(
    mon15response_v2 == "sCR or CR" & x15_month_mrd_clonoseq_positive_vs_negative == "Negative" ~ "sCR or CR, MRD-",
    mon15response_v2 == "sCR or CR" & x15_month_mrd_clonoseq_positive_vs_negative == "Positive" ~ "sCR or CR, MRD+",
    mon15response_v2 == "sCR or CR" & is.na(x15_month_mrd_clonoseq_positive_vs_negative) ~ "sCR or CR, MRD unknown",
    TRUE ~ mon15response_v2
  )) %>%
  mutate(mon18response_v2 = case_when(
    mon18response == "Died or progressed before 18 months"  ~ "PD",
    mon18response == "18 months not reached" | mon18response == "Not provided" | mon18response == "Not Evaluable" ~ NA_character_, 
    TRUE ~ mon18response
  )) %>%
  mutate(ORR_18mon = case_when(
    str_detect(mon18response_v2, "CR")   ~ "ORR",
    str_detect(mon18response_v2, "PR")   ~ "ORR",
    is.na(mon18response_v2) ~ NA_character_,
    TRUE                                                      ~ "SD/PD"
  )) %>%
  mutate(CRorbetter_18mon = case_when(
    str_detect(mon18response_v2, "CR")   ~ "CR or better",
    is.na(mon18response_v2)               ~ NA_character_,
    TRUE                                                        ~ "< CR"
  )) %>%
  mutate(mon18response_MRD = case_when(
    mon18response_v2 == "sCR or CR" & x18_month_mrd_clonoseq_positive_vs_negative == "Negative" ~ "sCR or CR, MRD-",
    mon18response_v2 == "sCR or CR" & x18_month_mrd_clonoseq_positive_vs_negative == "Positive" ~ "sCR or CR, MRD+",
    mon18response_v2 == "sCR or CR" & is.na(x18_month_mrd_clonoseq_positive_vs_negative) ~ "sCR or CR, MRD unknown",
    TRUE ~ mon18response_v2
  )) %>%
  mutate(best_ORR = case_when(
    ORR_30days == "ORR" | ORR_3mon == "ORR" | ORR_6mon == "ORR" | ORR_9mon == "ORR" | ORR_12mon == "ORR" | ORR_15mon == "ORR" | ORR_18mon == "ORR" ~ "ORR",
    is.na(ORR_30days) & is.na(ORR_3mon) & is.na(ORR_6mon) & is.na(ORR_9mon) & is.na(ORR_12mon) & is.na(ORR_15mon) & is.na(ORR_18mon) ~ NA_character_, 
    TRUE                                                               ~ "SD/PD"
  )) %>% 
  mutate(best_CRorbetter = case_when(
    CRorbetter_30days == "CR or better" | CRorbetter_3mon == "CR or better" | CRorbetter_6mon == "CR or better" | CRorbetter_9mon == "CR or better" | CRorbetter_12mon == "CR or better" | CRorbetter_18mon == "CR or better" | CRorbetter_18mon == "CR or better" ~ "CR or better",
    is.na(CRorbetter_30days) & is.na(CRorbetter_3mon) & is.na(CRorbetter_6mon) & is.na(CRorbetter_9mon) & is.na(CRorbetter_12mon) & is.na(CRorbetter_15mon) & is.na(CRorbetter_18mon) ~ NA_character_, 
    TRUE                                                               ~ "< CR"
  ))  %>%
  mutate(best_response = case_when (
    day30response_v2 == "sCR or CR" | mon3response_v2 == "sCR or CR" | mon6response_v2 == "sCR or CR" | mon9response_v2 == "sCR or CR" | mon12response_v2 == "sCR or CR" | mon15response_v2 == "sCR or CR" | mon15response_v2 == "sCR or CR" ~ "sCR or CR",
    day30response_v2 == "VGPR" | mon3response_v2 == "VGPR" | mon6response_v2 == "VGPR" | mon9response_v2 == "VGPR" | mon12response_v2 == "VGPR" | mon15response_v2 == "VGPR" | mon18response_v2 == "VGPR" ~ "VGPR",
    day30response_v2 == "PR" | mon3response_v2 == "PR" | mon6response_v2 == "PR" | mon9response_v2 == "PR" | mon12response_v2 == "PR" | mon15response_v2 == "PR" | mon18response_v2 == "PR"~ "PR",
    day30response_v2 == "SD" | mon3response_v2 == "SD" | mon6response_v2 == "SD" | mon9response_v2 == "SD" | mon12response_v2 == "SD" | mon15response_v2 == "SD" | mon18response_v2 == "SD"~ "SD",
    day30response_v2 == "PD" | mon3response_v2 == "PD" | mon6response_v2 == "PD" | mon9response_v2 == "PD" | mon12response_v2 == "PD" | mon15response_v2 == "PD" | mon18response_v2 == "PD" ~ "PD",
    TRUE ~ NA_character_
  )) %>%
  mutate(best_response_MRD = case_when (
    day30response_MRD == "sCR or CR, MRD-" | mon3response_MRD == "sCR or CR, MRD-" | mon6response_MRD == "sCR or CR, MRD-" | mon9response_MRD == "sCR or CR, MRD-" | mon12response_MRD == "sCR or CR, MRD-" | mon15response_MRD == "sCR or CR, MRD-" | mon18response_MRD == "sCR or CR, MRD-" ~ "sCR or CR, MRD-",
    day30response_MRD == "sCR or CR, MRD+" | mon3response_MRD == "sCR or CR, MRD+" | mon6response_MRD == "sCR or CR, MRD+" | mon9response_MRD == "sCR or CR, MRD+" | mon12response_MRD == "sCR or CR, MRD+" | mon15response_MRD == "sCR or CR, MRD+" | mon18response_MRD == "sCR or CR, MRD+" ~ "sCR or CR, MRD+",
    day30response_MRD == "sCR or CR, MRD unknown" | mon3response_MRD == "sCR or CR, MRD unknown" | mon6response_MRD == "sCR or CR, MRD unknown" | mon9response_MRD == "sCR or CR, MRD unknown" | mon12response_MRD == "sCR or CR, MRD unknown" | mon15response_MRD == "sCR or CR, MRD unknown" | mon18response_MRD == "sCR or CR, MRD unknown" ~ "sCR or CR, MRD unknown",
    TRUE ~ best_response
  ))

immune_data3 <- immune_data2 %>% 
  # Labs
  mutate(baseline_ferritin = case_when(
    baseline_ferritin >= 400                              ~ "≥ ULN at LD (≥ 400)",
    baseline_ferritin < 400                               ~ "Normal (< 400)"
  )) %>% 
  mutate(baseline_crp = case_when(
    baseline_crp >= 0.5                                   ~ "≥ ULN at LD (≥ 0.5)",
    baseline_crp < 0.5                                    ~ "Normal (< 0.5)"
  )) %>%
  mutate(baseline_B2M = case_when(
    baseline_b2m >= 5.5                                   ~ "≥ 5.5",
    baseline_b2m < 5.5                                    ~ "< 5.5"
  )) %>%
  # mutate(renal_insufficiency_cr_cl_45_m_l_min_yes_no = str_to_upper(renal_insufficiency_cr_cl_45_m_l_min_yes_no)) %>%
  mutate(excl_orgdys = case_when(
    lvef_45_percent_yes_no == "Yes" | lvef_45_percent_yes_no == "Yes" |
      ast_alt_2_5_uln_yes_no == "Yes" | ast_alt_2_5_uln_yes_no == "Yes" |
      total_bilirubin_1_5uln_yes_no == "Yes" | total_bilirubin_1_5uln_yes_no == "Yes" |
      renal_insufficiency_cr_cl_45_m_l_min_yes_no == "Yes" |
      str_detect(h_o_class_iii_or_iv_chf_or_severe_nonischemic_cardiomyopathy_unstable_or_poorly_controlled_angina_myocardial_infarction_or_ventricular_arrhythmia_within_6_months_yes_no, "Yes") | 
      str_detect(h_o_class_iii_or_iv_chf_or_severe_nonischemic_cardiomyopathy_unstable_or_poorly_controlled_angina_myocardial_infarction_or_ventricular_arrhythmia_within_6_months_yes_no, "Yes")~ "Yes",
    lvef_45_percent_yes_no == "Unknown" |
      ast_alt_2_5_uln_yes_no == "Unknown" |
      total_bilirubin_1_5uln_yes_no == "Unknown" |
      renal_insufficiency_cr_cl_45_m_l_min_yes_no == "Unknown" |
      str_detect(h_o_class_iii_or_iv_chf_or_severe_nonischemic_cardiomyopathy_unstable_or_poorly_controlled_angina_myocardial_infarction_or_ventricular_arrhythmia_within_6_months_yes_no, "Unknown") ~ NA_character_,
    TRUE ~ "No"
  )) %>%
  mutate(subtype = case_when(
    non_secretory_mm_yes_no == "Yes" ~ "Non-secretory",
    myeloma_sub_type == "kappa" | myeloma_sub_type =="lambda" ~ "Light chain only",
    non_secretory_mm_yes_no == "No" & is.na(myeloma_sub_type) ~ NA_character_,
    TRUE ~ "Intact immunoglobulin (Heavy + Light)" 
  )) %>%
  # mutate(eth_clean = case_when(
  #   ethnicity_yes_no_of_hispanic_latin_x_or_spanish_origin_unknown == "Unknown" ~ NA_character_,
  #   TRUE ~ ethnicity_yes_no_of_hispanic_latin_x_or_spanish_origin_unknown
  # )) %>%
  mutate(race_eth = case_when(
    ethnicity_yes_no_of_hispanic_latin_x_or_spanish_origin_unknown == "Yes" ~ "Hispanic (All races)",
    TRUE ~ race_white_black_asian_pacific_islander_american_indian_alaskan_native_other_unknown
  )) #%>% 
  # mutate(race_white_black_asian_pacific_islander_american_indian_alaskan_native_other_unknown = case_when(
  #   race_white_black_asian_pacific_islander_american_indian_alaskan_native_other_unknown == "Unknown" ~ NA_character_,
  #   TRUE ~ race_white_black_asian_pacific_islander_american_indian_alaskan_native_other_unknown
  # ))

immune_data <- immune_data3 %>% 
  select(-c(image_tag : area_analyzed_mm2), 
         image_tag : area_analyzed_mm2)


write_rds(immune_data, paste0(here::here(), "/immune_data.rds"))
write_csv(immune_data, paste0(here::here(), "/immune_data.csv"))


check_data <- function(data){
  
  for(i in 1:length(colnames(data))) {
    
    if(class(data[[i]]) == "factor" | class(data[[i]]) == "character") {
      
      # print(data[i])
      print(colnames(data[i]))
      print(table(data[[i]]))
    }
    
  }
}

check_data(immune_data[1:50])
check_data(immune_data[51:100])
check_data(immune_data[101:150])
check_data(immune_data[151:200])
check_data(immune_data[201:234])






