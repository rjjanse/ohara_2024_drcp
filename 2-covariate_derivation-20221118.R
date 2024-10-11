#----------------------------------------------------------------------------------------#
# Treatment prescription patterns for SGLT2i, GLP-1, and DPP4
# Roemer Janse - 10/11/2022
# Code for covariate derivation in SCREAM3
#----------------------------------------------------------------------------------------#

# 0. Set-up ----
# Load packages
pacman::p_load("dplyr"      # Data manipulation
)

# Clean up environment
gc(verbose = FALSE)

# Set working directionary
setwd("C:/Users/rjjanse/OneDrive - LUMC/Research/Projects/8. GDT_SGLT2i/Codes/SCREAM3/dataframes")

# Load cohort with relevant variables
load("cohort.Rdata") 

# Keep only relevant variables
cohort.base <- dplyr::select(cohort, lopnr, drug, dt_index)

# 1. Comorbidities ----
# Load data
load("kon.Rdata"); load("ovr.Rdata"); load("slv.Rdata")
 
# Combine data
diagnoses <- 
    rbind(# Only lopnrs from kon that are also in the cohort
          filter(kon, lopnr %in% cohort.base[["lopnr"]]),
          # Only lopnrs from ovr that are also in the cohort
          filter(ovr, lopnr %in% cohort.base[["lopnr"]]),
          # Only lopnrs from slv that are also in the cohort
          filter(slv, lopnr %in% cohort.base[["lopnr"]]))

# Remove singular diagnoses data frames
rm(kon, ovr, slv)

# Create function for a diagnostic code with specifiable lookback (in years)
fun.com <- function(icd, lookback = NULL){
    # Lookback is 200 years if nothing specified, else as specified
    lb <- ifelse(is.null(lookback), 365.25 * 200, 365.25 * lookback)
    
    # Derive comorbidity
    dat.tmp <- cohort.base %>%
        # Join diagnoses data
        left_join(diagnoses, "lopnr") %>%
        # Create variable for comorbidity presence
        mutate(com = ifelse(grepl(icd, diagnosis) & !is.na(bdat) & as.Date(bdat) <= dt_index & as.Date(bdat) >= (dt_index - lb), 1, 0)) %>%
        # Sort on lopnr & drug and then descending on comorbidity (putting presence on top)
        arrange(lopnr, drug, desc(com)) %>%
        # Group on person and drug
        group_by(lopnr, drug) %>%
        # Keep first row per person
        slice(1L) %>%
        # Ungroup again
        ungroup() %>%
        # Keep only relevant variables
        dplyr::select(lopnr, drug, com)
    
    # Return data
    return(dat.tmp)
}

## 1.1. Acute coronary syndrome ----
dat.acs <- fun.com("^I200|^I21|^I22") %>% rename(com_acs = com)

## 1.2. Other ischemic heart disease ----
dat.ihd <- fun.com("^I201|^I208|^I209|^I24|^I25") %>% rename(com_ihd = com)

## 1.3. Hypertension ----
dat.hyp <- fun.com("^I1[0-5]") %>% rename(com_hyp = com)

## 1.4. Diabetic complications ----
dat.dmc <- fun.com("^E11[01234568]|^E13[01234568]|^E14[0-8]|^E16[012]|^G990|^G590|^G632|^H280|^H358|^H360|^M142|^M146") %>% rename(com_dmc = com)

## 1.5. Heart failure ----
dat.hfl <- fun.com("^I110|^I130|^I132|^I50") %>% rename(com_hfl = com)

## 1.6. Valve disorders ----
dat.vds <- fun.com("^I3[4-7]") %>% rename(com_vds = com)

## 1.7. Stroke ----
dat.str <- fun.com("^I6[0-4]|^I693|^I698|^I694") %>% rename(com_str = com)

## 1.8. Other cerebrovascular disease ----
dat.cbv <- fun.com("^I6[5-9]|^G45[012389]|^G46") %>% rename(com_cbv = com)

## 1.9. Atrial fibrillation ----
dat.afb <- fun.com("^I48") %>% rename(com_afb = com)

## 1.10. Other arrhythmia ----
dat.arr <- fun.com("^I4[4-7]|^I49") %>% rename(com_arr = com)

## 1.11. Peirpheral vascular disease ----
dat.pvd <- fun.com("^I70|^I72|^I73") %>% rename(com_pvd = com)

## 1.12. Chronic obstructive pulmonary disease ----
dat.cpd <- fun.com("J44") %>% rename(com_cpd = com)

## 1.13. Other lung disease ----
dat.old <- fun.com("^I27[0289]|^J4[0123567]|^J6[0-9]|^J70|^J84|^J92|^J96|^J982|^J983") %>% rename(com_old = com)

## 1.14. Venous thromboembolism ----
dat.vte <- fun.com("^I26|^I80[12389]|^I81|^I820|^I822|^I823|^I828|^I829") %>% rename(com_vte = com)

## 1.15. Cancer (excluding non-melanoma skin cancer) in the previous year ----
dat.can <- fun.com("^C0[1-9]|^C1[0-9]|^C2[0-6]|^C3[01234789]|^C4[0156789|^C5[0-8]|^C6[0-9]|C7[0-6]|^C8[1234568]|^C9[0-7]", lookback = 1) %>% rename(com_can = com)

## 1.16. Liver disease ----
dat.lvd <- fun.com("^B18|^I850|^I859|^I982|^K7[0-7]") %>% rename(com_lvd = com)

## 1.17. Fracture in the previous year ----
dat.frc <- fun.com("^S02[012346789]|^S12|^S22|^S32|^S42|^S52|^S62|^S72|^S82|^S92|^T02|^T08|T10|^T12|^M484|^M485|^M843", lookback = 1) %>% rename(com_frc = com)

## 1.18. Join all comorbidities together ----
# Create base data frame
dat.com <- dplyr::select(cohort.base, lopnr, drug)

for(i in c("dat.acs", "dat.ihd", "dat.hyp", "dat.dmc", "dat.hfl", "dat.vds", "dat.str", "dat.cbv", "dat.afb", 
           "dat.arr", "dat.pvd", "dat.cpd", "dat.old", "dat.vte", "dat.can", "dat.lvd", "dat.frc")){
    # Join current object to previous data by lopnr and drug
    dat.com <- left_join(dat.com, get(i), c("lopnr", "drug"))
}

# Calculate cardiovascular disease and atherosclerotic CVD comorbidity
dat.com <- dat.com %>%
    # Calculate variable
    mutate(com_cvd = ifelse(com_acs + com_ihd + com_hfl + com_str + com_cbv + com_pvd >= 1, 1, 0),
           com_acv = ifelse(com_acs + com_ihd + com_str + com_cbv + com_pvd >= 1, 1, 0))
               
# Remove comorbidity data frames
rm(dat.acs, dat.ihd, dat.hyp, dat.dmc, dat.hfl, dat.vds, dat.str, dat.cbv, dat.afb, 
   dat.arr, dat.pvd, dat.cpd, dat.old, dat.vte, dat.can, dat.lvd, dat.frc, diagnoses)

# Save data
save(dat.com, file = "dat_com.Rdata")

# 2. Medication ----
# Load drugs data
load("lmed_1.Rdata"); load("lmed_2.Rdata"); load("lmed_3.Rdata"); load("lmed_4.Rdata")

# Combine data
medication <- 
    rbind(# Only lopnrs from lmed_1 that are also in the cohort
          filter(lmed_1, lopnr %in% cohort.base[["lopnr"]]),
          # Only lopnrs from lmed_2 that are also in the cohort
          filter(lmed_2, lopnr %in% cohort.base[["lopnr"]]),
          # Only lopnrs from lmed_3 that are also in the cohort
          filter(lmed_3, lopnr %in% cohort.base[["lopnr"]]),
          # Only lopnrs from lmed_4 that are also in the cohort
          filter(lmed_4, lopnr %in% cohort.base[["lopnr"]]))

# Remove singular lmed data frames
rm(lmed_1, lmed_2, lmed_3, lmed_4)

# Load all dispensations without other information
load("disps.Rdata")

# Create function for a drug code in the past half year
fun.med <- function(atcc, any_dm = FALSE){
    # If we are not looking at any diabetes drugs, we can look including index date
    if(!any_dm){
        # Derive drug
        dat.tmp <- cohort.base %>%
            # Join diagnoses data
            left_join(medication, "lopnr") %>%
            # Create variable for drug presence
            mutate(med = ifelse(grepl(atcc, atc) & !is.na(edatum) & as.Date(edatum) <= dt_index & as.Date(edatum) >= (dt_index - 180), 1, 0)) %>%
            # Sort on lopnr & drug and then descending on drug (putting presence on top)
            arrange(lopnr, drug, desc(med)) %>%
            # Group on person and drug
            group_by(lopnr, drug) %>%
            # Keep first row per person
            slice(1L) %>%
            # Ungroup again
            ungroup() %>%
            # Keep only relevant variables
            dplyr::select(lopnr, drug, med)
    }
    
    # If we are looking at any diabetes drug, we must excluse index date
    if(any_dm){
        # Derive drug
        dat.tmp <- cohort.base %>%
            # Join diagnoses data
            left_join(medication, "lopnr") %>%
            # Create variable for drug presence
            mutate(med = ifelse(grepl(atcc, atc) & !is.na(edatum) & as.Date(edatum) < dt_index & as.Date(edatum) >= (dt_index - 180), 1, 0)) %>%
            # Sort on lopnr & drug and then descending on drug (putting presence on top)
            arrange(lopnr, drug, desc(med)) %>%
            # Group on person and drug
            group_by(lopnr, drug) %>%
            # Keep first row per person
            slice(1L) %>%
            # Ungroup again
            ungroup() %>%
            # Keep only relevant variables
            dplyr::select(lopnr, drug, med)
    }
 
    # Return data
    return(dat.tmp)
}

## 2.1. Beta-blockers ----
dat.bbs <- fun.med("^C07") %>% rename(med_bbs = med)

## 2.2. Calcium channel blockers ----
dat.ccb <- fun.med("^C08") %>% rename(med_ccb = med)

## 2.3. Diuretics ----
dat.dur <- fun.med("^C03") %>% rename(med_dur = med)

## 2.4. RASi ----
dat.ras <- fun.med("^C09A|^C09B|^C09C|^C09D") %>% rename(med_ras = med)

## 2.5. Statins ----
dat.sta <- fun.med("^C10") %>% rename(med_sta = med)

## 2.6. Digoxin ----
dat.dgx <- fun.med("^C01AA05") %>% rename(med_dgx = med)

## 2.7. Nitrates ----
dat.nit <- fun.med("^C01DA") %>% rename(med_nit = med)

## 2.8. Antiplatelets ----
dat.aps <- fun.med("^B01AC") %>% rename(med_aps = med)

## 2.9. Anticoagulants ----
dat.acs <- fun.med("^B01AA|^B01AE07|^B01AF|^B01AX05") %>% rename(med_acs = med)

## 2.10. B2-agonists (inhalant) ----
dat.bai <- fun.med("^R03AC") %>% rename(med_bai = med)

## 2.11. Anticholinergic inhalants ----
dat.aci <- fun.med("^R03BB") %>% rename(med_aci = med)

## 2.12. Glucocorticoid inhalants ----
dat.gci <- fun.med("^R03BA|^R03AK") %>% rename(med_gci = med)

## 2.13. Glucocorticoid oral ----
dat.gco <- fun.med("^H02AB") %>% rename(med_gco = med)

## 2.14. NSAIDs ----
dat.nsd <- fun.med("^M01A") %>% rename(med_nsd = med)

## 2.15. Opioids ----
dat.ops <- fun.med("^N02A") %>% rename(med_ops = med)

## 2.16. Any diabetes drug in the last 6 months ----
dat.add <- fun.med("^A10", any_dm = TRUE) %>% rename(med_add = med)

## 2.17. Metformin ----
dat.met <- fun.med("^A10BA02|^A10BD0[23578]|^A10BD10|^A10BD11|^A10BD13|^A10BD14|^A10BD15|^A10BD16|^A10BD20") %>% rename(med_met = med)

## 2.18. Sulfonylureas ----
dat.sus <- fun.med("^A10BB|^A10BD01|^A10BD02|^A10BD04|^A10BD06") %>% rename(med_sus = med)

## 2.19. GLP1-RAs ----
dat.glp <- fun.med("A10BJ|^A10BJ0[12356]|^A10AE56") %>% rename(med_glp = med)

## 2.20. SGLT2is ----
dat.sgl <- fun.med("A10BK|^A10BK0[1234]|^A10BD15|^A10BD19|^A10BD2[0134]") %>% rename(med_sgl = med)

## 2.21. DPP4i ----
dat.dpp <- fun.med("A10BH|^A10BH0[1235]|^A10BD0[78]|^A10BD1[01]") %>% rename(med_dpp = med)

## 2.22. Insulin ----
dat.ins <- fun.med("^A10AB|^A10AC|^A10AD|^A10AE") %>% rename(med_ins = med)

## 2.23. Other antidiabetics ----
dat.oad <- fun.med("^A10BF01|^A10BG|^A10BD03|^A10BD04|^A10BD05|^A10BD06|^A10BD09|^A10BD14|^A10BX") %>% rename(med_oad = med)

## 2.24. Time since first diabetes drug
dat.fdd <- cohort.base %>%
    # Join all diabetes drugs
    left_join(filter(medication, grepl("^A10", atc)), "lopnr") %>%
    # Keep only drugs before index
    filter(as.Date(edatum) <= dt_index) %>%
    # Sort on person, drug and drug date
    arrange(lopnr, drug, as.Date(edatum)) %>%
    # Group on person and drug
    group_by(lopnr, drug) %>%
    # Keep first row per person
    slice(1L) %>%
    # Ungroup again
    ungroup() %>%
    # Calculate new variables
    mutate(fdd_yrs = round(as.numeric(dt_index - as.Date(edatum)) / 365.25),                                  # Years since first diabetes drug
           fdd_cat = cut(fdd_yrs, c(0, 1, 3, 5, 7, 10000), labels = c(1:5)),                                  # Categorize years
           fdd_cat = factor(fdd_cat, levels = c(1:5), labels = c("<1", "1-3", ">=3-5", ">=5-7", ">=7"))) %>%  # Change fdd category to factor
    # Keep only relevant variables
    dplyr::select(lopnr, drug, fdd_yrs, fdd_cat)

## 2.25. Total number of drugs used ----
dat.tot <- cohort.base %>%
    # Join dispensations
    left_join(disps, "lopnr") %>%
    # Keep only dispensations in the past year
    filter(as.Date(edatum) >= (dt_index - 365.25) & as.Date(edatum) <= dt_index) %>%
    # Sort by person and drug
    arrange(lopnr, drug) %>%
    # Group per person and drug
    group_by(lopnr, drug) %>%
    # Count amount of drugs
    mutate(med_tot = length(unique(atc))) %>%
    # Keep one row per person & drug
    slice(1L) %>%
    # Ungroup again
    ungroup() %>%
    # Keep only relevant variables
    dplyr::select(lopnr, drug, med_tot) %>%
    # Calculate categories
    mutate(med_toc = case_when(med_tot <= 5 ~ "<=5",                                   # Define categories
                               med_tot > 5 & med_tot <= 10 ~ ">5-<=10",
                               med_tot > 10 & med_tot <= 15 ~ ">10-<=15",
                               med_tot > 15 ~ ">15"),
           med_toc = factor(med_toc, levels = c("<=5", ">5-<=10", ">10-<=15", ">15"))) # Factorize categories

## 2.26. Join all comorbidities together ----
# Create base data frame
dat.med <- dplyr::select(cohort.base, lopnr, drug)

for(i in c("dat.bbs", "dat.ccb", "dat.dur", "dat.ras", "dat.sta", "dat.dgx", "dat.nit", "dat.aps", "dat.acs", "dat.bai", "dat.aci",
           "dat.gci", "dat.gco", "dat.nsd", "dat.ops", "dat.add", "dat.met", "dat.sus", "dat.glp", "dat.sgl", "dat.dpp", "dat.ins", 
           "dat.oad", "dat.tot", "dat.fdd")){
    # Join current object to previous data by lopnr and drug
    dat.med <- left_join(dat.med, get(i), c("lopnr", "drug"))
}

# Remove comorbidity data frames
rm(dat.bbs, dat.ccb, dat.dur, dat.ras, dat.sta, dat.dgx, dat.nit, dat.aps, dat.acs, dat.bai, dat.aci, dat.tot, disps, medication,
   dat.gci, dat.nsd, dat.ops, dat.add, dat.met, dat.sus, dat.glp, dat.sgl, dat.dpp, dat.ins, dat.oad, dat.fdd, dat.gco)

# Save data
save(dat.med, file = "dat_med.Rdata")

# 3. Healthcare access in the past year ----
# Load healthcare access data
load("slv.Rdata"); load("kon.Rdata"); load("ovr.Rdata")

# Combine kon and ovr for outpatient
otp <- rbind(kon, ovr); rm(kon, ovr)

# Create function to derive healthcare access
fun.hca <- function(icd, df){
    # Derive data
    dat.hca <- cohort.base %>%
        # Join healthcare data
        left_join(df %>%
                      # Only first diagnosis
                      filter(diag_no == "diag1"),
                  "lopnr") %>%
        # Check if code of interest is present
        mutate(hca = ifelse(grepl(icd, diagnosis) & as.Date(bdat) <= dt_index & as.Date(bdat) >= (dt_index - 365.25) & !is.na(bdat), 1, 0)) %>%
        # Sort on identifier and drug and place presence of healthcare access on top
        arrange(lopnr, drug, desc(hca)) %>%
        # Group on identifier and drug
        group_by(lopnr, drug) %>%
        # Keep one row per person & drug
        slice(1L) %>%
        # Ungroup again
        ungroup() %>% 
        # Keep only relevant variables
        dplyr::select(lopnr, drug, hca)
    
    # Return data
    return(dat.hca)
}
        
## 3.1. Hospitalization due to cardiovascular causes ----
dat.hcv <- fun.hca("^I|^G45|^G46|^H341", slv) %>% rename(hca_hcv = hca)

## 3.2. Hospitalization due to diabetes mellitus ----
dat.hdm <- fun.hca("^E11", slv) %>% rename(hca_hdm = hca)

## 3.3. Hospitalization due to other reasons
dat.hot <- cohort.base %>%
    # Join hospital data
    left_join(slv %>%
                  # Keep only first diagnoses
                  filter(diag_no == "diag1"),
              "lopnr") %>%
    # Create new variables
    mutate(valid = ifelse(as.Date(bdat) <= dt_index & as.Date(bdat) >= (dt_index - 365.25) & !is.na(bdat), 1, 0),
           hca_hot = ifelse(valid == 1 & !(grepl("^I|^G45|^G46|^H341|^E11", diagnosis)), 1, 0)) %>%
    # Sort on identifier and drug and put hca.hot == 1 on top
    arrange(lopnr, drug, desc(hca_hot)) %>%
    # Group on identifer and drug
    group_by(lopnr, drug) %>%
    # Keep one row per person and drug
    slice(1L) %>%
    # Ungroup again
    ungroup() %>%
    # Keep only relevant variables
    dplyr::select(lopnr, drug, hca_hot)

## 3.4. Hospitalization due to cardiovascular causes ----
dat.ocv <- fun.hca("^I|^G45|^G46|^H341", otp) %>% rename(hca_ocv = hca)

## 3.5. Hospitalization due to diabetes mellitus ----
dat.odm <- fun.hca("^E11", otp) %>% rename(hca_odm = hca)

## 3.6. Hospitalization due to other reasons
dat.oot <- cohort.base %>%
    # Join hospital data
    left_join(otp %>%
                  # Keep only first diagnoses
                  filter(diag_no == "diag1"),
              "lopnr") %>%
    # Create new variables
    mutate(valid = ifelse(as.Date(bdat) <= dt_index & as.Date(bdat) >= (dt_index - 365.25) & !is.na(bdat), 1, 0),
           hca_oot = ifelse(valid == 1 & !(grepl("^I|^G45|^G46|^H341|^E11", diagnosis)), 1, 0)) %>%
    # Sort on identifier and drug and put hca.hot == 1 on top
    arrange(lopnr, drug, desc(hca_oot)) %>%
    # Group on identifer and drug
    group_by(lopnr, drug) %>%
    # Keep one row per person and drug
    slice(1L) %>%
    # Ungroup again
    ungroup() %>%
    # Keep only relevant variables
    dplyr::select(lopnr, drug, hca_oot)

## 3.7. Add everything together
# Create base data frame
dat.hca <- dplyr::select(cohort.base, lopnr, drug)

for(i in c("dat.hcv", "dat.hdm", "dat.hot", "dat.ocv", "dat.odm", "dat.oot")){
    # Join current object to previous data by lopnr and drug
    dat.hca <- left_join(dat.hca, get(i), c("lopnr", "drug"))
}

# Remove comorbidity data frames
rm(dat.hcv, dat.hdm, dat.hot, dat.ocv, dat.odm, dat.oot, slv, otp)

# Save data
save(dat.hca, file = "dat_hca.Rdata")

# 4. Socioeconomic data ----
# Load socioeconomic data
load("ses.Rdata")

# Calculate variables
dat.ses <- cohort.base %>%
    # Join socioeconomic data
    left_join(ses, "lopnr") %>%
    # Change year in SES data to halfway the year
    mutate(dt_ses = as.Date(paste0(year, "-07-01"), origin = "1970-01-01")) %>%
    # Keep only data before or at the index date
    filter(dt_ses <= dt_index) %>%
    # Sort on person, drug, and descending date (latest above)
    arrange(lopnr, drug, desc(dt_ses)) %>%
    # Group on person and drug
    group_by(lopnr, drug) %>%
    # Keep only most recent date per person 
    slice(1L) %>% 
    # Ungroup again
    ungroup() %>%
    # Sort on drug
    arrange(drug) %>%
    # Group by drug
    group_by(drug) %>%
    # Calculate the median income per group
    mutate(income_median = median(scb_dispincome)) %>%
    # Ungroup again
    ungroup() %>%
    # Create income variable
    mutate(ses_inc = ifelse(scb_dispincome >= income_median, "above median", "below median")) %>%
    # Rename family status and education
    rename(ses_fam = scb_fam, ses_edu = scb_education) %>%
    # Keep only relevant variables
    dplyr::select(lopnr, drug, ses_inc, ses_fam, ses_edu)

# As SES did not include everyone, join again to cohort base
dat.ses <- cohort.base %>%
    # Join data
    left_join(dat.ses, c("lopnr", "drug")) %>%
    # Add NA categories to ses variables
    mutate(ses_inc = ifelse(is.na(ses_inc), "missing", ses_inc),
           ses_fam = ifelse(is.na(ses_fam), "missing", ses_fam),
           ses_edu = ifelse(is.na(ses_edu), "missing", ses_edu)) %>%
    # Factorize SES variables
    mutate(ses_inc = factor(ses_inc, levels = c("below median", "above median", "missing")),
           ses_fam = factor(ses_fam, levels = c("Living alone", "Cohabitating", "missing")),
           ses_edu = factor(ses_edu, levels = c("Compulsory school", "Secondary school", "University", "missing"))) %>%
    # Remove index date
    dplyr::select(-dt_index)

# Save data
save(dat.ses, file = "dat_ses.Rdata")

# Remove data
rm(ses)

# 5. Laboratory data ----
## 5.1. Albumin-creatinine ratio ----
# Load albuminuria and dipstick data
load("uacr.Rdata"); load("dip.Rdata")

# Get ACR data
dat.acr <- cohort.base %>%
    # Join ACR data
    left_join(uacr, "lopnr") %>%
    # Keep only ACRs in the last year
    filter(as.Date(datum) >= (dt_index - 365.25) & as.Date(datum) <= dt_index & !is.na(datum)) %>%
    # Arrange on person, drug, and most recent measurement
    arrange(lopnr, drug, desc(as.Date(datum))) %>%
    # Group on person and drug
    group_by(lopnr, drug) %>%
    # Keep one row per person and drug
    slice(1L) %>%
    # Ungroup again
    ungroup() %>%
    # Create categories for ACR
    mutate(lab_acc = case_when(resultat < 3 ~ "A1",
                               resultat >= 3 & resultat <= 30 ~ "A2",
                               resultat > 30 ~ "A3")) %>%
    # Rename acr
    rename(lab_acr = resultat) %>%
    # Keep only relevant variables
    dplyr::select(lopnr, drug, lab_acr, lab_acc) %>%
    # Create indicator for source of data
    mutate(src = "lab")

# Get dipstick data for enrichment of ACR categories
dat.dip <- cohort.base %>%
    # Join dipstick data
    left_join(dip, "lopnr") %>%
    # Keep only dipsticks in the last year
    filter(as.Date(datum) >= (dt_index - 365.25) & as.Date(datum) <= dt_index & !is.na(datum)) %>%
    # Arrange on person, drug, and most recent measurement
    arrange(lopnr, drug, desc(as.Date(datum))) %>%
    # Group on person and drug
    group_by(lopnr, drug) %>%
    # Keep one row per person and drug
    slice(1L) %>%
    # Ungroup again
    ungroup() %>%
    # Create variables
    mutate(lab_acc = case_when(resultat == 1 | resultat == 2 ~ "A1",    # ACR categories based on dipstick
                               resultat == 2 ~ "A2",
                               resultat == 3 ~ "A3"),
           lab_acr = NA) %>%                                            # Empty variable for ACR to allow rbinding of data
    # Keep only relevant variables
    dplyr::select(lopnr, drug, lab_acr, lab_acc) %>%
    # Create indicator for source of data
    mutate(src = "dip")

# Combine data 
dat.acr <-
    # Bind data together
    rbind(
        # ACR data
        dat.acr,
        # Dipstick data
        dat.dip) %>%
    # Arrange on person, drug and descending on source (putting lab in the first row)
    arrange(lopnr, drug, desc(src)) %>%
    # Group on person an ddrug
    group_by(lopnr, drug) %>%
    # Keep only first row per person
    slice(1L) %>%
    # Ungroup again
    ungroup() %>%
    # Remove source column
    dplyr::select(-src)

# As some people might have been filtered out, add measurements to base cohort again
dat.acr <- cohort.base %>%
    # Join ACR data
    left_join(dat.acr, c("lopnr", "drug")) %>%
    # Create missing indicator for categories and then factorize categories
    mutate(lab_acc = ifelse(is.na(lab_acc), "missing", lab_acc),
           lab_acc = factor(lab_acc, levels = c("A1", "A2", "A3", "missing")))

## 5.2. HbA1c ----
# Load data
load("hba1c.Rdata")

# Derive HbA1c in past year
dat.hbc <- cohort.base %>%
    # Join data
    left_join(hba1c, "lopnr") %>%
    # Convert % to mmol/mol
    mutate(resultat = ifelse(enhet != "mmol/mol", 10.93 * resultat - 23.5, resultat)) %>%
    # Set results and date outside of the period of interest (year before index) to NA
    mutate(resultat = ifelse(as.Date(datum) >= (dt_index - 365.25) & as.Date(datum) <= dt_index & !is.na(datum), resultat, NA),
           datum = as.Date(ifelse(is.na(resultat), NA, datum), origin = "1970-01-01")) %>%
    # Arrange on person, drug
    arrange(lopnr, drug) %>%
    # Group by person and drug
    group_by(lopnr, drug) %>%
    # Per person calculate average lab value
    mutate(lab_hbc = mean(resultat, na.rm = TRUE)) %>%
    # If average is NaN, change to NA
    mutate(lab_hbc = ifelse(is.nan(lab_hbc), NA, lab_hbc)) %>%
    # Keep first row per person
    slice(1L) %>%
    # Ungroup again 
    ungroup() %>%
    # Keep only relevant variables
    dplyr::select(lopnr, drug, lab_hbc) %>%
    # Create new variables
    mutate(lab_hcc = case_when(lab_hbc < 53 ~ "<53",                                           # HbA1c categories
                               lab_hbc >= 53 & lab_hbc < 58 ~ ">=53-<58",
                               lab_hbc >= 58 & lab_hbc < 64 ~ ">=58-<64",
                               lab_hbc >= 64 & lab_hbc < 69 ~ ">=64-<69",
                               lab_hbc >= 69 & lab_hbc < 75 ~ ">=69-<75",
                               lab_hbc >= 75 ~ ">=75"),
           lab_hcc = ifelse(is.na(lab_hcc), "missing", lab_hcc),                               # Missing indicator for HbA1c categories
           lab_hcc = factor(lab_hcc, levels = c("<53", ">=53-<58", ">=58-<64", ">=64-<69",     # Factorize HbA1c categories
                                                ">=69-<75", ">=75", "missing")))

## 5.3. Creatinine and eGFR ----
# Load data
load("creat.Rdata"); load("artcreat.Rdata"); load("vencreat.Rdata"); load("demo.Rdata")

# Prepare data
creas <-
    # Bind data
    rbind(
        # Normal creatinines
        creat,
        # Venous creatinines
        vencreat,
        # Arterial creatinines
        artcreat) %>%
    # Keep only individuals in cohort
    filter(lopnr %in% cohort.base[["lopnr"]]) %>%
    # Remove inpatinet
    filter(ip == 0) %>%
    # Join date of birth
    left_join(demo %>%
                  # Keep only date of birth
                  dplyr::select(lopnr, dob_s3), 
              "lopnr") %>%
    # Calculate age at lab measurement
    mutate(test_age = as.numeric((as.Date(datum) - as.Date(dob_s3)) / 365.25)) %>%
    # Keep only relevant variables
    dplyr::select(lopnr, datum, resultat, test_age)

# Remove data
rm(creat, artcreat, vencreat, demo)

# Create function to calculate eGFR
ckd_epi <- function(creatinine, age, female){
    # Determine k coefficient
    k <- ifelse(female == 1, 62, 80)
    
    # Determine alpha coefficient
    alpha <- ifelse(female == 1, -0.329, -0.411)

    # Calculate eGFR if participant is female
    egfr <- ifelse(female == 1, 
                   141 * (pmin(creatinine/k, 1) ^ alpha) * (pmax(creatinine/k, 1) ^ (-1.209)) * (0.993 ^ age) * 1.018,
                   141 * (pmin(creatinine/k, 1) ^ alpha) * (pmax(creatinine/k, 1) ^ (-1.209)) * (0.993 ^ age))
    
    # Return eGFR
    return(egfr)
}

# Calculate creatinine and eGFR
dat.crt <- cohort %>%
    # Keep only relevant variables
    dplyr::select(lopnr, drug, dt_index, female) %>%
    # Join data
    left_join(creas, "lopnr") %>%
    # Set results and date outside of the period of interest (year before index) to NA
    mutate(resultat = ifelse(as.Date(datum) >= (dt_index - 365.25) & as.Date(datum) <= dt_index & !is.na(datum), resultat, NA),
           datum = as.Date(ifelse(is.na(resultat), NA, datum), origin = "1970-01-01")) %>%
    # Arrange on person, drug
    arrange(lopnr, drug) %>%
    # Group by person and drug
    group_by(lopnr, drug) %>%
    # Per person calculate new variables
    mutate(egfr = ckd_epi(resultat, test_age, female),         # Calculate eGFR for each creatinine    
           lab_gfr = mean(egfr, na.rm = TRUE),                 # Calculate average eGFR per person per drug
           lab_crt = mean(resultat, na.rm = TRUE),             # Calculate average creatine per person per drug
           lab_gfr = ifelse(is.nan(lab_gfr), NA, lab_gfr),     # Change NaN to NA
           lab_crt = ifelse(is.nan(lab_crt), NA, lab_crt)) %>% # Change NaN to NA
    # Keep first row per person
    slice(1L) %>%
    # Ungroup again 
    ungroup() %>%
    # Keep only relevant variables
    dplyr::select(lopnr, drug, lab_gfr, lab_crt) %>%
    # Create new variables
    mutate(lab_gfc = case_when(lab_gfr >= 90 ~ ">=90",                                                                          # Categorize eGFR
                               lab_gfr < 90 & lab_gfr >= 60 ~ ">=60-<90",                   
                               lab_gfr < 60 & lab_gfr >= 45 ~ ">=45-<60",
                               lab_gfr < 45 & lab_gfr >= 30 ~ ">=30-<45",
                               lab_gfr < 30 & lab_gfr >= 15 ~ ">=15-<30",
                               lab_gfr <15 ~ "<15"),
           lab_gfc = ifelse(is.na(lab_gfc), "missing", lab_gfc),                                                                # Create missing indicator
           lab_gfc = factor(lab_gfc, levels = c(">=90", ">=60-<90", ">=45-<60", ">=30-<45", ">=15-<30", "<15", "missing")))     # Factorize eGFR categories
    
## 5.4. Join lab data together ----
dat.lab <- dat.acr %>%
    # Join HbA1c data
    left_join(dat.hbc, c("lopnr", "drug")) %>%
    # Join creat & eGFR data
    left_join(dat.crt, c("lopnr", "drug")) %>%
    # Remove index date
    dplyr::select(-dt_index)

# Save data
save(dat.lab, file = "dat_lab.Rdata")
    
# Remove data
rm(dat.crt, dat.acr, dat.hbc, creas, hba1c, dip, uacr)

# 6. Merging and finalizing covariates data ----
# Merge all covariates data
covariates <- cohort %>%
    # Keep only identifier variables
    dplyr::select(lopnr, drug, age, dt_index) %>%
    # Join all data
    left_join(dat.com, c("lopnr", "drug")) %>%   # Comorbidities
    left_join(dat.med, c("lopnr", "drug")) %>%   # Medication
    left_join(dat.hca, c("lopnr", "drug")) %>%   # Healthcare access
    left_join(dat.ses, c("lopnr", "drug")) %>%   # Socialeconomic information
    left_join(dat.lab, c("lopnr", "drug")) %>%   # Laboratory data
    # Create a set of variables to be used during modelling
    mutate(mod_nnd = med_bbs + med_ccb + med_dur + med_ras + med_sta + med_dgx + med_nit +      # Number of non-diabetic drugs
               med_aps + med_acs + med_bai + med_aci + med_gci + med_gco + med_nsd + med_ops,
           mod_gfr = cut(lab_gfr, breaks = seq(0, 180, 10), labels = 17:0),                     # eGFR decrease in steps of 10 mL/min/1.73m2
           mod_gfr = as.numeric(ifelse(is.na(mod_gfr), 18, mod_gfr)),                           # Add missing category to eGFR decrease and change to numeric
           mod_hbc = cut(lab_hbc, breaks = seq(0, 130, 10), labels = 0:12),                     # HbA1c increase in steps of 10 mmol/mol
           mod_hbc = as.numeric(ifelse(is.na(mod_hbc), 13, mod_hbc)),                           # Add missing category to HbA1c decrease and change to numeric
           mod_age = as.numeric(cut(age, seq(15, 105, 10), labels = 0:8))) %>%                  # Age increase in steps of 10 years
    # Categorize age and calculate calendar year             
    mutate(age_cat = case_when(age < 50 ~ "<50",                                                # Age categories
                               age >= 50 & age < 60 ~ "50-59",
                               age >= 60 & age < 69 ~ "60-69",
                               age >= 70 & age < 79 ~ "70-79",
                               age >= 80 ~ ">=80"),
           age_cat = factor(age_cat, levels = c("<50", "50-59", "60-69", "70-79", ">=80")),     # Factorize age categories
           year = format(dt_index, "%Y")) %>%                                                   # Get calendar year
    # Remove age and index date variables again
    dplyr::select(-c(age, dt_index))

# Save covariate data
save(covariates, file = "covariates.Rdata")
