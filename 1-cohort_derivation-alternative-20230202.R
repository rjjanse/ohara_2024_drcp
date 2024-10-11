#----------------------------------------------------------------------------------------#
# Treatment prescription patterns for SGLT2i, GLP-1, and DPP4
# Roemer Janse - 02/02/2023
# Alternative code for data derivation in SCREAM3 to get unique number of people per step
#----------------------------------------------------------------------------------------#

# 0. Set-up ----
# Load packages
pacman::p_load("dplyr",       # Data manipulation
               "readr"        # Reading data
)

# Set working directionary
setwd("C:/Users/rjjan/Onedrive/Bureaublad/8. GDT_SGLT2i/Codes/SCREAM3/dataframes/")

# 1. Select new-users of the exposure drugs ----
# Load data
load("sglt2i.Rdata")
load("glp1.Rdata")
load("dpp4i.Rdata")

# Create function to deal with any corrections
fun.crc <- function(df){
    # Perform correction
    dat.tmp <- df %>%
        # Remove drugs for which we don't know antal or antnum
        filter(!is.na(antal) & !is.na(antnum)) %>%
        # Arrange on lopnr, atc, date, and descending on antal
        arrange(lopnr, atc, edatum, desc(antal)) %>%
        # Group on drug per person
        group_by(lopnr, atc) %>%
        # Create new variables
        mutate(lead_antal = lead(antal),                                    # Antal in next row
               lead_antal = ifelse(is.na(lead_antal), 0, lead_antal),       # If last in row, set to 0 instead of NA
               antal = ifelse(lead_antal < 0, antal + lead_antal, antal),   # If the next antal is negative, it's a correction of the previous antal, we add those together
               edatum = as.Date(ifelse(antal <= 0 & !is.na(antal),          # Set drug date to missing for correction rows
                                       NA, edatum), origin = "1970-01-01")) %>%
        # Ungroup again
        ungroup() %>%
        # Keep only observations with available drug date
        filter(!is.na(edatum)) %>%
        # Remove unnecessary variables
        dplyr::select(-c(forpddd, otyp, ar, pop, edatum_c, lead_antal))
    
    # Return data
    return(dat.tmp)
}

# Correct data
sglt2i <- fun.crc(sglt2i)
glp1 <- fun.crc(glp1)
dpp4i <- fun.crc(dpp4i)

# Function for counting number of individuals
fun.noi <- function(df){
    # Return length of unique identifiers
    return(length(unique(df[["lopnr"]])))
}

# Function for counting total number of unique individuals among all data 
fun.tui <- function(){
    # Return length of unique identifiers in combined data
    return(length(unique(rbind(sglt2i, glp1, dpp4i)[["lopnr"]])))
}

# Total number of users
fun.noi(sglt2i) # n = 33,012
fun.noi(glp1) # n = 33,765
fun.noi(dpp4i) # n = 24,837
fun.tui() # n = 66,946

# Function to keep only new-users in relevant time period
fun.nus <- function(df){
    # Code
    dat.tmp <- df %>%
        # Change edatum to date class
        mutate(edatum = as.Date(edatum, origin = "1970-01-01")) %>%
        # Sort on identifier and drug date
        arrange(lopnr, edatum) %>%
        # Group per person
        group_by(lopnr) %>%
        # Keep only first observation per person
        slice(1L) %>%
        # Ungroup again
        ungroup() %>%
        # Keep only individuals with first drug date between 2015 and 2020
        filter(edatum >= as.Date("2015-01-01") & edatum <= as.Date("2020-12-31")) %>%
        # Rename drug date to index date
        rename(dt_index = edatum)
    
    # Return data
    return(dat.tmp)
}

# Keep only SGLT2i users from relevant time period
sglt2i <- fun.nus(sglt2i)
glp1 <- fun.nus(glp1)
dpp4i <- fun.nus(dpp4i)

# Total number of new-users in relevant time period
fun.noi(sglt2i) # n = 15,815
fun.noi(glp1) # n = 17,032
fun.noi(dpp4i) # n = 14,128
fun.tui() # n = 37,445

# 2. Apply eligibility criteria ----
## 2.1. Age >= 18 years ----
# Create function
fun.age <- function(df){
    # Load demographics data
    load("demo.Rdata")
    
    # Code
    dat.tmp <- df %>%
        # Add demographics data
        left_join(demo, "lopnr") %>%
        # Calculate age
        mutate(age = round(as.numeric(dt_index - as.Date(dob_s3)) / 365.25)) %>%
        # Keep only individuals >= 18 years old
        filter(age >= 18) %>%
        # Drop unnecessary variables
        dplyr::select(-c(female:dob_s3)) %>%
        # Rename variables
        rename(female = female_s3)
    
    # Remove demographics data
    rm(demo)
    
    # Return data
    return(dat.tmp)
}

# Apply criterium
sglt2i <- fun.age(sglt2i)
glp1 <- fun.age(glp1)
dpp4i <- fun.age(dpp4i)

# Count individuals
fun.noi(sglt2i) # n = 15,814
fun.noi(glp1) # n = 17,017
fun.noi(dpp4i) # n = 14,126
fun.tui() # n = 37,427

## 2.2. Residency at baseline ----
# Note: in this section, we immediately determine the censor date, which is either death, migration, or 2021-12-31

# Create function
fun.mgr <- function(df){
    # Load migration and death data
    load("migration.Rdata"); load("death.Rdata")
    
    # For death data, keep only one record
    death <- distinct(death, lopnr, .keep_all = TRUE) %>%
        # Keep only necessary variables
        dplyr::select(lopnr, dodsdat)
    
    # Step 1: determine what individuals have an emigration as last migration point before index
    dat.exc <- df %>%
        # Add migration data
        left_join(migration, "lopnr") %>%
        # Remove any individuals without hdat
        filter(!is.na(hdat)) %>%
        # Keep only migration points before index
        filter(as.Date(hdat) <= dt_index) %>%
        # Sort on individual and then descending migration date
        arrange(lopnr, desc(as.Date(hdat))) %>%
        # Group per person
        group_by(lopnr) %>%
        # Keep only first observation (i.e., last migration point)
        slice(1L) %>% 
        # Ungroup again
        ungroup() %>%
        # Keep individuals who migrated before index
        filter(hkod == "U")
    
    # Step 2: applying criterium
    dat.tmp <- df %>%
        # Remove individuals who are in dat.exc (not meeting migration criterium)
        filter(!(lopnr %in% dat.exc[["lopnr"]]))
    
    # Step 3: calculating censoring date per person
    dat.cen <- df %>%
        left_join(migration, "lopnr") %>%
        # Remove any individuals without hdat
        filter(!is.na(hdat)) %>%
        # Keep only migration points after index
        filter(as.Date(hdat) > dt_index) %>%
        # Keep only migration points that indicate emigration
        filter(hkod == "U") %>%
        # Sort on individual and ascending migration date
        arrange(lopnr, as.Date(hdat)) %>%
        # Group per person
        group_by(lopnr) %>%
        # Keep only first observation (i.e., first emigration after index)
        slice(1L) %>%
        # Ungroup again
        ungroup() %>%
        # Keep only lopnr and emigration date
        dplyr::select(lopnr, hdat)
    
    # Step 4: calculate censoring date
    dat.cen <- dat.tmp %>%
        # Keep only necessary variables
        dplyr::select(lopnr, dt_index) %>%
        # Add death date
        left_join(death, "lopnr") %>%
        # Add migration date
        left_join(dat.cen, "lopnr") %>%
        # Compute variables
        mutate(dt_censor = pmin(as.Date(hdat), as.Date(dodsdat), as.Date("2021-12-31"), na.rm = TRUE),                       # Censor date
               censor_reason = ifelse(!is.na(dodsdat), "death", ifelse(!is.na(hdat), "migration", "administrative"))) %>%    # Reason for censoring
        # Keep only relevant variables
        dplyr::select(lopnr, dt_censor, censor_reason)
    
    # Step 5: add cenesoring dates to data
    dat.tmp <- dat.tmp %>%
        # Add censoring dates
        left_join(dat.cen, "lopnr")
    
    # Return data
    return(dat.tmp)
}

# Apply criterium
sglt2i <- fun.mgr(sglt2i)
glp1 <- fun.mgr(glp1)
dpp4i <- fun.mgr(dpp4i)

# Count individuals
fun.noi(sglt2i) # n = 14,435
fun.noi(glp1) # n = 15,754
fun.noi(dpp4i) # n = 12,771
fun.tui() # n = 34,484

## 2.3. Gestational diabetes or diabetes mellitus type I ----
# Create function
fun.wdt <- function(df){
    # Load data
    load("kon.Rdata"); load("ovr.Rdata"); load("slv.Rdata")
    
    # Combine data and keep only relevant diagnoses
    diags <- rbind(filter(kon, grepl("^O244|^E10", diagnosis)),     # Gestational diabetes and DM type I from kon
                   filter(ovr, grepl("^O244|^E10", diagnosis)),     # Gestational diabetes and DM type I from ovr
                   filter(slv, grepl("^O244|^E10", diagnosis))) %>% # Gestational diabetes and DM type I from slv
        # Change diagnosis date from character to date class
        mutate(bdat = as.Date(bdat)) %>%
        # Arrange on person and diagnosis date
        arrange(lopnr, bdat)
    
    # Determine exclusions for diabetes type 1
    dat.exc.d1 <- df %>%
        # Join diagnoses data
        left_join(diags, "lopnr") %>%
        # Keep only diabetes type 1 codes
        filter(grepl("^E10", diagnosis)) %>%
        # Keep only diagnoses up until baseline
        filter(bdat <= dt_index) %>%
        # Keep one row per person
        distinct(lopnr, .keep_all = TRUE)
    
    # Determine exclusions for gestational diabetes (2 recorded diagnoses)
    dat.exc.gd <- df %>%
        # Join diagnoses data
        left_join(diags, "lopnr") %>%
        # Keep only gestational diabetes codes
        filter(grepl("^O244", diagnosis)) %>%
        # Keep only diagnoses up until baseline
        filter(bdat <= dt_index) %>%
        # Sort by person
        arrange(lopnr) %>%
        # Group per person
        group_by(lopnr) %>%
        # Keep only people who have two or more diagnoses
        filter(max(row_number()) >= 2) %>%
        # Keep one row per person
        distinct(lopnr, .keep_all = TRUE)
    
    # Keep only individuals not in dat.exc.d1 or dat.exc.gd
    dat.tmp <- df %>% 
        # Exclude persons
        filter(!(lopnr %in% dat.exc.d1[["lopnr"]]) & !(lopnr %in% dat.exc.gd[["lopnr"]]))
    
    # Return data
    return(dat.tmp)
}

# Apply criterium
sglt2i <- fun.wdt(sglt2i)
glp1 <- fun.wdt(glp1)
dpp4i <- fun.wdt(dpp4i)

# Count individuals
fun.noi(sglt2i) # n = 11,449
fun.noi(glp1) # n = 12,378
fun.noi(dpp4i) # n = 10,452
fun.tui() # n = 27,719

## 2.4. Prior kidney replacement therapy ----
# Create function
fun.krt <- function(df){
    # Load data
    load("krt.Rdata")
    
    # Determine exclusions
    dat.exc <- df %>%
        # Join krt data
        left_join(krt, "lopnr") %>%
        # Keep only KRT
        filter(rrt == 1) %>%
        # Keep only observations with KRT date up to index date
        filter(as.Date(rrt_date) <= dt_index) %>%
        # Keep one row per person
        distinct(lopnr, .keep_all = TRUE)
    
    # Keep only individuals that are not excluded
    dat.tmp <- df %>%
        # Exclude people
        filter(!(lopnr %in% dat.exc[["lopnr"]]))
    
    # Return data
    return(dat.tmp)
}

# Apply criterium
sglt2i <- fun.krt(sglt2i)
glp1 <- fun.krt(glp1)
dpp4i <- fun.krt(dpp4i)

# Count individuals
fun.noi(sglt2i) # n = 11,440
fun.noi(glp1) # n = 12,345
fun.noi(dpp4i) # n = 10,358
fun.tui() # n = 27,597

## 2.5. Use of liraglutide with obesity indication (Saxenda, product number 513490) ----
# Load raw GLP1-RA data
load("glp1_raw.Rdata")

# Derive all individuals using Saxenda
dat.exc.sax <- filter(glp1_raw, VARUNR == 513490)

# Exclude individuals & remove data
glp1 <- filter(glp1, !(lopnr %in% dat.exc.sax[["LopNr"]])); rm(glp1_raw, dat.exc.sax)

# Count individuals
fun.noi(glp1) # n = 11,096
fun.tui() # n = 26,362

## 2.6. End-stage disease ----
# Create function
fun.esd <- function(df){
    # Load diagnoses data
    load("kon.Rdata"); load("ovr.Rdata"); load("slv.Rdata")
    
    # Combine data and keep only relevant people
    diags <- 
        # Combine data
        rbind(kon, ovr, slv) %>%
        # Keep only lopnrs in data
        filter(lopnr %in% df[["lopnr"]])
    
    # Load dementia drugs
    load("dementia_drugs.Rdata")
    
    # Create exclusion based on diagnoses
    dat.exc.dgs <- df %>%
        # Add diagnoses data
        left_join(diags, "lopnr") %>%
        # Keep only relevant diagnoses
        filter(grepl("^E4[0-3]|^F0[0-3]|^G30|^R402|^R64", diagnosis)) %>%
        # Remove diagnoses issued over 10 years ago
        filter(as.Date(bdat) >= (dt_index - (365.25 * 10))) %>%
        # Remove diagnoses issued later than baseline
        filter(as.Date(bdat) <= dt_index) %>%
        # Keep one row per person
        distinct(lopnr, .keep_all = TRUE)
    
    # Create exclusion based on dementia drugs
    dat.exc.dem <- df %>%
        # Add dementia data
        left_join(dem, "lopnr") %>%
        # Remove diagnoses issued over 10 years ago
        filter(as.Date(edatum) >= (dt_index - (365.25 * 10))) %>%
        # Remove diagnoses issued later than baseline
        filter(as.Date(edatum) <= dt_index) %>%
        # Keep one row per person
        distinct(lopnr, .keep_all = TRUE)
    
    # Keep only individuals that are not excluded
    dat.tmp <- df %>% 
        # Exclude persons
        filter(!(lopnr %in% dat.exc.dgs[["lopnr"]]) & !(lopnr %in% dat.exc.dem[["lopnr"]]))
    
    # Return data
    return(dat.tmp)
}

# Apply criterium
sglt2i <- fun.esd(sglt2i)
glp1 <- fun.esd(glp1)
dpp4i <- fun.esd(dpp4i)

# Count individuals
fun.noi(sglt2i) # n = 11,327
fun.noi(glp1) # n = 11,001
fun.noi(dpp4i) # n = 10,060
fun.tui() # n = 25,917

## 2.7. Drug misuse ----
# Create function
fun.dmu <- function(df){
    # Load drug data
    # Medication data was too big to extract as one file, so it is split into four files.
    # lmed_1 contains N07BB and N07BC for SCREAM2 and lmed_4 contains these for SCREAM3
    load("lmed_1.Rdata"); load("lmed_4.Rdata")
    
    # Combine data and keep only relevant people
    med <- 
        # Combine data
        rbind(lmed_1, lmed_4) %>%
        # Keep only lopnrs in data
        filter(lopnr %in% df[["lopnr"]]) %>%
        # Rename atc column to atc_new as atc is already present in other data
        rename(atc_new = atc)
    
    # Load diagnoses data
    load("kon.Rdata"); load("ovr.Rdata"); load("slv.Rdata")
    
    # Combine data and keep only relevant people
    diags <- 
        # Combine data
        rbind(kon, ovr, slv) %>%
        # Keep only lopnrs in data
        filter(lopnr %in% df[["lopnr"]])
    
    # Exclusion based on drugs
    dat.exc.drg <- df %>%
        # Join medication data
        left_join(med, "lopnr") %>%
        # Keep only relevant drugs
        filter(grepl("^N07B[BC]", atc_new)) %>%
        # Keep only drugs up to baseline
        filter(as.Date(edatum) <= dt_index) %>%
        # Keep only drugs up to 10 years before baseline
        filter(as.Date(edatum) >= (dt_index - (365.25 * 10))) %>%
        # Keep only one row per person
        distinct(lopnr, .keep_all = TRUE)
    
    # Exclusion based on diagnoses
    dat.exc.dgs <- df %>%
        # Add diagnoses data
        left_join(diags, "lopnr") %>%
        # Keep only relevant diagnoses
        filter(grepl("^F1[12345689]|^R78[1-5]|^T40", diagnosis)) %>%
        # Remove diagnoses issued over 10 years ago
        filter(as.Date(bdat) >= (dt_index - (365.25 * 10))) %>%
        # Remove diagnoses issued later than baseline
        filter(as.Date(bdat) <= dt_index) %>%
        # Keep one row per person
        distinct(lopnr, .keep_all = TRUE)
    
    # Keep only individuals that are not excluded
    dat.tmp <- df %>% 
        # Exclude persons
        filter(!(lopnr %in% dat.exc.drg[["lopnr"]]) & !(lopnr %in% dat.exc.dgs[["lopnr"]]))
    
    # Return data
    return(dat.tmp)
}

# Apply criterium
sglt2i <- fun.dmu(sglt2i)
glp1 <- fun.dmu(glp1)
dpp4i <- fun.dmu(dpp4i)

# Count individuals
fun.noi(sglt2i) # n = 10,959
fun.noi(glp1) # n = 10,523
fun.noi(dpp4i) # n = 9,722
fun.tui() # n = 24,970

## 2.8. Severe pancreatic disorder ----
# Create function
fun.spd <- function(df){
    # Load diagnoses data
    load("kon.Rdata"); load("ovr.Rdata"); load("slv.Rdata")
    
    # Combine data and keep only relevant people
    diags <- 
        # Combine data
        rbind(kon, ovr, slv) %>%
        # Keep only lopnrs in data
        filter(lopnr %in% df[["lopnr"]])
    
    # Load procedure data
    load("ovr_procedurecodes.Rdata"); load("slv_procedurecodes.Rdata")
    
    # Combine data and keep only relevant people
    procs <- 
        # Combine data
        rbind(ovr_proc, slv_proc) %>%
        # Keep only lopnrs in data
        filter(lopnr %in% df[["lopnr"]])
    
    # Load drug data
    # Medication data was too big to extract as one file, so it is split into four files.
    # lmed_1 contains N07BB and N07BC for SCREAM2 and lmed_4 contains these for SCREAM3
    load("lmed_1.Rdata"); load("lmed_4.Rdata")
    
    # Combine data and keep only relevant people
    med <- 
        # Combine data
        rbind(lmed_1, lmed_4) %>%
        # Keep only lopnrs in data
        filter(lopnr %in% df[["lopnr"]]) %>%
        # Rename atc column to atc_new as atc is already present in other data
        rename(atc_new = atc)
    
    # Exclusion based on diagnoses
    dat.exc.dgs <- df %>%
        # Add diagnoses data
        left_join(diags, "lopnr") %>%
        # Keep only relevant diagnoses
        filter(grepl("^C25|^K86[01]", diagnosis)) %>%
        # Remove diagnoses issued over 10 years ago
        filter(as.Date(bdat) >= (dt_index - (365.25 * 10))) %>%
        # Remove diagnoses issued later than baseline
        filter(as.Date(bdat) <= dt_index) %>%
        # Keep one row per person
        distinct(lopnr, .keep_all = TRUE)
    
    # Exclusion based on procedures
    dat.exc.prc <- df %>%
        # Join procedure data
        left_join(procs, "lopnr") %>%
        # Keep only relevant drugs
        filter(grepl("^JL[CE]", opk)) %>%
        # Keep only drugs up to baseline
        filter(as.Date(datum) <= dt_index) %>%
        # Keep only drugs up to 10 years before baseline
        filter(as.Date(datum) >= (dt_index - (365.25 * 10))) %>%
        # Keep only one row per person
        distinct(lopnr, .keep_all = TRUE)
    
    # Exclusion based on drugs
    dat.exc.drg <- df %>%
        # Join medication data
        left_join(med, "lopnr") %>%
        # Keep only relevant drugs
        filter(grepl("^A09AA02", atc_new)) %>%
        # Keep only drugs up to baseline
        filter(as.Date(edatum) <= dt_index) %>%
        # Keep only drugs up to 10 years before baseline
        filter(as.Date(edatum) >= (dt_index - (365.25 * 10))) %>%
        # Keep only one row per person
        distinct(lopnr, .keep_all = TRUE)
    
    # Keep only individuals that are not excluded
    dat.tmp <- df %>% 
        # Exclude persons
        filter(!(lopnr %in% dat.exc.drg[["lopnr"]]) & !(lopnr %in% dat.exc.dgs[["lopnr"]]) & !(lopnr %in% dat.exc.prc[["lopnr"]]))
    
    # Return data
    return(dat.tmp)
}

# Apply criterium
sglt2i <- fun.spd(sglt2i)
glp1 <- fun.spd(glp1)
dpp4i <- fun.spd(dpp4i)

# Count individuals
fun.noi(sglt2i) # n = 10,745
fun.noi(glp1) # n = 10,315
fun.noi(dpp4i) # n = 9,488
fun.tui() # n = 24,472

## 2.9. No specialist care / prescription drugs in the last year ----
# Create function
fun.nsp <- function(df){
    # Load specialist data
    load("ovr.Rdata")
    
    # Load dispensation data
    load("disps.Rdata")
    
    # Derive individuals who had specialist care and should thus stay included
    dat.inc.spc <- df %>%
        # Join specialist data
        left_join(ovr, "lopnr") %>%
        # Keep only care within 1 year of baseline
        filter(as.Date(bdat) >= (dt_index - 365.25)) %>%
        # Keep only care up until baseline
        filter(as.Date(bdat) <= dt_index) %>%
        # Keep only one row per person
        distinct(lopnr, .keep_all = TRUE)
    
    # Derive individuals who had prescription drugs and should thus stay included
    dat.inc.pds <- df %>%
        # Join drug data
        left_join(disps, "lopnr") %>%
        # Keep only drugs within 1 year of baseline
        filter(as.Date(edatum) >= (dt_index - 365.25)) %>%
        # Keep only drugs up until baseline
        filter(as.Date(edatum) <= dt_index) %>%
        # Keep only one row per person
        distinct(lopnr, .keep_all = TRUE)
    
    # Keep only individuals that are in either one of the two inclusion data frames
    dat.tmp <- df %>%
        # Include individuals
        filter((lopnr %in% dat.inc.spc[["lopnr"]]) | (lopnr %in% dat.inc.pds[["lopnr"]]))
    
    # Return data
    return(dat.tmp)
}

# Apply criterium
sglt2i <- fun.nsp(sglt2i)
glp1 <- fun.nsp(glp1)
dpp4i <- fun.nsp(dpp4i)

# Count individuals
fun.noi(sglt2i) # n = 10,745
fun.noi(glp1) # n = 10,315
fun.noi(dpp4i) # n = 9,488
fun.tui() # n = 24,472

## 2.10. Prescription on day of death or erroneous prescription after death ----
# Create function
fun.err <- function(df){
    # Remove individuals
    dat.tmp <- filter(df, dt_censor > dt_index)
    
    # Return data
    return(dat.tmp)
}

# Apply criterium
sglt2i <- fun.err(sglt2i)
glp1 <- fun.err(glp1)
dpp4i <- fun.err(dpp4i)

# Count individuals
fun.noi(sglt2i) # n = 10,743
fun.noi(glp1) # n = 10,315
fun.noi(dpp4i) # n = 9,488
fun.tui() # n = 24,470

# 3. Combine all data ----
# Create final cohort data frame
cohort <- rbind(# Further modify SGLT2i data
    sglt2i %>%
        # Remove unnecessary variables
        dplyr::select(-c(antal:antnum)) %>%
        # Create drug indicator
        mutate(drug = "sglt2i"),
    # Further modify GLP1-RA data
    glp1 %>%
        # Remove unnecessary variables
        dplyr::select(-c(antal:antnum)) %>%
        # Create drug indicator
        mutate(drug = "glp1"),
    # Further modify DPP4i data
    dpp4i %>%
        # Remove unnecessary variables
        dplyr::select(-c(antal:antnum)) %>%
        # Create drug indicator
        mutate(drug = "dpp4i"))

# Save data
save(cohort, file = "cohort.Rdata")
