#----------------------------------------------------------------------------------------#
# Treatment prescription patterns for SGLT2i, GLP-1, and DPP4
# Roemer Janse - 10/11/2022
# Code for adherence and persistence calculation in SCREAM3
#----------------------------------------------------------------------------------------#

# 0. Set-up ----
# Load packages
pacman::p_load("dplyr",     # Data manipulation
               "readxl"     # Read .xlsx files
)

# Clean up environment
gc(verbose = FALSE); rm(list = ls())

# Set working directionary
setwd("C:/Users/rjjan/OneDrive/Bureaublad/SCREAM3/dataframes")

# Load cohort
load("cohort.Rdata")

# Get base cohort
cohort.base <- dplyr::select(cohort, lopnr, drug, dt_index, dt_censor) 

# 1. Prepare functions ----
## 1.1. Create function to deal with any corrections ----
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

## 1.2. Create function for adherence ----
fun.adh <- function(df_base, drg, df_drug, yr){
    # To get PDC:
    # First get days we followed a person
    # Then determine for each dispensation when it supposedly ended and how many days it covered
    # For end dates after the censor date, we determine the days between the start date and censor date
    # Now we have days covered and total days, so we can determine percentage and cut-off
    
    # Apply corrections in data
    dat.drg <- fun.crc(df_drug)
    
    # Prepare drug dataset
    dat.drg <- df_base %>%
        # Keep only drug of interest
        filter(drug == drg) %>% 
        # Create variables
        mutate(yr = yr,                                                                                            # Create column with years
               dt_year = as.Date(ifelse(yr == 0, dt_censor, (dt_index + (yr * 365.25))),  origin = "1970-01-01"),  # End date for adherence calculation
               dt_year = as.Date(ifelse(dt_year > dt_censor, dt_censor, dt_year), origin = "1970-01-01")) %>%      # Set to censor date if later than censoring
        # Join dispensation data
        left_join(dat.drg, "lopnr") %>%
        # Keep only drugs prescribed before or on the last day of observation
        filter(edatum <= dt_year) %>% 
        # Remove unnecessary variables
        dplyr::select(-(drug:dt_year)) %>% 
        # Change edatum to date
        mutate(edatum = as.Date(edatum, origin = "1970-01-01")) %>%
        # Sort on person and date of prescription
        arrange(lopnr, edatum)
    
    # Prepare cohort data
    dat.adh <- df_base %>%
        # Keep only drug of interest
        filter(drug == drg) %>% 
        # Create new variables
        mutate(yr = yr,                                                                                          # Create column with years
               dt_year = as.Date(ifelse(yr == 0, dt_censor, (dt_index + (yr * 365.25))), origin = "1970-01-01"), # End date for adherence calculation
               dt_year = as.Date(ifelse(dt_year > dt_censor, dt_censor, dt_year), origin = "1970-01-01"),        # Set to censor date if later than censoring 
               totday = as.numeric(dt_year - dt_index)) %>%                                                      # Total amount of days (denominator for proportion of days covered)
        # Join prepared drugs data
        left_join(dat.drg, "lopnr") %>%
        # Create base variables
        mutate(days = ifelse(drug == "glp1", (antal * dur), (antal * antnum)),     # Calculate amount of days the dispensation lasts
               dt_end = as.Date(edatum + days),                                    # Date where current dispensation ends
               missed_days = 0,                                                    # Indicator for days without pills
               dt_out_of_pills = dt_end) %>%                                       # Date that pills run out
        # Arrange on person
        arrange(lopnr) %>% 
        # Group per person
        group_by(lopnr) %>% 
        # Calculate new variables 
        mutate(rownr = row_number(),                                                                           # Create a row number per dispensation
               nrows = max(rownr),                                                                             # Amount of rows (i.e., dispensations) each person has
               dt_out_of_pills = as.Date(ifelse(rownr == 1, dt_out_of_pills, NA), origin = "1970-01-01")) %>%  # For now, we only know the out of pills date for the first dispensation
        # Ungroup again
        ungroup()
    
    # Take into account stockpiling
    for(i in 2:max(dat.adh[["nrows"]])){
        # Loop on each dispensation
        dat.adh <- dat.adh %>% 
            # Sort on person
            arrange(lopnr) %>% 
            # Group by person
            group_by(lopnr) %>% 
            # Create and update variables
            mutate(dt_out_of_pills_lag = lag(dt_out_of_pills),                                                               # Determine last out of pills date
                   missed_days = ifelse(rownr == i & edatum > dt_out_of_pills_lag, as.numeric(edatum - dt_out_of_pills_lag), # If new edatum is later than out of pills date (previous observation), 
                                        ifelse(rownr == i & edatum <= dt_out_of_pills_lag, 0, missed_days)),                 # calculate missed days, otherwise 0
                   dt_out_of_pills = as.Date(ifelse(rownr == i & missed_days == 0, dt_out_of_pills_lag + days,               # If no days missed, out of pills date is previous out of pills date + days
                                                    ifelse(rownr == i & missed_days != 0, edatum + days, dt_out_of_pills)),  # If days were missed, out of pills date starts anew from current dispensation
                                             origin = "1970-01-01")) %>%
            # Ungroup again
            ungroup()
    }

    # Make final calculations for adherence                                                  
    dat.adh <- dat.adh %>% 
        # Sort on person and put last dispensation on top
        arrange(lopnr, desc(rownr)) %>% 
        # Group per person
        group_by(lopnr) %>%
        # Calculate and update variables
        mutate(dt_out_of_pills = as.Date(ifelse(dt_out_of_pills > dt_year, dt_year, dt_out_of_pills), origin = "1970-01-01"),  # Limit out of pills date at end of observation date to prevent PDC > 100%
               md = sum(missed_days),                                                                                          # Calculate all missed days
               md = md + as.numeric(dt_year - dt_out_of_pills)) %>%                                                            # Add any days missed between last dispensation and end of observation
        # Keep one row per person
        slice(1L) %>% 
        # Ungroup again
        ungroup() %>%
        # Make final calculatoins
        mutate(dc = totday - md,                      # Days covered
               pdc = round(dc / totday * 100),        # Proportion of days covered
               pdc_ind = ifelse(pdc >= 80, 0, 1)) %>% # Indicator for PDC >= 80%
        # Select only relevant variables
        dplyr::select(lopnr, drug, pdc, pdc_ind)
    
    # Return data
    return(dat.adh)
}

## 1.3. Create function for persistence ----
fun.per <- function(df_base, drg, df_drug, yr, gp){
    # To get persistence (60 day gap)
    # First determine when a prescription ends
    # Then calculate the gap between each end date and subsequent start date
    
    # Apply corrections in data
    dat.drg <- fun.crc(df_drug)
    
    # Prepare drug dataset
    dat.drg <- df_base %>%
        # Keep only drug of interest
        filter(drug == drg) %>% 
        # Create variables
        mutate(yr = yr,                                                                                            # Create column with years
               dt_year = as.Date(ifelse(yr == 0, dt_censor, (dt_index + (yr * 365.25))),  origin = "1970-01-01"),  # End date for adherence calculation
               dt_year = as.Date(ifelse(dt_year > dt_censor, dt_censor, dt_year), origin = "1970-01-01")) %>%      # Set to censor date if later than censoring
        # Join dispensation data
        left_join(dat.drg, "lopnr") %>%
        # Keep only drugs prescribed before or on the last day of observation
        filter(edatum <= dt_year) %>% 
        # Remove unnecessary variables
        dplyr::select(-(drug:dt_year)) %>% 
        # Change edatum to date
        mutate(edatum = as.Date(edatum, origin = "1970-01-01")) %>%
        # Sort on person and date of prescription
        arrange(lopnr, edatum)
    
    # Prepare cohort data
    dat.per <- df_base %>%
        # Keep only drug of interest
        filter(drug == drg) %>% 
        # Create new variables
        mutate(yr = yr,                                                                                          # Create column with years
               dt_year = as.Date(ifelse(yr == 0, dt_censor, (dt_index + (yr * 365.25))), origin = "1970-01-01"), # End date for adherence calculation
               dt_year = as.Date(ifelse(dt_year > dt_censor, dt_censor, dt_year), origin = "1970-01-01"),        # Set to censor date if later than censoring 
               totday = as.numeric(dt_year - dt_index)) %>%                                                      # Total amount of days (denominator for proportion of days covered)
        # Join prepared drugs data
        left_join(dat.drg, "lopnr") %>%
        # Create base variables
        mutate(days = ifelse(drug == "glp1", (antal * dur), (antal * antnum)),     # Calculate amount of days the dispensation lasts
               dt_end = as.Date(edatum + days)) %>%                                # Date where current dispensation ends
        # Arrange on person and drug date
        arrange(lopnr, edatum) %>% 
        # Group per person
        group_by(lopnr) %>% 
        # Calculate new variables 
        mutate(rownr = row_number(),                                                              # Create a row number per dispensation
               nrows = max(rownr),                                                                # Amount of rows (i.e., dispensations) each person has
               dt_out_of_pills = as.Date(ifelse(rownr == 1, dt_end, NA), origin = "1970-01-01"),  # For now, we only know when the current dispensation has ended
               gap = NA) %>%                                                                      # For now, we don't know the gap between dispensations
        # Ungroup again
        ungroup()

    # Take into account stockpiling
    for(i in 1:max(dat.per[["rownr"]])){
        # Apply first calculations only to first row
        if(i == 1){
            # Calculate gaps
            dat.per <- dat.per %>% 
                # Arrange on person
                arrange(lopnr) %>% 
                # Group by person
                group_by(lopnr) %>%
                # Create new varaibles
                mutate(lead_edatum = lead(edatum),                                                                            # Determine date of next dispensation
                       lead_edatum = as.Date(ifelse(is.na(lead_edatum), dt_year, lead_edatum), origin = "1970-01-01"),        # If next dispensation is missing (i.e., only one dispensation ever),
                                                                                                                              # then next date is the end of observation
                       gap = ifelse(rownr == 1 & lead_edatum > dt_out_of_pills, as.numeric(lead_edatum - dt_out_of_pills),    # If the new drug is later than out of pills date, calculate gap between
                                    ifelse(rownr == 1 & lead_edatum <= dt_out_of_pills, 0, gap))) %>% 
                # Ungroup again
                ungroup()
        }
        
        # Now calculate other rows
        else{
            # Calculate gaps
            dat.per <- dat.per %>% 
                # Arrange on person
                arrange(lopnr) %>% 
                # Group by person
                group_by(lopnr) %>%
                # Create new variables
                mutate(lag_gap = lag(gap),                                                                                 # Get gap from previous dispensation
                       dt_out_of_pills_lag = lag(dt_out_of_pills),                                                         # Get date previous dispensation ran out of pills
                       dt_out_of_pills = as.Date(ifelse(rownr == i & lag_gap == 0, dt_out_of_pills_lag + days,             # If there was a gap, no matter how long, out of pills date is reset to end 
                                                        ifelse(rownr == i & lag_gap != 0, dt_end, dt_out_of_pills)),       # date of current dispensation, otherwise add days of current dispensation to 
                                                 origin = "1970-01-01"),                                                   # out of pills date
                       lead_edatum = lead(edatum),                                                                         # Determine when next dispensation is    
                       lead_edatum = as.Date(ifelse(is.na(lead_edatum), dt_year, lead_edatum), origin = "1970-01-01"),     # If no next dispensation, it's the last row: we look until end of observation
                       gap = ifelse(rownr == i & lead_edatum > dt_out_of_pills, as.numeric(lead_edatum - dt_out_of_pills), # If new edatum is later than out of pills date, calculate gap, otherwise 0
                                    ifelse(rownr == i & lead_edatum <= dt_out_of_pills, 0, gap))) %>% 
                # Ungroup again
                ungroup()
        }
    }
    
    # Make final calculations
    dat.per <- dat.per %>%
        # Create new variables
        mutate(gap_ind = ifelse(gap >= gp, 1, 0),                                                                   # Gap indicator
               gap_dt = as.Date(ifelse(gap_ind == 1, (dt_out_of_pills + gp), dt_year), origin = "1970-01-01"),      # Gap date for people with a gap is end date plus grace period, otherwise end date
               gap_ind = ifelse(gap_dt > dt_year & gap_ind == 1, 0, gap_ind),                                       # Gap indicators for people with grace period  after the year date is 0 (no gap)
               gap_dt = as.Date(ifelse(gap_dt > dt_year, dt_year, gap_dt), origin = "1970-01-01")) %>%              # End date for end dates past year date is year date. This also sets end dates with 
                                                                                                                    # grace period after a year for people who went from gap to no gap back to right date
        # Sort on person, put gaps on top and then sort on first gap
        arrange(lopnr, desc(gap_ind), gap_dt) %>% 
        # Group per person
        group_by(lopnr) %>% 
        # Keep first row per person
        slice(1L) %>% 
        # Ungropu again
        ungroup() %>%
        # Calculate time to gap
        mutate(time2gap = ifelse(gap_ind == 1, as.numeric(gap_dt - dt_index), as.numeric(dt_year - dt_index))) %>%
        # Keep only relevant varaibles
        dplyr::select(lopnr, drug, gap_ind, time2gap, gap_dt)
    
    # Return data
    return(dat.per)
}

## 1.4. Create function to combine all data ----
fun.com <- function(yr, gp){
    # Run adherence code for SGLT2i
    dat.adh.sglt2i <- fun.adh(cohort.base, "sglt2i", sglt2i, yr = yr)
    
    # Run adherence code for GLP1-RA
    dat.adh.glp1 <- fun.adh(cohort.base, "glp1", glp1, yr = yr)
    
    # Run adherence code for DPP4i
    dat.adh.dpp4i <- fun.adh(cohort.base, "dpp4i", dpp4i, yr = yr)
    
    # Run persistence code for SGLT2i
    dat.per.sglt2i <- fun.per(cohort.base, "sglt2i", sglt2i, yr = yr, gp = gp)
    
    # Run persistence code for GLP1-RA
    dat.per.glp1 <- fun.per(cohort.base, "glp1", glp1, yr = yr, gp = gp)
    
    # Run persistence code for DPP4i
    dat.per.dpp4i <- fun.per(cohort.base, "dpp4i", dpp4i, yr = yr, gp = gp)
    
    # Combine adherence data
    dat.adh <- 
        # Rowbind data
        rbind(
            # SGLT2i
            dat.adh.sglt2i,
            # GLP1-RA
            dat.adh.glp1,
            # DPP4i
            dat.adh.dpp4i)
    
    # Combine persistence data
    dat.per <- 
        # Rowbind data
        rbind(
            # SGLT2i
            dat.per.sglt2i,
            # GLP1-RA
            dat.per.glp1,
            # DPP4i
            dat.per.dpp4i)
    
    # Combine adherence and persistence data
    dat.com <- dat.adh %>%
        # Join persistence data
        left_join(dat.per, c("lopnr", "drug"))
    
    # Return data
    return(dat.com)
}

## 1.5. Create function for adherence in 3 month time periods ----
fun.adh.period <- function(df_base, drg, df_drug){
    # To get PDC per period:
    # Start a repeat loop in which we take the drugs up until the end of the period
    # Check if the period has people in it (if everyone is censored stop repeat loop)
    # Then calculate missed days up until the end of the period
    # Then calculate only the days missed in the period of interest
    # Calculate days in that period (could be < 91 if someone is censored)
    # Calculate PDC for that period
    
    # Apply corrections in data
    dat.drg <- fun.crc(df_drug)
    
    # Start lapply
    dat.fin <- do.call("rbind", lapply(0:29, \(x){
        # Prepare temporary dataset to check how many individuals are left in the period
        n <- df_base %>%
            # Keep only drug of interest
            filter(drug == drg) %>%
            # Arrange for grouping
            arrange(lopnr) %>%
            # Group on person (necessary to calculate correct dt_start)
            group_by(lopnr) %>%
            # Create new variables
            mutate(# Start date of period is equal to date index if first period or on 91 * i + i days after index date (this way a period is always 91 days)
                dt_start = as.Date(ifelse(x == 0, dt_index, dt_index + 91 * x + x), origin = "1970-01-01"),
                # End date is 91 days after start date
                dt_end = dt_start + 91,
                # If end date is after censoring, reset end date to censor date
                dt_end = as.Date(ifelse(dt_censor < dt_end, dt_censor, dt_end), origin = "1970-01-01")) %>%
            # Ungroup again
            ungroup() %>%
            # Remove individuals where the start date is on or after the censor date
            filter(!dt_start >= dt_censor) %>%
            # Keep only distinct identifiers
            distinct(lopnr) %>%
            # Get length
            nrow() %>%
            # As numeric
            as.numeric()

        if(n > 0){
            # Prepare drug dataset
            dat.drg.i <- df_base %>%
                # Keep only drug of interest
                filter(drug == drg) %>% 
                # Arrange for grouping
                arrange(lopnr) %>%
                # Group on person (necessary to calculate correct dt_start)
                group_by(lopnr) %>%
                # Create new variables
                mutate(# Start date of period is equal to date index if first period or on 91 * i + i days after index date (this way a period is always 91 days)
                       dt_start = as.Date(ifelse(x == 0, dt_index, dt_index + 91 * x + x), origin = "1970-01-01"),
                       # End date is 91 days after start date
                       dt_end = dt_start + 91,
                       # If end date is after censoring, reset end date to censor date
                       dt_end = as.Date(ifelse(dt_censor < dt_end, dt_censor, dt_end), origin = "1970-01-01")) %>%
                # Ungroup again
                ungroup() %>%
                # Remove individuals where the start date is on or after the censor date
                filter(!dt_start >= dt_censor) %>%
                # Join dispensation data
                left_join(dat.drg, "lopnr") %>%
                # Keep only drugs prescribed before or on the last day of observation 
                filter(edatum <= dt_end) %>% 
                # Remove unnecessary variables
                dplyr::select(-(drug:dt_end)) %>% 
                # Change edatum to date
                mutate(edatum = as.Date(edatum, origin = "1970-01-01")) %>%
                # Sort on person and date of prescription
                arrange(lopnr, edatum)
            
            # Prepare cohort data
            dat.adh.i <- df_base %>%
                # Keep only drug of interest
                filter(drug == drg) %>%
                # Arrange for grouping
                arrange(lopnr) %>%
                # Group on person (necessary to calculate correct dt_start)
                group_by(lopnr) %>%
                # Create new variables
                mutate(# Start date of period is equal to date index if first period or on 91 * i + i days after index date (this way a period is always 91 days)
                    dt_start = as.Date(ifelse(x == 0, dt_index, dt_index + 91 * x + x), origin = "1970-01-01"),
                    # End date is 91 days after start date
                    dt_end = dt_start + 91,
                    # If end date is after censoring, reset end date to censor date
                    dt_end = as.Date(ifelse(dt_censor < dt_end, dt_censor, dt_end), origin = "1970-01-01")) %>%
                # Ungroup again
                ungroup() %>%
                # Join prepared drugs data
                left_join(dat.drg.i, "lopnr") %>%
                # Create base variables
                mutate(# Amount of days the dispensation lasts
                       days = ifelse(drug == "glp1", (antal * dur), (antal * antnum)),
                       # End of dispensation
                       dt_disp_end = edatum + days,
                       # Indicator for days without pills
                       missed_days = 0,
                       # Date that pills run out
                       dt_out_of_pills = dt_disp_end) %>%
                # Arrange for grouping
                arrange(lopnr) %>%
                # Group per person
                group_by(lopnr) %>%
                # Calcualte new variables
                mutate(# Create row number per dispensation
                       rownr = row_number(),
                       # Calculate max amount of dispensations each person has
                       nrows = max(rownr),
                       # Keep only dt_out_of_pills for first dispensation, as this is the only one we know now
                       dt_out_of_pills = as.Date(ifelse(rownr == 1, dt_out_of_pills, NA), origin = "1970-01-01")) %>%
                # Ungroup again
                ungroup()
   
            # Iterate through dispensatiosn taking into account stockpiling
            for(i in 2:max(dat.adh.i[["nrows"]])){
                # Loop on each dispensation
                dat.adh.i <- dat.adh.i %>%
                    # Sort for grouping
                    arrange(lopnr) %>%
                    # Group by person
                    group_by(lopnr) %>%
                    # Create and update variables
                    mutate(# Placeholder out of pills date is previous one
                           dt_out_of_pills_lag = lag(dt_out_of_pills),
                           # If new edatum is later than out of pills date of previous observation, calcualte missed days, otherwise 0
                           missed_days = ifelse(rownr == i & edatum > dt_out_of_pills_lag, as.numeric(edatum - dt_out_of_pills_lag),
                                                ifelse(rownr == i & edatum <= dt_out_of_pills_lag, 0, missed_days)),
                           # If no days were missed, out of pills date is previous out of pills date + days, if days missed, out of pills date start over from current dispensation
                           dt_out_of_pills = as.Date(ifelse(rownr == i & missed_days == 0, dt_out_of_pills_lag + days,
                                                            ifelse(rownr == i & missed_days != 0, edatum + days, dt_out_of_pills)), origin = "1970-01-01")) %>%
                    # Ungroup again
                    ungroup()
            }
                
            # Make final calculations for adherence
            dat.adh.i <- dat.adh.i %>%
                # Keep only dispensations within relevant time frame
                filter(edatum >= dt_start & edatum <= dt_end) %>%
                # Arrange for grouping
                arrange(lopnr) %>%
                # Group per person
                group_by(lopnr) %>%
                # Create new variables
                mutate(# New rownumber
                       rownr = row_number(),
                       # If it is the first period, dt_out_of_pills_lag is NA, set to index date
                       dt_out_of_pills_lag = as.Date(ifelse(x == 0 & rownr == 1, dt_index, dt_out_of_pills_lag), origin = "1970-01-01"),
                       # In the first row of every person check the distance between last out of pills date and edatum to accurately calculate missed days from start of period
                       missed_days = case_when(rownr == 1 & edatum > dt_out_of_pills_lag & dt_out_of_pills_lag < dt_start ~ as.numeric(edatum - dt_start),
                                               rownr == 1 & edatum > dt_out_of_pills_lag & dt_out_of_pills_lag >= dt_start ~ as.numeric(edatum - dt_out_of_pills_lag),
                                               rownr == 1 & edatum <= dt_out_of_pills_lag ~ 0,
                                               rownr != 1 ~ missed_days),
                       # Limit out of pills date to end of observational period if larger
                       dt_out_of_pills = as.Date(ifelse(dt_out_of_pills > dt_end, dt_end, dt_out_of_pills), origin = "1970-01-01"),
                       # Calculate total of missed days
                       md = sum(missed_days)) %>%
                # Ungroup again
                ungroup() %>%
                # Rearrange, putting last dispensation at the top
                arrange(lopnr, desc(rownr)) %>%
                # Group on person
                group_by(lopnr) %>%
                # Keep only first row per person (i.e., only last dispensation)
                slice(1L) %>%
                # Create and update variables
                mutate(# Add to missings days any remaining days between last dispensation's out of pills date and end of period
                       md = round(md + as.numeric(dt_end - dt_out_of_pills)),
                       # Calculate total days in period per person
                       totdays = as.numeric(dt_end - dt_start),
                       # Days covered
                       dc = totdays - md,
                       # Calculate PDC
                       pdc = round(dc / totdays * 100),
                       # Create PDC indicator for low PDC
                       pdc_ind = ifelse(pdc >= 80, 0, 1)) %>%
                # Ungroup again
                ungroup()
            
            # Get group summary
            dat.fin.i <- 
                # Create data frame with basic information
                data.frame(period = x + 1,
                           n = n,
                           mean_pdc = mean(dat.adh.i[["pdc"]]),
                           sd = sd(dat.adh.i[["pdc"]]),
                           min_pdc = min(dat.adh.i[["pdc"]]),
                           max_pdc = max(dat.adh.i[["pdc"]]),
                           prop_low_pdc = sum(dat.adh.i[["pdc_ind"]]) / n * 100) %>%
                # Calculate further variables
                mutate(# Standard error
                       se = sd / sqrt(n),
                       # Lower limit of mean
                       ll = mean_pdc - 1.96 * se,
                       # Upper limit of mean
                       ul = mean_pdc + 1.96 * se)
            
            # Return data
            return(dat.fin.i)
        }
    }))
}
  
## 1.6. Code for adding last variables and renaming
fun.fin <- function(df, name_addon){
    # Join data
    dat.fin <- cohort %>%
        # Keep only relevant variables
        dplyr::select(lopnr, drug, dt_index, dt_censor, censor_reason) %>%
        # Join adherence and persistence data
        left_join(df, c("lopnr", "drug"))
    
    # Determine restart; this data frame contains indvidiuals who had a gap, but who had a new dispensation within 3 months
    dat.res <- dat.fin %>%
        # Keep only relevant variables
        dplyr::select(lopnr, drug, gap_dt, gap_ind) %>%
        # Join drug data
        left_join(
            # Combine drug data
            rbind(
                # SGLT2i with drug indicator
                mutate(fun.crc(sglt2i), drug = "sglt2i"),
                # GLP1-RA with drug indicator
                mutate(fun.crc(glp1), drug = "glp1") %>%
                    # Remove duration variable
                    dplyr::select(-dur),
                # DPP4i with drug indicator
                mutate(fun.crc(dpp4i), drug = "dpp4i")
            ), c("lopnr", "drug")) %>%
        # Keep only dispensations after date of gap but not later than 3 months after
        filter(as.Date(edatum) > gap_dt & as.Date(edatum) <= (gap_dt + 90)) %>%
        # Keep only individuals who had a gap
        filter(gap_ind == 1) %>%
        # Arrange on lopnr, drug and earliest drug date
        arrange(lopnr, drug, as.Date(edatum)) %>%
        # Group by lopnr and drug
        group_by(lopnr, drug) %>%
        # Keep first row per person and drug
        slice(1L) %>%
        # Ungroup again
        ungroup() %>%
        # Calculate time to restart and add indicator for restart
        mutate(time2restart = as.numeric(as.Date(edatum) - gap_dt),   # Time to restart
               res_ind = 1) %>%                                       # Restart indicator
        # Keep only relevant variables
        dplyr::select(lopnr, drug, res_ind, time2restart)
    
    # Finalize data
    dat.fin <- dat.fin %>%
        # Join restart data
        left_join(dat.res, c("lopnr", "drug")) %>%
        # If there was a gap but no restart, set indicator to 0 and calculate time2restart as time to censoring
        mutate(res_ind = ifelse(gap_ind == 1 & is.na(res_ind), 0, res_ind),                            # Setting indicator to 0
               time2restart = ifelse(res_ind == 0, as.numeric(dt_censor - gap_dt), time2restart)) %>%  # Calculating time to restart for non-restarters
        # Create variable for cumulative incidence function (i.e., multi-level indicator for event with mortality as censoring)
        mutate(cif_ind = ifelse(censor_reason == "death" & dt_censor <= gap_dt, 2,               # Indicator for event or censoring
                                ifelse(gap_ind == 1, 1, 0)),
               time2cif = case_when(cif_ind == 0 ~ as.numeric(dt_censor - dt_index),             # Time to event or censoring
                                    cif_ind == 1 ~ as.numeric(gap_dt - dt_index),
                                    cif_ind == 2 ~ as.numeric(dt_censor - dt_index))) %>%
        # Keep only relevant variables
        dplyr::select(lopnr, drug, pdc, pdc_ind, gap_ind, gap_dt, time2gap, res_ind, time2restart, cif_ind, time2cif)
    
    # Define column names
    colnames(dat.fin) <- c("lopnr", "drug", paste0(colnames(dat.fin)[3:11], "_", name_addon))
    
    # Return data
    return(dat.fin)
}
                   
# 2. Run codes ----
## 2.0. Set-up ----
# Load data. We use raw  data for GLP1-RA to be able to link duration to brand name
load("sglt2i.Rdata"); load("dpp4i.Rdata"); load("glp1_raw.Rdata"); load("glp1.Rdata")

# Load duration for GLP1-RA data as this is injection
glp1_durations <- read_xlsx("C:/Users/rjjan/OneDrive/Bureaublad/8. GDT_SGLT2i/GLP1-RA.xlsx")

# Add duration to GLP1-RA raw data
glp1_interim <- glp1_raw %>%
    # Keep only relevant variables
    dplyr::select(LopNr, ATC, EDATUM, OTYP, ANTAL, forpddd, antnum, AR, produkt) %>%
    # Create dummy variables that can be removed later (necessary to make fun.crc() work
    mutate(pop = NA, edatum_c = NA) %>%
    # Rename variables
    rename(lopnr = 1, atc = 2, edatum = 3, otyp = 4, antal = 5, ar = 8) %>%
    # Join durations
    left_join(dplyr::select(glp1_durations, `ATC code`, `Brand name`, `Package duration (days)`) %>% 
                      # Rename columns
                      rename(atc = 1, produkt = 2, dur = 3), 
                  c("atc", "produkt")) %>%
    # Keep only relevant variables and set order at the same time
    dplyr::select(lopnr, atc, edatum, dur)
    
# Add duration to GLP1-RA data 
glp1 <- glp1 %>%
    # Join data
    left_join(glp1_interim, c("lopnr", "atc", "edatum"))

## 2.1. 1 year & 60 day grace period ----
# Calculate adherence and persistence
dat.y1gp60 <- fun.com(yr = 1, gp = 60)

# Save data
save(dat.y1gp60, file = "dat_y1gp60.Rdata")

## 2.2. 2 year & 60 day grace period ----
# Calculate adherence and persistence
dat.y2gp60 <- fun.com(yr = 2, gp = 60)

# Save data
save(dat.y2gp60, file = "dat_y2gp60.Rdata")

## 2.3. No limit & 60 day grace period ----
# Calculate adherence and persistence
dat.yigp60 <- fun.com(yr = 0, gp = 60)

# Save data
save(dat.yigp60, file = "dat_yigp60.Rdata")

## 2.4. 1 year & 30 day grace period ----
# Calculate adherence and persistence
dat.y1gp30 <- fun.com(yr = 1, gp = 30)

# Save data
save(dat.y1gp30, file = "dat_y1gp30.Rdata")

# 3. Finalize data ----
# Load covariates
load("covariates.Rdata")

# Load adherence and persistence data (optional)
load("dat_y1gp60.Rdata"); load("dat_y2gp60.Rdata"); load("dat_yigp60.Rdata"); load("dat_y1gp30.Rdata")

# Combine everything
pop <- cohort %>%
    # Join covariates
    left_join(covariates, c("lopnr", "drug")) %>%
    # Join adherence and persistence measures after finalizing
    left_join(fun.fin(dat.y1gp60, name_addon = "y1gp60"), c("lopnr", "drug")) %>%       # 1 year and 60 day grace period
    left_join(fun.fin(dat.y2gp60, name_addon = "y2gp60"), c("lopnr", "drug")) %>%       # 2 years and 60 day grace period
    left_join(fun.fin(dat.yigp60, name_addon = "yigp60"), c("lopnr", "drug")) %>%       # No year limit and 60 day grace period
    left_join(fun.fin(dat.y1gp30, name_addon = "y1gp30"), c("lopnr", "drug"))           # 1 year and 30 day grace period
    
# Save data
save(pop, file = "pop.Rdata")
start <- Sys.time()
# 4. Create separate data for adherence per period ----
## 4.1. Run code for SGLT2i ----
# Run code
dat.sglt2i.adh.period <- fun.adh.period(cohort.base, "sglt2i", sglt2i)%>%
    # Add drug column
    mutate(drug = "sglt2i")

# Save data
save(dat.sglt2i.adh.period, file = "ap_sglt2i.Rdata")

## 4.2. Run code for GLP1-RA ----
# Run code
dat.glp1.adh.period <- fun.adh.period(cohort.base, "glp1", glp1) %>%
    # Add drug column
    mutate(drug = "glp1")

# Save data
save(dat.glp1.adh.period, file = "ap_glp1.Rdata")

## 4.3. Run code for DPP4i ----
# Run code
dat.dpp4i.adh.period <- fun.adh.period(cohort.base, "dpp4i", dpp4i)%>%
    # Add drug column
    mutate(drug = "dpp4i")

# Save data
save(dat.dpp4i.adh.period, file = "ap_dpp4i.Rdata")

## 4.4. Combine all data ----
# Combine data
dat.adh.period <- 
    # Bind data
    rbind(# SGLT2i
          dat.sglt2i.adh.period,
          # GLP1-RA
          dat.glp1.adh.period,
          # DPP4i
          dat.dpp4i.adh.period)

# Save data
save(dat.adh.period, file = "ap.Rdata")
finish <- Sys.time()
finish-start