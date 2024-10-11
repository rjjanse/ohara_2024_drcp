##############################################################################################################################################
# Functions to determine adherence and persistence
# Roemer Janse - 27/01/2022
# -------------------------------------------------------------------------------------------------------------------------------------------#


##### Function for adherence #####
adh <- function(drg, drug_df, yr){
    # Prepare drugs data
    drugs_df <- filter(sample, drug == drg) %>% mutate(yr = yr,
                                                       year_dt = as.Date(ifelse(yr == 0, censor_dt_cod, (index_dt + (yr * 365.25))),
                                                                         origin = "1970-01-01"),
                                                       year_dt = as.Date(ifelse(year_dt > censor_dt_cod, censor_dt_cod, year_dt),
                                                                         origin = "1970-01-01")) %>% left_join(drug_df, "lopnr") %>%
        filter(edatum <= year_dt) %>% dplyr::select(-(drug:year_dt)) %>% arrange(lopnr, edatum)
    
    # Prepare data for the loop
    adh_drug <- filter(sample, drug == drg) %>% mutate(yr = yr,
                                                       year_dt = as.Date(ifelse(yr == 0, censor_dt_cod, (index_dt + (yr * 365.25))),
                                                                         origin = "1970-01-01"),
                                                       year_dt = as.Date(ifelse(year_dt > censor_dt_cod, censor_dt_cod, year_dt),
                                                                         origin = "1970-01-01"),
                                                       totday = as.numeric(year_dt - index_dt)) %>% left_join(drugs_df, "lopnr") %>%
        mutate(days = ifelse(drug == "glp1", (antal * dur), (antal * antnum)), 
               end_dt = as.Date(edatum + days),
               # Days without pills
               missed_days = 0,
               # Out of pills date
               out_of_pills_dt = end_dt) %>% 
        arrange(lopnr) %>% group_by(lopnr) %>% mutate(rownr = row_number(),
                                                      nrow = max(rownr),
                                                      out_of_pills_dt = as.Date(ifelse(rownr == 1, out_of_pills_dt, NA), 
                                                                                origin = "1970-01-01")) %>% ungroup()
    
    # Loop to also take into account stockpiling
    for(i in 2:max(adh_drug$rownr)){
        adh_drug <- adh_drug %>% arrange(lopnr) %>% group_by(lopnr) %>% 
            mutate(lag_out_of_pills_dt = lag(out_of_pills_dt),
                   # If new edatum is later than out of pills date (previous observation), calculate missed days, otherwise 0
                   missed_days = ifelse(rownr == i & edatum > lag_out_of_pills_dt, as.numeric(edatum - lag_out_of_pills_dt), 
                                        ifelse(rownr == i & edatum <= lag_out_of_pills_dt, 0, missed_days)),
                   # If no days missed, out of pills date is lag out of pills date + days
                   # If days were missed, out of pills date starts anew from edatum
                   out_of_pills_dt = as.Date(ifelse(rownr == i & missed_days == 0, lag_out_of_pills_dt + days,
                                                    ifelse(rownr == i & missed_days != 0, edatum + days, out_of_pills_dt)), 
                                             origin = "1970-01-01")) %>%
            ungroup()
    }
    
    # Determine adherence                                                  
    adh_drug <- adh_drug %>% arrange(lopnr, desc(rownr)) %>% group_by(lopnr) %>% 
        mutate(out_of_pills_dt = as.Date(ifelse(out_of_pills_dt > year_dt, year_dt, out_of_pills_dt), origin = "1970-01-01"),
               md = sum(missed_days),
               md = md + as.numeric(year_dt - out_of_pills_dt)) %>% 
        slice(1L) %>% ungroup() %>%
        mutate(dc = totday - md,
               pdc = round(dc / totday * 100),
               pdc_ind = ifelse(pdc >= 80, 0, 1)) %>% dplyr::select(lopnr, drug, pdc, pdc_ind)
    
    return(adh_drug)
}


##############################################################################################################################################


##### Function for persistence #####
per <- function(drg, drug_df, yr, gp){
    # To get persistence (60 day gap)
    # First determine when a prescription ends
    # Then calculate the gap between each end date and subsequent start date
    drugs_df <- filter(sample, drug == drg) %>% mutate(yr = yr,
                                                       year_dt = as.Date(ifelse(yr == 0, censor_dt_cod, (index_dt + (yr * 365.25))),
                                                                         origin = "1970-01-01"),
                                                       year_dt = as.Date(ifelse(year_dt > censor_dt_cod, censor_dt_cod, year_dt),
                                                                         origin = "1970-01-01")) %>% left_join(drug_df, "lopnr") %>%
        filter(edatum <= year_dt) %>% dplyr::select(-(drug:year_dt)) %>% arrange(lopnr, edatum)
    
    per_drug <- filter(sample, drug == drg) %>% mutate(yr = yr,
                                                       year_dt = as.Date(ifelse(yr == 0, censor_dt_cod, (index_dt + (yr * 365.25))),
                                                                         origin = "1970-01-01"),
                                                       year_dt = as.Date(ifelse(year_dt > censor_dt_cod, censor_dt_cod, year_dt),
                                                                         origin = "1970-01-01")) %>% left_join(drugs_df, "lopnr") %>%
        mutate(days = ifelse(drug == "glp1", (antal * dur), (antal * antnum)),
               end_dt = as.Date(edatum + days)) %>%
        arrange(lopnr, edatum) %>% group_by(lopnr) %>% mutate(rownr = row_number(),
                                                              nrow = max(rownr),
                                                              out_of_pills_dt = as.Date(ifelse(rownr == 1, end_dt, NA), 
                                                                                        origin = "1970-01-01"),
                                                              gap = NA) %>% ungroup()
    
    # Take into account stockpiling
    for(i in 1:max(per_drug$rownr)){
        if(i == 1){
            per_drug <- per_drug %>% arrange(lopnr) %>% group_by(lopnr) %>%
                mutate(lead_edatum = lead(edatum),
                       lead_edatum = as.Date(ifelse(is.na(lead_edatum), year_dt, lead_edatum), origin = "1970-01-01"),
                       # If new edatum is later than out of pills date, calculate gap
                       gap = ifelse(rownr == 1 & lead_edatum > out_of_pills_dt, as.numeric(lead_edatum - out_of_pills_dt), 
                                    ifelse(rownr == 1 & lead_edatum <= out_of_pills_dt, 0, gap))) %>% ungroup()
        }
        
        else{
            per_drug <- per_drug %>% arrange(lopnr) %>% group_by(lopnr) %>%
                mutate(lag_gap = lag(gap),
                       lag_out_of_pills_dt = lag(out_of_pills_dt),
                       # If there was a gap, no matter how long, out of pills date is reset to end date of current edatum, otherwise add
                       # days of current dispensatio nto out of pills date
                       out_of_pills_dt = as.Date(ifelse(rownr == i & lag_gap == 0, lag_out_of_pills_dt + days, 
                                                        ifelse(rownr == i & lag_gap != 0, end_dt, out_of_pills_dt)), origin = "1970-01-01"),
                       lead_edatum = lead(edatum),
                       lead_edatum = as.Date(ifelse(is.na(lead_edatum), year_dt, lead_edatum), origin = "1970-01-01"),
                       # If new edatum is later than out of pills date, calculate gap, otherwise 0
                       gap = ifelse(rownr == i & lead_edatum > out_of_pills_dt, as.numeric(lead_edatum - out_of_pills_dt), 
                                    ifelse(rownr == i & lead_edatum <= out_of_pills_dt, 0, gap))) %>% ungroup()
        }
    }
    
    per_drug <- per_drug %>%
        mutate(gap_ind = ifelse(gap >= gp, 1, 0),
               # Gap date for people with a gap is end date plus grace period, otherwise end date
               gap_dt = as.Date(ifelse(gap_ind == 1, (out_of_pills_dt + gp), year_dt), origin = "1970-01-01"),
               # Gap indicators for people of whom the grace period ends after the year date are 0 (no gap)
               gap_ind = ifelse(gap_dt > year_dt & gap_ind == 1, 0, gap_ind),
               # End date for end dates past year date is year date
               # This also sets end dates with grace period after one year for people who went from gap (1) to no gap (0) back to right date
               gap_dt = as.Date(ifelse(gap_dt > year_dt, year_dt, gap_dt), origin = "1970-01-01")) %>%
        arrange(lopnr, desc(gap_ind), gap_dt) %>% group_by(lopnr) %>% slice(1L) %>% ungroup() %>%
        mutate(time2gap = ifelse(gap_ind == 1, as.numeric(gap_dt - index_dt), as.numeric(year_dt - index_dt))) %>%
        dplyr::select(lopnr, drug, gap_ind, time2gap, gap_dt)
    
    return(per_drug)
}