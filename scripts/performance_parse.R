## ---------------------------
## Script name: performance_parse.R
## Author: McKenna Farmer
## Date Created: 2022-12-06
## ---------------------------
## Notes:
##   Parsing performance data from Excel file into R data file
##
## ---------------------------

# packages -------
library(readr)
library(dplyr)
library(readxl)


# LIMS
lims <- read_excel(file.path("data","Seneca WRF5071 Data 221206.xlsx"), sheet="Raw LIMS", range=cell_cols("A:V"))

lims_to_keep <- c("Aeration Basin 1", "Aeration Basin 2", "Secondary Effluent 1", "Secondary Effluent 2", 
                 "Raw Influent", "Raw Influent + Recycle", "AT 1", "AT 2")

lims_select <- lims %>% filter(`LIMS Sample` %in% lims_to_keep)

saveRDS(lims_select, file.path("data", "lims_data.rds"))

# FIT 

fit <- read_excel(file.path("data","Seneca WRF5071 Data 221206.xlsx"), sheet="Raw FIT data", range=cell_cols("AS:BB"))

fit_to_keep <- c("AT 1", "AT 2", "FE", "Raw+Recycle", "RAS", "Raw")

fit_select <- fit %>% filter(LocationAbbrev %in% fit_to_keep)

saveRDS(fit_select, file.path("data", "fit_data.rds"))
