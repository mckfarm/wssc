## ---------------------------
## Script name: glycogen.R
## Author: McKenna Farmer
## Date Created: 2022-12-12
## ---------------------------
## Notes:
##   glycogen parsing and data analysis
##
## ---------------------------

# packages -------
library(tidyverse)
library(readxl)
library(lubridate)
library(MetBrewer)

source(file.path("scripts","plotting.R"))

# data read in and parsing --------
in_path <- "C:/Users/mckyf/Northwestern University/Wells Research Group - WRF 5071 (WSSC resource efficient nutrient removal)/4 Data and results/Glycogen quantification"

df <- read_excel(path=file.path(in_path,"all_results.xlsx"), sheet="data")
op <- read_excel(path=file.path(in_path,"all_results.xlsx"), sheet="op_log")

# df parsing - glycogen results
df$date <- as_date(df$date)

df$battery <- recode(df$battery, "AT1" = "test", "AT2" = "control")
df$battery <- factor(df$battery, levels=c("test","control","RAS"))


location_list <- c("AN1", "AX3", "OX4", "SZ1", "SZ2", "RAS")
df$location <- factor(df$location, levels=location_list)

df$mass1 <- df$conc1 * df$volume / df$dry_mass
df$mass2 <- df$conc2 * df$volume / df$dry_mass

df <- df %>% rowwise() %>%
  mutate(mass_ave=mean(c(mass1, mass2)),
         mass_sd=sd(c(mass1,mass2)))

exclude_gly_dates <- c(ymd("2021-09-15"), ymd("2022-06-13"))

all_gly_dates <- unique(df$date)

gly_dates <- setdiff(all_gly_dates, exclude_gly_dates)

saveRDS(gly_dates, file=file.path("data","gly_dates.rds"))


# op log parsing

op$date <- as_date(op$date)
op$battery <- "test"

df_select <- df %>% 
  filter(!date %in% exclude_gly_dates) %>%
  left_join(., op, by=c("date","battery"))

# testing correlations
df_select %>% filter(battery=="test") %>%
ggplot(., aes(x=DO_SZ1, y=mass_ave)) +
  geom_point() +
  stat_cor(method="pearson") +
  theme_classic()

# summary stats
df_select %>% filter(battery != "RAS") %>% group_by(battery, date) %>%
  summarise(diff = max(mass_ave) - min(mass_ave)) %>%
  pivot_wider(names_from=battery, values_from=diff)

df_select %>% filter(battery != "RAS") %>% group_by(battery, date) %>%
  summarise(ave=mean(mass_ave),
            sd=sd(mass_ave))

df_select %>% filter(battery == "RAS") %>%
  select(c(date, mass_ave, mass_sd))

# plotting
# control and test
ggplot(df_select, aes(x=location, y=mass_ave, 
                     color=battery, group=battery, shape=battery)) +
  facet_wrap(~date, nrow=2) + 
  geom_point(alpha=0.5) + 
  geom_line() + 
  geom_errorbar(aes(ymin=mass_ave-mass_sd, ymax=mass_ave+mass_sd),
                alpha=0.5) +
  labs(x="Sampling location", y="Glycogen \n[mg glucose/mg dry biomass]") +
  theme_bw() +
  color_battery_gly +
  shape_battery_gly +
  theme(legend.position="top")
ggsave(filename=file.path("results","glycogen_all.png"),
       height=4.5, width=8, units="in", dpi=600)


# test only
df_select2 <- df_select %>% filter(battery=="test")

ggplot(df_select2, aes(x=location, y=mass_ave, group=date)) +
  facet_wrap(~date, nrow=2) + 
  geom_point(alpha=0.5) + 
  geom_line() + 
  geom_errorbar(aes(ymin=mass_ave-mass_sd, ymax=mass_ave+mass_sd),
                alpha=0.5) +
  labs(x="Sampling location", y="Glycogen \n[mg glucose/mg dry biomass]") +
  theme_bw()
ggsave(filename=file.path("results","glycogen_test_all.png"),
       height=4.5, width=6.5, units="in", dpi=600)






# relationship to performance variables ----
lims <- readRDS(file.path("data","lims_data.rds"))
fit <- readRDS(file.path("data","fit_data.rds"))
gly_dates <- sort(gly_dates)

# add and subtract a weeks worth of dates

gly_dates_add <- gly_dates
gly_dates_ref <- gly_dates
for (x in 1:20) {
  gly_dates_add <- append(gly_dates_add, gly_dates - x)
  gly_dates_ref <- append(gly_dates_ref, gly_dates)
}

df_gly_dates <- tibble(date=as_date(gly_dates_add), 
                       ref_date=factor(as_date(gly_dates_ref)))


inf_eff <- lims %>% filter(`LIMS Sample` %in% c("Raw Influent + Recycle", "Secondary Effluent 1")) %>%
  select(c("Collect Date", "Parameter Abbrev", "LIMS Sample", "Result"))

colnames(inf_eff) <- c("date", "parameter", "location", "value")
inf_eff$date <- as_date(inf_eff$date)
inf_eff <- inf_eff[!duplicated(inf_eff), ] # remove duplicate rows

# adding glycogen reference dates
inf_eff <- left_join(inf_eff, df_gly_dates, by="date")

inf_eff_ratio <- inf_eff %>% filter(location=="Raw Influent + Recycle") %>%
  filter(parameter %in% c("BOD","TP")) %>%
  pivot_wider(id_cols=c(date, ref_date), names_from=parameter, values_from=value)

inf_eff_ratio$BOD_TP <- inf_eff_ratio$BOD / inf_eff_ratio$TP


y_gly <- 15
inf_eff %>% filter(date %in% gly_dates_add) %>% 
  filter(parameter=="BOD") %>%
  ggplot(., aes(x=date, y=value, color=ref_date)) +
  geom_point() + 
  annotate(geom = "text", x = as_date(gly_dates[1]),
           y = y_gly, label = "1") +
  annotate(geom = "text", x = as_date(gly_dates[2]),
           y = y_gly, label = "2") +
  annotate(geom = "text", x = as_date(gly_dates[3]),
           y = y_gly, label = "3") +
  annotate(geom = "text", x = as_date(gly_dates[4]),
           y = y_gly, label = "4") +
  annotate(geom = "text", x = as_date(gly_dates[5]),
           y = y_gly, label = "5") +
  annotate(geom = "text", x = as_date(gly_dates[6]),
           y = y_gly, label = "6") +
  labs(y="BOD [mg/L]", x="Date") +
  theme_bw()


y_gly <- 9
inf_eff %>% filter(date %in% gly_dates_add) %>% 
  filter(parameter=="TP") %>%
  ggplot(., aes(x=date, y=value, color=ref_date)) +
  geom_point() + 
  annotate(geom = "text", x = as_date(gly_dates[1]),
           y = y_gly, label = "1") +
  annotate(geom = "text", x = as_date(gly_dates[2]),
           y = y_gly, label = "2") +
  annotate(geom = "text", x = as_date(gly_dates[3]),
           y = y_gly, label = "3") +
  annotate(geom = "text", x = as_date(gly_dates[4]),
           y = y_gly, label = "4") +
  annotate(geom = "text", x = as_date(gly_dates[5]),
           y = y_gly, label = "5") +
  annotate(geom = "text", x = as_date(gly_dates[6]),
           y = y_gly, label = "6") +
  labs(y="Total P [mgP/L]", x="Date") +
  theme_bw()



ggplot(inf_eff_ratio %>% filter(date %in% gly_dates_add), 
       aes(x=date, y=BOD_TP, color=ref_date)) +
  geom_boxplot() + 
  geom_point() + 
  annotate(geom = "text", x = as_date(gly_dates[1]),
           y = y_gly, label = "1") +
  annotate(geom = "text", x = as_date(gly_dates[2]),
           y = y_gly, label = "2") +
  annotate(geom = "text", x = as_date(gly_dates[3]),
           y = y_gly, label = "3") +
  annotate(geom = "text", x = as_date(gly_dates[4]),
           y = y_gly, label = "4") +
  annotate(geom = "text", x = as_date(gly_dates[5]),
           y = y_gly, label = "5") +
  annotate(geom = "text", x = as_date(gly_dates[6]),
           y = y_gly, label = "6") +
  labs(y="BOD:TP ratio", x="Date") +
  theme_bw() +
  theme(legend.position="none")

ggsave(filename=file.path("results","bod_tp_gly.png"), 
       height=4, width=5, units="in", dpi=600)


