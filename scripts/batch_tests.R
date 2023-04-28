## ---------------------------
## Script name: batch_tests.R
## Author: McKenna Farmer
## Date Created: 2022-12-06
## ---------------------------
## Notes:
##   Figures and calcs for batch testing results
##
## ---------------------------

# packages -------
library(tidyverse)
library(readxl)
library(lubridate)

library(ggplot2)
library(cowplot)
library(MetBrewer)

# plotting params -------


list_nut <- c("op","ammonia","nitrate","nitrite")
list_cod <- c("sCOD","acetate","propionate")

color_nut <- scale_color_manual(
  limits = c("op","ammonia","nitrate","nitrite"),
  labels = c("OP","Ammonia","Nitrate","Nitrite"),
  values = met.brewer("Egypt", 4),
  name = "Parameter"
)

shape_nut <- scale_shape_manual(
  limits = c("op","ammonia","nitrate","nitrite"),
  labels = c("OP","Ammonia","Nitrate","Nitrite"),
  values = c(15,16,17,18),
  name = "Parameter"
)

color_cod <- scale_color_manual(
  limits = c("sCOD","acetate","propionate"), 
  labels = c("sCOD","Acetate","Propionate"), 
  values = met.brewer("Lakota", 3),
  name = "Parameter"
)

shape_cod <- scale_shape_manual(
  limits = c("sCOD","acetate","propionate"), 
  labels = c("sCOD","Acetate","Propionate"), 
  values = c(15,16,17),
  name = "Parameter"
)

phases <- data.frame(x1=90, x2=390, x3=460,
                     lab1=40, lab2=240, lab3=420, lab4=485)

phase_lines <- list(
  geom_vline(xintercept=phases$x1, color="grey", alpha=0.5),
  geom_vline(xintercept=phases$x2, color="grey", alpha=0.5),
  geom_vline(xintercept=phases$x3, color="grey", alpha=0.5))


# data read in ------
in_path <- "C:/Users/mckyf/Northwestern University/Wells Research Group - WRF 5071 (WSSC resource efficient nutrient removal)/4 Data and results/Ex-Situ Batch Assays"
low <- read_excel(file.path(in_path,"batch_compiled.xlsx"), sheet="low")
high <- read_excel(file.path(in_path,"batch_compiled.xlsx"), sheet="high")

low <- low %>% pivot_longer(!c(inf, date, do_sp, time_tot, time_phase),
                                 names_to="parameter", values_to="value") %>% drop_na()

high <- high %>% pivot_longer(!c(inf, date, do_sp, time_tot, time_phase),
                            names_to="parameter", values_to="value")%>% drop_na()

# nutrient plots ----

y_lim <- 80

low_syn <- low %>% filter(parameter %in% list_nut) %>% filter(inf=="syn") %>% 
  ggplot(., aes(x=time_tot, y=value, color=parameter, shape=parameter)) +
  phase_lines + 
  geom_point() +
  geom_line() +
  theme_classic() +
  theme(legend.position = "none") + 
  ylim(0,y_lim) +
  labs(x="", y="Concentration\n[mgN/L or mgP/L]", title="Synthetic feed, DO setpoint = 0.2 mg/L") +
  color_nut + 
  shape_nut +
  annotate(geom="text", x=phases$lab1, y=y_lim, label="Anaerobic", size=3) +
  annotate(geom="text", x=phases$lab2, y=y_lim, label="Aerobic", size=3) +
  annotate(geom="text", x=phases$lab3, y=y_lim, label="Anoxic", size=3) +
  annotate(geom="text", x=phases$lab4, y=y_lim, label="Aerobic", size=3)


y_lim <- 100
low_raw <- low %>% filter(parameter %in% list_nut) %>% filter(inf=="raw") %>% 
  ggplot(., aes(x=time_tot, y=value, color=parameter, shape=parameter)) +
  phase_lines + 
  geom_point() +
  geom_line() +
  theme_classic() +
  theme(legend.position = "bottom") + 
  ylim(0,y_lim) +
  labs(x="Time [min]", y="Concentration\n[mgN/L or mgP/L]", title="Raw influent, DO setpoint = 0.2 mg/L") +
  color_nut + 
  shape_nut +
  annotate(geom="text", x=phases$lab1, y=y_lim, label="Anaerobic", size=3) +
  annotate(geom="text", x=phases$lab2, y=y_lim, label="Aerobic", size=3) +
  annotate(geom="text", x=phases$lab3, y=y_lim, label="Anoxic", size=3) +
  annotate(geom="text", x=phases$lab4, y=y_lim, label="Aerobic", size=3)

legend_nut <- get_legend(low_raw)

low_raw <- low_raw + theme(legend.position="none")

plot_grid(legend_nut, low_syn, low_raw, 
          nrow=3, rel_heights=c(0.1,1,1), axis="tblr", align="hv")

ggsave(file.path("results","lowdo_batch.jpg"),
       height=7, width=6, dpi=600)




y_lim <- 90
high_syn <- high %>% filter(parameter %in% list_nut) %>% filter(inf=="syn") %>% 
  ggplot(., aes(x=time_tot, y=value, color=parameter, shape=parameter)) +
  phase_lines + 
  geom_point() +
  geom_line() +
  theme_classic() +
  theme(legend.position = "none") + 
  ylim(0,y_lim) +
  labs(x="", y="Concentration\n[mgN/L or mgP/L]", title="Synthetic feed, DO setpoint = 2 mg/L") +
  color_nut + 
  shape_nut +
  annotate(geom="text", x=phases$lab1, y=y_lim, label="Anaerobic", size=3) +
  annotate(geom="text", x=phases$lab2, y=y_lim, label="Aerobic", size=3) +
  annotate(geom="text", x=phases$lab3, y=y_lim, label="Anoxic", size=3) +
  annotate(geom="text", x=phases$lab4, y=y_lim, label="Aerobic", size=3)


y_lim <- 50
high_raw <- high %>% filter(parameter %in% list_nut) %>% filter(inf=="raw") %>% 
  ggplot(., aes(x=time_tot, y=value, color=parameter, shape=parameter)) +
  phase_lines + 
  geom_point() +
  geom_line() +
  theme_classic() +
  theme(legend.position = "none") + 
  ylim(0,y_lim) +
  labs(x="Time [min]", y="Concentration\n[mgN/L or mgP/L]", title="Raw influent, DO setpoint = 2 mg/L") +
  color_nut + 
  shape_nut +
  annotate(geom="text", x=phases$lab1, y=y_lim, label="Anaerobic", size=3) +
  annotate(geom="text", x=phases$lab2, y=y_lim, label="Aerobic", size=3) +
  annotate(geom="text", x=phases$lab3, y=y_lim, label="Anoxic", size=3) +
  annotate(geom="text", x=phases$lab4, y=y_lim, label="Aerobic", size=3)

plot_grid(legend_nut, high_syn, high_raw, 
          nrow=3, rel_heights=c(0.1,1,1), axis="tblr", align="hv")

ggsave(file.path("results","highdo_batch.jpg"),
       height=7, width=6, dpi=600)




# carbon profiles ------

y_lim <- 350
high_cod <- high %>% filter(parameter %in% list_cod) %>% filter(inf=="syn") %>% 
  ggplot(., aes(x=time_tot, y=value, color=parameter, shape=parameter)) +
  phase_lines + 
  geom_point() +
  geom_line() +
  theme_classic() +
  theme(legend.position = "none") + 
  ylim(0,y_lim) +
  shape_cod +
  color_cod + 
  labs(x="", y="Concentration [mgCOD/L]", title="Synthetic feed, DO setpoint = 2 mg/L") +
  annotate(geom="text", x=phases$lab1, y=y_lim, label="Anaerobic", size=3) +
  annotate(geom="text", x=phases$lab2, y=y_lim, label="Aerobic", size=3) +
  annotate(geom="text", x=phases$lab3, y=y_lim, label="Anoxic", size=3) +
  annotate(geom="text", x=phases$lab4, y=y_lim, label="Aerobic", size=3)

y_lim <- 400
low_cod <- low %>% filter(parameter %in% list_cod) %>% filter(inf=="syn") %>% 
  ggplot(., aes(x=time_tot, y=value, color=parameter, shape=parameter)) +
  phase_lines + 
  geom_point() +
  geom_line() +
  theme_classic() +
  theme(legend.position = "bottom") + 
  ylim(0,y_lim) +
  labs(x="Time [min]", y="Concentration [mgCOD/L]", title="Synthetic feed, DO setpoint = 0.2 mg/L") +
  shape_cod +
  color_cod + 
  annotate(geom="text", x=phases$lab1, y=y_lim, label="Anaerobic", size=3) +
  annotate(geom="text", x=phases$lab2, y=y_lim, label="Aerobic", size=3) +
  annotate(geom="text", x=phases$lab3, y=y_lim, label="Anoxic", size=3) +
  annotate(geom="text", x=phases$lab4, y=y_lim, label="Aerobic", size=3)

legend_cod <- get_legend(low_cod)
low_cod <- low_cod + theme(legend.position = "none")

plot_grid(legend_cod, high_cod, low_cod, 
          nrow=3, rel_heights=c(0.1,1,1), axis="tblr", align="hv")

ggsave(file.path("results","syn_cod_batch.jpg"),
       height=7, width=6, dpi=600)



y_lim <- 75
high_cod <- high %>% filter(parameter %in% list_cod) %>% filter(inf=="raw") %>% 
  ggplot(., aes(x=time_tot, y=value, color=parameter, shape=parameter)) +
  phase_lines + 
  geom_point() +
  geom_line() +
  theme_classic() +
  theme(legend.position = "none") + 
  ylim(0,y_lim) +
  shape_cod +
  color_cod + 
  labs(x="", y="Concentration [mgCOD/L]", title="Raw influent, DO setpoint = 2 mg/L") +
  annotate(geom="text", x=phases$lab1, y=y_lim, label="Anaerobic", size=3) +
  annotate(geom="text", x=phases$lab2, y=y_lim, label="Aerobic", size=3) +
  annotate(geom="text", x=phases$lab3, y=y_lim, label="Anoxic", size=3) +
  annotate(geom="text", x=phases$lab4, y=y_lim, label="Aerobic", size=3)

y_lim <- 75
low_cod <- low %>% filter(parameter %in% list_cod) %>% filter(inf=="raw") %>% 
  ggplot(., aes(x=time_tot, y=value, color=parameter, shape=parameter)) +
  phase_lines + 
  geom_point() +
  geom_line() +
  theme_classic() +
  theme(legend.position = "bottom") + 
  ylim(0,y_lim) +
  labs(x="Time [min]", y="Concentration [mgCOD/L]", title="Raw influent, DO setpoint = 0.2 mg/L") +
  shape_cod +
  color_cod + 
  annotate(geom="text", x=phases$lab1, y=y_lim, label="Anaerobic", size=3) +
  annotate(geom="text", x=phases$lab2, y=y_lim, label="Aerobic", size=3) +
  annotate(geom="text", x=phases$lab3, y=y_lim, label="Anoxic", size=3) +
  annotate(geom="text", x=phases$lab4, y=y_lim, label="Aerobic", size=3)

legend_cod <- get_legend(low_cod)
low_cod <- low_cod + theme(legend.position = "none")

plot_grid(legend_cod, high_cod, low_cod, 
          nrow=3, rel_heights=c(0.1,1,1), axis="tblr", align="hv")

ggsave(file.path("results","raw_cod_batch.jpg"),
       height=7, width=6, dpi=600)