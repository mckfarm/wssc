## ---------------------------
## Script name: plotting.R
## Author: McKenna Farmer
## Date Created: 2022-09-14
## ---------------------------
## Notes:
## copied from Calumet 
##
## ---------------------------

# microbe screening
pao_gao_list <- c("Ca_Accumulibacter","Tetrasphaera", 
                  "Dechloromonas","Microlunatus", "Ca_Obscuribacter",
                  "Ca_Competibacter","Defluviicoccus",
                  "Micropruina","Ca_Contendobacter","Propionivibrio")

pao_gao_select <- c("Ca_Accumulibacter","Tetrasphaera","Ca_Competibacter")

nitrifier_list <- c("Nitrotoga","Nitrospira","Nitrobacter",
                    "Nitrosomonas","Nitrosospira","Nitrosococcus",
                    "Ca_Brocadia", "Ca_Anammoxmicrobium")


# plot defaults
phases <- data.frame(x1=ymd("2022-03-22"),
                     x2=ymd("2022-04-07"))

pause_lines <- list(
  geom_vline(xintercept=phases$x1, color="grey"),
  geom_vline(xintercept=phases$x2, color="grey"))


scale_x_main <- scale_x_date(
  breaks="2 months",
  date_labels="%b %y",
  limits = c(ymd("2021-06-15"), ymd("2022-09-10"))
)


fill_battery <- scale_fill_manual(
  limits = c("test","control"),
  labels = c("Test","Control"),
  values = c("red","blue"),
  name = "Battery"
)

color_battery <- scale_color_manual(
  limits = c("test","control"),
  labels = c("Test","Control"),
  values = c("red","blue"),
  name = "Battery"
)

shape_battery <- scale_shape_manual(
  limits = c("test","control"),
  labels = c("Test","Control"),
  values = c(16,17),
  name = "Battery"
)

color_battery_gly <- scale_color_manual(
  limits = c("test","control","RAS"),
  labels = c("Test","Control","RAS"),
  values = c("red","blue","darkgreen"),
  name = "Battery"
)

shape_battery_gly <- scale_shape_manual(
  limits = c("test","control","RAS"),
  labels = c("Test","Control","RAS"),
  values = c(16,17,18),
  name = "Battery"
)

scale_x_battery <- scale_x_discrete(
  limits = c("test","control"),
  labels = c("Test","Control"),
)

labels_basins <- labeller(battery=
  c("control"="Control",
    "test"="Test"))

labels_paogao <- labeller(Genus = c("Ca_Accumulibacter" = "Ca. Accumulibacter",
                                  "Tetrasphaera" = "Tetrasphaera",
                                  "Dechloromonas" = "Dechloromonas",
                                  "Microlunatus" = "Microlunatus",
                                  "Ca_Obscuribacter" = "Ca. Obscuribacter",
                                  "Ca_Competibacter" = "Ca. Competibacter",
                                  "Defluviicoccus" = "Defluviicoccus",
                                  "Micropruina" = "Micropruina",
                                  "Ca_Contendobacter" = "Ca. Contendobacter",
                                  "Propionivibrio" = "Propionivibrio"))

scale_x_paogao <- scale_x_discrete(labels=c("Ca_Accumulibacter" = "Ca. Accumulibacter",
                                            "Tetrasphaera" = "Tetrasphaera",
                                            "Dechloromonas" = "Dechloromonas",
                                            "Microlunatus" = "Microlunatus",
                                            "Ca_Obscuribacter" = "Ca. Obscuribacter",
                                            "Ca_Competibacter" = "Ca. Competibacter",
                                            "Defluviicoccus" = "Defluviicoccus",
                                            "Micropruina" = "Micropruina",
                                            "Ca_Contendobacter" = "Ca. Contendobacter",
                                            "Propionivibrio" = "Propionivibrio"))


scale_fill_paogao <- scale_fill_manual(
  values = met.brewer("Redon",10),
  limits = c("Ca_Accumulibacter","Tetrasphaera","Dechloromonas", "Microlunatus", "Ca_Obscuribacter",
             "Ca_Competibacter","Defluviicoccus", "Micropruina", "Ca_Contendobacter","Propionivibrio"),
  labels = c("Ca. Accumulibacter","Tetrasphaera","Dechloromonas","Microlunatus","Ca. Obscuribacter",
             "Ca. Competibacter","Defluviicoccus", "Micropruina", "Ca. Contendobacter","Propionivibrio"))

scale_color_paogao <- scale_color_manual(
  values = met.brewer("Redon",10),
  limits = c("Ca_Accumulibacter","Tetrasphaera","Dechloromonas", "Microlunatus", "Ca_Obscuribacter",
             "Ca_Competibacter","Defluviicoccus","Micropruina", "Ca_Contendobacter","Propionivibrio"),
  labels = c("Ca. Accumulibacter","Tetrasphaera","Dechloromonas","Microlunatus","Ca. Obscuribacter",
             "Ca. Competibacter","Defluviicoccus", "Micropruina", "Ca. Contendobacter","Propionivibrio"))

scale_fill_paogao_select <- scale_fill_manual(
  values = met.brewer("Cross",3),
  limits = c("Ca_Accumulibacter","Tetrasphaera", "Ca_Competibacter"),
  labels = c("Ca. Accumulibacter","Tetrasphaera","Ca. Competibacter"))

scale_color_paogao_select <- scale_color_manual(
  values = met.brewer("Cassatt1",3),
  limits = c("Ca_Accumulibacter","Tetrasphaera", "Ca_Competibacter"),
  labels = c("Ca. Accumulibacter","Tetrasphaera","Ca. Competibacter"))


