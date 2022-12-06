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

pao_gao_select <- c("Ca_Accumulibacter","Tetrasphaera", "Dechloromonas",
                    "Ca_Competibacter","Micropruina","Ca_Contendobacter","Propionivibrio")

nitrifier_list <- c("Nitrotoga","Nitrospira","Nitrobacter",
                    "Nitrosomonas","Nitrosospira","Nitrosococcus",
                    "Ca_Brocadia", "Ca_Anammoxmicrobium")

# plot defaults

scale_x_main <- scale_x_date(
  breaks="2 months",
  date_labels="%b %y",
  limits = c(ymd("2021-06-01"), ymd("2022-10-01"))
)


fill_battery <- scale_fill_manual(
  limits = c("test","control"),
  labels = c("Test","Control"),
  values = met.brewer("Nizami", 2),
  name = "Battery"
)

color_battery <- scale_color_manual(
  limits = c("test","control"),
  labels = c("Test","Control"),
  values = met.brewer("Nizami", 2),
  name = "Battery"
)

shape_battery <- scale_shape_manual(
  limits = c("test","control"),
  labels = c("Test","Control"),
  values = c(16,17),
  name = "Battery"
)

scale_x_battery <- scale_x_discrete(
  limits = c("test","control"),
  labels = c("Test","Control"),
)


labels_basins <- labeller(
  location = c("control"="Control",
               "test"="Test"))

labels_paogao <- labeller(Genus = c("Ca_Accumulibacter" = "Ca. Accumulibacter",
                                  "Tetrasphaera" = "Tetrasphaera",
                                  "Dechloromonas" = "Dechloromonas",
                                  "Ca_Competibacter" = "Ca. Competibacter",
                                  "Micropruina" = "Micropruina",
                                  "Ca_Contendobacter" = "Ca. Contendobacter",
                                  "Propionivibrio" = "Propionivibrio"))

scale_x_paogao <- scale_x_discrete(labels=c("Ca_Accumulibacter" = "Ca. Accumulibacter",
                                            "Tetrasphaera" = "Tetrasphaera",
                                            "Dechloromonas" = "Dechloromonas",
                                            "Ca_Competibacter" = "Ca. Competibacter",
                                            "Micropruina" = "Micropruina",
                                            "Ca_Contendobacter" = "Ca. Contendobacter",
                                            "Propionivibrio" = "Propionivibrio"))

scale_fill_paogao <- scale_fill_manual(
  values = met.brewer("Archambault",7),
  limits = c("Ca_Accumulibacter","Tetrasphaera","Dechloromonas", 
             "Ca_Competibacter","Micropruina", "Ca_Contendobacter","Propionivibrio"),
  labels = c("Ca. Accumulibacter","Tetrasphaera","Dechloromonas",
             "Ca. Competibacter","Micropruina", "Ca. Contendobacter","Propionivibrio"))


scale_color_paogao <- scale_color_manual(
  values = met.brewer("Archambault",7),
  limits = c("Ca_Accumulibacter","Tetrasphaera","Dechloromonas", 
             "Ca_Competibacter","Micropruina", "Ca_Contendobacter","Propionivibrio"),
  labels = c("Ca. Accumulibacter","Tetrasphaera","Dechloromonas",
             "Ca. Competibacter","Micropruina", "Ca. Contendobacter","Propionivibrio"))

