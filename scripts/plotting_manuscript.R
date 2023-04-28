pao_gao_list <- c("Ca_Accumulibacter","Tetrasphaera", 
                  "Dechloromonas","Microlunatus", "Ca_Obscuribacter",
                  "Ca_Competibacter","Defluviicoccus",
                  "Micropruina","Ca_Contendobacter","Propionivibrio")

pao_gao_select <- c("Ca_Accumulibacter","Tetrasphaera","Ca_Competibacter")

nitrifier_list <- c("Nitrotoga","Nitrospira","Nitrobacter",
                    "Nitrosomonas","Nitrosospira","Nitrosococcus",
                    "Ca_Brocadia", "Ca_Anammoxmicrobium")

scale_x_main <- scale_x_date(
  breaks="3 months",
  date_labels="%b %y",
  limits = c(ymd("2021-06-15"), ymd("2022-09-10"))
)


scale_x_paogao <- scale_x_discrete(labels=c("Ca_Accumulibacter" = "Ca. Accumulibacter",
                                            "Tetrasphaera" = "Tetrasphaera",
                                            "Dechloromonas" = "Dechloromonas",
                                            "Microlunatus" = "Microlunatus",
                                            "Ca_Obscuribacter" = "Ca. Obscuribacter",
                                            "Ca_Competibacter" = "Ca. Competibacter",
                                            "Defluviicoccus" = "Defluviicoccus",
                                            "Micropruina" = "Micropruina",
                                            "Ca_Contendobacter" = "Ca. Contendobacter",
                                            "Propionivibrio" = "Propionivibrio"), limits = rev)

scale_color_paogao <- scale_color_manual(
  values = met.brewer("VanGogh2", 10),
  limits = c("Ca_Accumulibacter","Tetrasphaera","Dechloromonas", "Microlunatus", "Ca_Obscuribacter",
             "Ca_Competibacter","Defluviicoccus","Micropruina", "Ca_Contendobacter","Propionivibrio"),
  labels = c("Ca. Accumulibacter","Tetrasphaera","Dechloromonas","Microlunatus","Ca. Obscuribacter",
             "Ca. Competibacter","Defluviicoccus", "Micropruina", "Ca. Contendobacter","Propionivibrio"))

scale_fill_paogao_select <- scale_fill_manual(
  values = c("firebrick", "darkorange1", "darkolivegreen"),
  limits = c("Ca_Accumulibacter","Tetrasphaera", "Ca_Competibacter"),
  labels = c("Ca. Accumulibacter","Tetrasphaera","Ca. Competibacter"),
  name = NULL)

scale_x_nitrospira <- scale_x_discrete(
  labels=c("Nitrospira_defluvii"="Nitrospira defluvii", 
           "Nitrospira_nitrosa" ="Nitrospira nitrosa", 
           "Unknown" = "Unknown"))

scale_fill_nitrospira <- scale_fill_manual(
  values = c("dodgerblue", "dodgerblue3", "dodgerblue4"),
  limits = c("Nitrospira_defluvii", "Nitrospira_nitrosa", "Unknown"),
  labels = c("Nitrospira defluvii", "Nitrospira nitrosa", "Unknown"))
