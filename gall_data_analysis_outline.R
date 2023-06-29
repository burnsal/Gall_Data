## Gall Data Analysis


#########
## -- Set Up
#########

# load libraries
library(tidyverse) # data formatting and plotting
library(kableExtra) # pretty tables
library(ggsci) # colors
library(emmeans) # treatment contrasts
library(rstatix) # dplyr-friendly stat  calculations
library(pscl) # zero-inflated Poisson regression
library(lmtest) # model comparisons

# set directory if working outside RProj
#setwd("path/to/gall/data/folder")

# load data
gall_data <- readxl::read_xlsx("./data/Gall Data '23.xlsx", sheet = "Gall Data")
density_df <- readxl::read_xlsx("./data/Gall Data '23.xlsx", sheet = "Summary gall data")
insect_ids <- readxl::read_xlsx("./data/Gall Data '23.xlsx", sheet = "Insect identifications")
lpi_data <- readxl::read_xlsx("./data/Gall Data '23.xlsx", sheet = "Line-Point Intercept Data")

# overview df
str(gall_data)

# remove special characters from column names
colnames(gall_data) <- gsub(" ", "", colnames(gall_data))
colnames(gall_data) <- gsub("\\(", "_", colnames(gall_data))
colnames(gall_data) <- gsub("\\)", "", colnames(gall_data))

## - format variables
# turn pasture id, transect side, host species into factors
gall_data <- dplyr::mutate(gall_data,
                           PastureID = as.factor(PastureID),
                           Transectside = as.factor(Transectside),
                           HostSpecies = as.factor(ifelse(HostSpecies == "hybrid", "Hybrid", HostSpecies)))

# add columns for each treatment type
gall_data <- dplyr::mutate(gall_data,
                           Fire = factor(ifelse(PastureID == "1B" | PastureID == "2B" | PastureID == "EX-2B" | PastureID == "EX-1B",
                                         "Burn", "NoBurn"), levels = c("NoBurn", "Burn")),
                           Graze = factor(ifelse(PastureID == "1B" | PastureID == "2A", "Spring",
                                          ifelse(PastureID == "1A" | PastureID == "2B", "Fall", "NoGraze")), 
                                          levels = c("NoGraze", "Spring", "Fall")),
                           Treatment = factor(paste0(Graze, "_", Fire), 
                                              levels = c("NoGraze_NoBurn", "Spring_NoBurn", "Fall_NoBurn",
                                                         "NoGraze_Burn", "Spring_Burn", "Fall_Burn")))
gall_data <- dplyr::relocate(gall_data, c(Fire, Graze, Treatment), .after = PastureID)

# replace NA values with 0 in gall counts 
gall_data <- dplyr::mutate_at(gall_data, vars(DaisyGall:Greenthorn), ~replace_na(., 0))

# add col for volume of plant
gall_data <- dplyr::mutate(gall_data, PlantVol_cm3 = Width_cm * Height_cm * Cross_cm * pi/6) 
gall_data <- dplyr::relocate(gall_data, PlantVol_cm3, .after = Height_cm)

# add gall totals
galls <- dplyr::select(gall_data, c(DaisyGall:Greenthorn))
gall_data$GallTotal <- rowSums(galls)

# add galls per plant volume
gall_data <- gall_data %>%
  dplyr::filter(PlantVol_cm3 != 0) %>% # remove two rows with vol=0
  dplyr::mutate(GallperVol = GallTotal/PlantVol_cm3) # calculate galls by plant vol
  
# add plant density
density_df <- density_df %>%
  dplyr::select(c(PastureID, Transect, PlantTotal, TransectArea_m2)) %>%
  dplyr::mutate(Plants_m2 = PlantTotal / TransectArea_m2)
gall_data <- gall_data %>%
  dplyr::left_join(density_df, by = c("PastureID", "Transect"))

# check the data again
str(gall_data)

# create pivot df for plotting
gall_long_df <- gall_data %>%
  pivot_longer(cols = c(DaisyGall:Greenthorn),
               names_to = "GallType",
               values_to = "GallCount") %>%
  mutate(GallCountperVol = GallCount / PlantVol_cm3,
         GallPercent = GallCount / GallTotal,
         GallPercentperVol = GallPercent * PlantVol_cm3,
         GallCount_m2 = GallCount * Plants_m2) %>%
  mutate(across(GallPercent:GallPercentperVol, ~ replace(., is.nan(.), 0)))


########
## -- Treatment Effects on Gall Abundance
########

## How many plants have galls?

# create a binary factor for presence/absence of galls on a plant
gall_binary <- gall_data %>%
   dplyr::mutate(GallsPresent = factor(ifelse(GallTotal == 0, "No", "Yes")))

# summarise presence/absence per individual treatment, and across combination treatment
gall_binary %>%
  dplyr::group_by(Fire, GallsPresent) %>%
  dplyr::summarize(PlantCount = n()) %>%
  dplyr::mutate(Prop = PlantCount / sum(PlantCount))

gall_binary %>%
  dplyr::group_by(Graze, GallsPresent) %>%
  dplyr::summarize(PlantCount = n()) %>%
  dplyr::mutate(Prop = PlantCount / sum(PlantCount))

gall_presence <- gall_binary %>%
  dplyr::group_by(Treatment, GallsPresent) %>%
  dplyr::summarize(PlantCount = n()) %>%
  dplyr::mutate(Prop = round(PlantCount / sum(PlantCount), 3)) %>%
  dplyr::mutate(GallsPresent = ifelse(GallsPresent ==0, "No", "Yes"))

# Example of making a pretty table. 
# See kableExtra user guide for more options: https://cran.r-project.org/web/packages/kableExtra/vignettes/awesome_table_in_html.html
gall_presence %>%
  kbl(caption = "Gall Presence by Treatment") %>% # create table and give it a title
  kable_classic_2(full_width = F) %>%  # formatting
  save_kable("./viz/gall_presence_table.png")  # save it as a .png file

# make plots of presence/absence
ggplot(gall_binary, aes(x = GallsPresent, fill = Fire)) + 
  geom_bar(position = "dodge") + 
  theme_bw() + scale_fill_startrek() + 
  facet_wrap(vars(Graze)) + 
  labs(x = "Galls Present", y = "Plant Count", title = "Gall Presence by Treatment")
ggplot(gall_binary, aes(x = GallsPresent, fill = Graze)) + 
  geom_bar(position = "dodge") + 
  theme_bw() + scale_fill_startrek() + 
  facet_wrap(vars(Fire)) + 
  labs(x = "Galls Present", y = "Plant Count", title = "Gall Presence by Treatment")
ggplot(gall_binary, aes(x = GallsPresent, fill = as.factor(Transect))) + 
  geom_bar(position = "dodge") + 
  theme_bw() + scale_fill_startrek() + 
  facet_wrap(vars(Treatment)) + 
  labs(x = "Galls Present", y = "Plant Count", title = "Gall Presence by Treatment, Transect")


#---#

## Does total number of galls per plant vary by treatment?
gall_totals <- gall_data %>% 
  dplyr::group_by(Fire, Graze, Treatment) %>% 
  dplyr::summarize(GallTotal = sum(GallTotal), PlantTotal = n(), 
                   TotalPlantVol = sum(PlantVol_cm3), MeanPlantVol = mean(PlantVol_cm3), sdPlantVol = sd(PlantVol_cm3),
                   MeanPlantDensity = mean(GallperVol), sdPlantDensity = sd(GallperVol)) %>%
  dplyr::mutate(MeanGallsperPlant = GallTotal / PlantTotal) # calculate galls per plant to account for dif sample sizes

# make table of gall totals
gall_totals %>%
  select(!Treatment) %>%
  relocate(MeanGallsperPlant, .after = PlantTotal) %>%
  kbl(caption = "Gall Summary by Treatment") %>%
  kable_minimal(full_width = F) %>%
  save_kable("./viz/gall_summary_table.png")

# perform Kruskal-Wallis rank sum test for significant differences in gall totals by treatment
# Note: Kruskal-Wallis test is a non-parametric test comparable to ANOVA, but does not make any underlying assumptions (ie linearity, normality) about the data
kruskal.test(GallTotal ~ Fire, data = gall_totals)
kruskal.test(GallTotal ~ Graze, data = gall_totals)
kruskal.test(GallTotal ~ Treatment, data = gall_totals)
kruskal.test(MeanGallsperPlant ~ Fire, data = gall_totals)
kruskal.test(MeanGallsperPlant ~ Graze, data = gall_totals)
kruskal.test(MeanGallsperPlant ~ Treatment, data = gall_totals)
# no significant differences identified

# frequency plot of gall totals by treatment   
ggplot(gall_data, aes(x = GallTotal, fill = Treatment)) + 
  geom_histogram(binwidth = 10, alpha = 0.5) + 
  theme_bw() + 
  labs(x = "Gall Total", y = "Plant Count", title = "Gall Total per Plant, by Treatment")

# density plot of gall totals by treatment (to account for different sample sizes)   
ggplot(gall_data, aes(x = GallTotal, after_stat(density), color = Treatment)) + 
  geom_freqpoly(binwidth = 50, lwd = 1.2, alpha = 0.5) + 
  theme_bw() + 
  labs(x = "Gall Total", y = "Plant Count Density", title = "Gall Total per Plant")

#---#

## Do gall totals vary between transects within treatment?
galltotals_transect <- gall_data %>%
  dplyr::select(c(Fire, Graze, Treatment, Transect, Transectside, PlantVol_cm3, GallTotal, GallperVol)) %>%
  group_by(Treatment, Transect) %>%  # you could add Transectside as a grouping variable to check if transect side is important
  dplyr::summarise(meanPlantVol = mean(PlantVol_cm3), 
                   meanGalls = mean(GallTotal), 
                   meanGallperVol = mean(GallperVol)) %>%
  dplyr::arrange(Treatment, Transect)

# Test for differences of gall totals within treatment between transects
gall_data %>% group_by(Treatment) %>% rstatix::kruskal_test(GallTotal ~ Transect) # Note: p-value of true significance is .05 / 6 = .0083 (6 levels)
gall_data %>% group_by(Fire) %>% rstatix::kruskal_test(GallTotal ~ Transect) # Note: p-value of true significance is .05 / 2 = .025 (2 levels)
gall_data %>% group_by(Graze) %>% rstatix::kruskal_test(GallTotal ~ Transect) # Note: p-value of true significance is .05 / 3 = .017 (3 levels)

# Test for differences of gall per plant volumne within treatment between transects
# Note: p-value of true significance is .05 / 6 = .0083 (6 treatments)
gall_data %>% group_by(Treatment) %>% rstatix::kruskal_test(GallperVol ~ Transect)
gall_data %>% group_by(Fire) %>% rstatix::kruskal_test(GallperVol ~ Transect)
gall_data %>% group_by(Graze) %>% rstatix::kruskal_test(GallperVol ~ Transect)

# patterns look similar between measures of total galls and galls per volume

# visualize gall totals within treatment by transect
ggplot(gall_data, aes(x = Transect, y = GallTotal, fill = Transectside)) + 
  geom_col() + 
  facet_wrap(vars(Treatment)) +
  scale_fill_aaas() + 
  theme_bw() + 
  ggtitle("Gall Total By Treatment, Transect")
# most galls come from east side of transects

#---#

# get number of gall types
n_galls <- length(unique(gall_long_df$GallType))

# average galls per plant by gall type, treatment
gall_type_counts <- gall_long_df %>%
  dplyr::group_by(Fire, Graze, GallType) %>%
  dplyr::summarize(TotalCount = sum(GallCount), MeanCount = mean(GallCount), sdCount = sd(GallCount),
                   TotalDensity = sum(GallCount)/sum(PlantVol_cm3), MeanDensity = mean(GallCountperVol), sdDensity = sd(GallCountperVol),
                   MeanPercent = mean(GallPercent), sdPercent = sd(GallPercent), 
                   MeanPercentperVol = mean(GallPercentperVol), sdPercentperVol = sd(GallPercentperVol))

gall_type_counts %>%
  kbl(caption = "Gall Type summary Table") %>%
  kable_classic_2() %>%
  save_kable("./viz/galltype_summary_table.png")

# make table of gall counts by type and treatment 
galltypes <- gall_long_df %>% 
  group_by(GallType, Treatment) %>% 
  summarize(GallTotals = sum(GallCount)) %>%
  pivot_wider(names_from = "Treatment", values_from = "GallTotals")

galltypes_mat <- as.matrix(galltypes[,2:7])
galltypes$TypeTotal <- rowSums(galltypes_mat)
TreatmentTotal <- data.frame(t(c("TreatmentTotal", colSums(galltypes_mat), sum(rowSums(galltypes_mat)))))
TreatmentTotal[,2:8] <- as.numeric(TreatmentTotal[,2:8])
colnames(TreatmentTotal) <- colnames(galltypes)
galltypes <- bind_rows(galltypes, TreatmentTotal)

galltypes %>%
  kbl(caption = "Gall Count by Gall Type, Treatment",
      col.names = c("Gall Type", rep(c("No Graze", "Spring", "Fall"), times = 2), "Total")) %>%
  kable_classic(full_width = F, html_font = "Cambria") %>%
  add_header_above(c(" " = 1, "No Burn" = 3, "Burn" = 3, " " = 1)) %>%
  save_kable("./viz/galltype_count_table.png")

# gall density by type and treatment 
galldensity <- gall_long_df %>% 
  group_by(GallType, Treatment) %>% 
  summarize(GallDensity = sum(GallCountperVol)) %>%
  pivot_wider(names_from = "Treatment", values_from = "GallDensity")

galldensity_mat <- as.matrix(galldensity[,2:7])
galldensity$TypeTotal <- rowSums(galldensity_mat)
TreatmentTotal <- data.frame(t(c("TreatmentTotal", colSums(galldensity_mat), sum(rowSums(galldensity_mat)))))
TreatmentTotal[,2:8] <- as.numeric(TreatmentTotal[,2:8])
colnames(TreatmentTotal) <- colnames(galldensity)
galldensity <- bind_rows(galldensity, TreatmentTotal)

galldensity %>%
  kbl(caption = "Gall Density (Galls per Plant Volume) by Gall Type, Treatment",
      col.names = c("Gall Type", rep(c("No Graze", "Spring", "Fall"), times = 2), "Total")) %>%
  kable_classic(full_width = F, html_font = "Cambria") %>%
  add_header_above(c(" " = 1, "No Burn" = 3, "Burn" = 3, " " = 1)) %>%
  save_kable("./viz/galltype_density_table.png")

# galltype by percent on plant, treatment
gallpercents <- gall_long_df %>% 
  group_by(GallType, Treatment) %>% 
  summarize(PercentCt = mean(GallPercent)) %>%
  mutate(PercentCt = round(PercentCt, 3)) %>%
  pivot_wider(names_from = "Treatment", values_from = PercentCt)
gallpercents %>%
  kbl(caption = "Mean Proportion of Galls by Gall Type, Treatment",
      col.names = c("Gall Type", rep(c("No Graze", "Spring", "Fall"), times = 2))) %>%
  kable_classic(full_width = F, html_font = "Cambria") %>%
  add_header_above(c(" " = 1, "No Burn" = 3, "Burn" = 3)) %>%
  save_kable("./viz/galltypepercentsbytrt_table.png")

# galltype by percent density on plant, treatment
gallpctdens <- gall_long_df %>% 
  group_by(GallType, Treatment) %>% 
  summarize(Density = mean(GallPercentperVol)) %>%
  pivot_wider(names_from = "Treatment", values_from = Density)
gallpctdens %>%
  kbl(caption = "Mean Percentage Gall Density by Gall Type, Treatment",
      col.names = c("Gall Type", rep(c("No Graze", "Spring", "Fall"), times = 2))) %>%
  kable_classic(full_width = F, html_font = "Cambria") %>%
  add_header_above(c(" " = 1, "No Burn" = 3, "Burn" = 3)) %>%
  save_kable("./viz/galltypepercentdensbytrt_table.png")


# playing with colors and themes in plots
# gall counts by treatment
ggplot(gall_long_df, aes(x = Graze, y = GallCount, fill = GallType)) + 
  geom_col() + 
  facet_grid(cols = vars(Fire), scales = "free_x") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust = 1)) +
  scale_fill_ucscgb() +
  ggtitle("Total Galls by Treatment and Gall Type")

ggplot(gall_long_df, aes(x = Graze, y = GallPercentperVol, fill = GallType)) + 
  geom_col() + 
  facet_grid(cols = vars(Fire), scales = "free_x") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust = 1)) +
  scale_fill_ucscgb() +
  labs(x = "Treament", y = "Gall Count per cm^3")+ 
  ggtitle("Gall per Plant Density by Treatment and Gall Type")


ggplot(gall_long_df, aes(x = Treatment, y = GallCount, fill = GallType)) + 
  geom_col() + 
  facet_grid(cols = vars(Graze), scales = "free_x") + 
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust = 1)) +
  scale_fill_ucscgb() +
  ggtitle("Total Galls by Treatment and Gall Type")

# call counts by fire treatment
ggplot(gall_long_df, aes(x = Fire, y = GallCount, fill = GallType)) + 
  geom_col() + 
  theme_minimal() +
  scale_fill_igv(palette = "default") +
  ggtitle("Total Galls per Fire Treatment, by Gall Type")

# gall counts by graze treatment
ggplot(gall_long_df, aes(x = Graze, y = GallCount, fill = GallType)) + 
  geom_col() + 
  theme_light() +
  scale_fill_d3(palette = "category20b") +
  ggtitle("Total Galls per Graze Treatment, by Gall Type")

# gall counts by treatment combo
ggplot(gall_long_df, aes(x = Treatment, y = GallCount, fill = GallType)) + 
  geom_col() + 
  theme_light() +
  scale_fill_d3(palette = "category20b") +
  theme(axis.text.x = element_text(angle=45, hjust = 1)) + 
  ggtitle("Total Galls per Graze:Fire Treatment, by Gall Type")


# look at average gall per plant volume by treatments
gall_data %>%
  group_by(Fire, Graze) %>%
  filter(PlantVol_cm3 != 0) %>%
  summarize(avg.gall = mean(GallperVol))

ggplot(gall_data, (aes(x = Fire, y = GallperVol, fill = Graze))) + 
  geom_violin(draw_quantiles = TRUE) + 
  ggtitle("Total Galls per Plant Volume by Treatment")


########
## -- Effects of Available Host Plant Material on Gall Abundance
########

# Are there more galls on larger plants?
# look at correlation between galltotal and plant volume
cor(gall_data$PlantVol_cm3, gall_data$GallTotal)
cor(gall_data$PlantVol_cm3, gall_data$GallTotal, method = "kendall")
cor(gall_data$PlantVol_cm3, gall_data$GallTotal, method = "spearman")
# seems like there is a moderate positive correlation

# Are galls more dense on larger plants?
cor(gall_data$GallperVol, gall_data$PlantVol_cm3)
cor(gall_data$PlantVol_cm3, gall_data$GallperVol, method = "kendall")
cor(gall_data$PlantVol_cm3, gall_data$GallperVol, method = "spearman")
# no correlation of note

# summary table of galls by plant volume
gall_data %>%
  group_by(Fire, Graze) %>%
  

# visualize relationship by Treatment
ggplot(gall_data, aes(x = PlantVol_cm3, y=GallTotal)) + 
  geom_point(aes(col = Treatment), size = 3, alpha = 0.6) + 
  theme_minimal() + lims(x = c(0, 300000), y = c(0, 300)) + 
  labs(x = expression("Plant Volume cm"^3), title = "Gall Total by Plant Volume, Treatment")

# visualize relationship by Treatment, Gall Type
ggplot(gall_long_df, aes(x = PlantVol_cm3, y=GallPercent)) + 
  geom_point(aes(col = Treatment), size = 3, alpha = 0.6) + 
  scale_color_nejm() +
  theme_bw() + lims(x = c(0, 300000), y = c(0.001, 1)) + 
  facet_wrap(vars(GallType)) + 
  theme(axis.text.x = element_text(angle=45, hjust = 1)) + 
  labs(x = expression("Plant Volume cm"^3), title = "Gall Percentage by Plant Volume, Treatment, Gall Type \n Excludes 0-Counts")

# visualize relationship by Treatment, Gall Type
ggplot(gall_long_df, aes(x = PlantVol_cm3, y=GallPercent)) + 
  geom_point(aes(col = GallType), size = 3, alpha = 0.6) + 
  scale_color_igv() +
  theme_bw() + lims(x = c(0, 300000), y = c(.001, 1)) + 
  facet_wrap(vars(Treatment)) + 
  theme(axis.text.x = element_text(angle=45, hjust = 1)) + 
  labs(x = expression("Plant Volume cm"^3), title = "Gall Percentage by Plant Volume, Treatment, Gall Type \n Excludes 0-Counts")

ggplot(gall_long_df, aes(x = Graze, y = GallCountperVol, fill = GallType)) + 
  geom_col() + facet_grid(cols = vars(Fire)) + 
  scale_fill_ucscgb() +
  theme_bw() + 
  ggtitle("Total Galls per Plant Volume by Fire, Graze Treatments and Gall Type")

ggplot(gall_long_df, aes(x = Graze, y = GallPercent, fill = GallType)) + 
  geom_col() + facet_grid(cols = vars(Fire)) + 
  ggtitle("Total Gall Counts by Fire, Graze Treatments and Gall Type")



########
## -- Effects of Plant Community on Gall Abundance
########

# Just look at plant density by Treatment
plant_density <- gall_data %>%
  select(c(Fire, Graze, Transect, PlantVol_cm3, PlantTotal, Plants_m2)) %>%
  distinct()

ggplot(plant_density, aes(x = Plants_m2, y = PlantTotal)) + 
  geom_point(aes(col = Graze), size = 3, alpha = 0.5) + facet_wrap(vars(Fire)) + 
  theme_bw() + ggtitle("Plant Total by Plant Density, Treatment")
ggplot(plant_density, aes(x = Plants_m2, y = PlantTotal)) + 
  geom_point(aes(col = Fire), size = 3, alpha = 0.5) + facet_wrap(vars(Graze)) + 
  theme_bw() + ggtitle("Plant Total by Plant Density, Treatment")

ggplot(plant_density, aes(x = Plants, y = Plants_m2)) + 
  geom_point(aes(col = Fire), size = 3, alpha = 0.5) + 
  facet_wrap(vars(Graze)) + 
  theme_bw() + ggtitle("Plant Density by Plant Volume, Treatment")

## -- Galls per Plant vs Plant Density
galltype_plant <- gall_long_df %>%
  group_by(GallType, Treatment) %>%
  summarize(meanGallsper_m2 = mean(GallCount_m2),
            sdGallsper_m2 = sd(GallCount_m2))

# create a table of mean counts
galltype_plant %>%
  select(!sdGallsper_m2) %>%
  mutate(meanGallsper_m2 = round(meanGallsper_m2, 2)) %>%
  pivot_wider(names_from = "Treatment", values_from = "meanGallsper_m2") %>%
  kbl(caption = "Mean Gall Count per Plant Density by Gall Type, Treatment",
      col.names = c("Gall Type", rep(c("No Graze", "Spring", "Fall"), times = 2))) %>%
  kable_classic(full_width = F, html_font = "Cambria") %>%
  add_header_above(c(" " = 1, "No Burn" = 3, "Burn" = 3)) %>%
  save_kable("./viz/galltypeplantdens_table.png")

########
## -- Modeling Gall Abundance by Treatment, Gall Type
########


# Estimate gall abundance using a zero-inflated poisson regression model
mod0 <- pscl::zeroinfl(GallCount ~ Fire * Graze + Plants_m2 + GallType, data = gall_long_df, dist = "poisson")
mod1 <- pscl::zeroinfl(GallCount ~ Fire * Graze + Plants_m2, data = gall_long_df, dist = "poisson")
mod2 <- pscl::zeroinfl(GallCount ~ Fire * Graze, data = gall_long_df, dist = "poisson")

lmtest(mod2, mod1, mod0)

# TODO: plot models on top of empirical data

