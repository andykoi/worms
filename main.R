library("DescTools")
library("ggplot2")
library("geomtextpath")
library("glmmTMB")
library("DHARMa")
library("ggdist")
library("usethis")
library("effects")
library("car")
library("emmeans")
library("png")
library("ggpubr")
library("tibble")
library("grid")
library("tidytable")
library("parameters")
#install.packages("parameters")

# 1. set up data
# 2. visualizations
# 3. determine best model
# 4. use best model to run glmm
# 5. conduct ANOVA/tueky's?


##### QUESTIONS:
# - do I need to find a diff best fit model for location vs sp?
# ---- ie tweedie for sp, zigam for locations

# - what do the ANOVA and tukey tests tell us relative to the glmm output?
# ---- ie if the glmm has ** next to the items, but the anova/tukey
#       return >0.05, does that mean the results are insigificant?

##################################################
# 1 - FILE MANAGEMENT & DATA
##################################################

# Read Files for Gardner
wfish <- read.csv(header = TRUE, "dec16revis.csv")
wfish

# Renaming/Ordering/Formatting Fish Names
wfish$SPECIES[wfish$SPECIES=="Red"] <- "Sockeye"
wfish$SPECIES <- factor(wfish$SPECIES, levels=c("Pink", "Sockeye", "Chum"))
wfish$year <- c(rep(as.numeric(2024)))

wfishes <- lm(wfish$W.G.INDEX ~ unlist(wfish$SPECIES))

# Read + Format Files for Mastick
mast <- read.csv("Can Data_Jan23.csv")
mast <- mast[mast$box.no.!="Practice",]
mast <- mast[mast$can.size!="small" & mast$can.size!="med",]
mast$salmon.species[mast$salmon.species=="coho"] <- "Coho"
mast$salmon.species[mast$salmon.species=="pink"] <- "Pink"
mast$salmon.species[mast$salmon.species=="red"] <- "Sockeye"
mast$salmon.species[mast$salmon.species=="chum"] <- "Chum"
mast$salmon.species <- factor(mast$salmon.species, levels=c("Coho", "Pink", "Sockeye", "Chum"))
mast$wgindex <- as.numeric(mast$wgindex)
# Testing and removing outlier chum
zc <- (0.2750127 - mean(mast$wgindex))/sd(mast$wgindex)
zc
mast <- mast[mast$wgindex<="0.2750127",]

mast$whg <- as.numeric(mast$whg)
mast$year <- as.numeric(mast$year)

wfish$set <- c("gardner")
mast$set <- c("mastick")

combo <- data.frame(
  species <- c(mast$salmon.species, wfish$SPECIES),
  whg <- c(mast$whg, wfish$WHG),
  year <- c(mast$year, wfish$year),
  set <- c(mast$set, wfish$set)
  
)
names(combo)[names(combo) == 'species....c.mast.salmon.species..wfish.SPECIES.'] <- 'species'
names(combo)[names(combo) == 'year....c.mast.year..wfish.year.'] <- 'year'
names(combo)[names(combo) == 'whg....c.mast.whg..wfish.WHG.'] <- 'whg'
names(combo)[names(combo) == 'set....c.mast.set..wfish.set.'] <- 'set'
combo$whg <- as.numeric(combo$whg)

comboModern = combo[combo$year >= "2000",]

comboHist = combo[combo$year < "2000",]



##################################################
# 2 - VISUALIZATIONS
##################################################

########################### LOCATION

### FIGURE - Regional Comparison of 2024 Samples
im <- readPNG("akt1.png")

im2 <- matrix(rgb(im[,,1],im[,,2],im[,,3], im[,,4] * 0.5), nrow=dim(im)[1])

boxplot(wfish$WHG ~ wfish$REGION,
        main="Regional Comparison of 2024 Samples",
        xlab="", ylab="Nematodes per 100g",
        col=c("slateblue" , "royalblue", "purple"),
        border=c("slateblue4", "royalblue4", "purple4"))

ggplot(wfish, aes(x=REGION, y=WHG, fill=REGION)) +
  #background_image(akm1) +
  #geom_point(color = "red", size = 5) +
  geom_boxplot() +
  scale_fill_manual(values=c("slateblue" , "royalblue", "purple") )
#annotation_custom(rasterGrob(im2,  width = unit(1,"npc"),  height = unit(1,"npc")), -Inf, Inf, -Inf, Inf))

grid()


# VISUALIZATION - Half-eye for Regional Distribution
ggplot (wfish, aes(y = REGION, x = WHG)) +
  stat_halfeye()

########################################################


### MODEL - GLMM using ziGamma
gModelRegion = glmmTMB(
  WHG ~ REGION + SPECIES,
  data = wfish,
  ziformula = ~.,
  family=ziGamma("log"),
)
print(summary(gModelRegion),show.residuals=TRUE)

gMRanova_table_zipart = glmmTMB:::Anova.glmmTMB(gModelRegion,type=3,component="zi")
gMRanova_table_zipart

######################

gModelSpecies = glmmTMB(
  WHG ~ SPECIES + REGION,
  data = wfish,
  ziformula = ~.,
  family=ziGamma("log"),
)
print(summary(gModelSpecies),show.residuals=TRUE)

##########################################################


### RESULTS: No values under 0.05, no significant difference.


######################### SPECIES 1979-2020+2024

# VISUALZATION - Boxplot of Species for 2024
boxplot(wfish$WHG ~ wfish$SPECIES,
        main="Species Comparison of 2024 Samples",
        xlab="Species", ylab="Nematodes per 100g",
        #xlab="Species",
        #ylab="Infection Index (nematodes/gram)",
        col=c("pink" , "tomato", "#a7cdd6"),
        border=c("#f06790", "red", "#5894a3")
)

# SUMMARY - Mastick
as.numeric(mast$whg)
MPinks<-mast[mast$salmon.species=="Pink",]
MReds<-mast[mast$salmon.species=="Sockeye",]
MChum<-mast[mast$salmon.species=="Chum",]
MCoho<-mast[mast$salmon.species=="Coho",]

MPinks$whg<-as.numeric(MPinks$whg)
MReds$whg<-as.numeric(MReds$whg)
MChum$whg<-as.numeric(MChum$whg)
MCoho$whg<-as.numeric(MCoho$whg)

MStatsSum <- data.frame(
  Species = c("Pink", "Sockeye", "Chum", "Coho"),
  SampleSize = c(nrow(MPinks), nrow(MReds), nrow(MChum), nrow(MCoho)),
  Mean = c(mean(MPinks$whg), mean(MReds$whg), mean(MChum$whg), mean(MCoho$whg)),
  StDev = c(sd(MPinks$whg), sd(MReds$whg), sd(MChum$whg), sd(MCoho$whg))
)
MStatsSum

# SUMMARY - Gardner
as.numeric(wfish$WHG)
GPinks<-wfish[wfish$SPECIES=="Pink",]
GReds<-wfish[wfish$SPECIES=="Sockeye",]
GChums<-wfish[wfish$SPECIES=="Chum",]

GStatsSum <- data.frame(
  Species = c("Pink", "Sockeye", "Chum"),
  SampleSize = c(nrow(GPinks), nrow(GReds), nrow(GChums)),
  Mean = c(mean(GPinks$WHG), mean(GReds$WHG), mean(GChums$WHG)),
  StDev = c(sd(GPinks$WHG), sd(GReds$WHG), sd(GChums$WHG))
)
GStatsSum

# VISUALIZATION - Boxplot of Mastick data 
boxplot(mast$whg ~ mast$salmon.species,
        main="Species Comparison of Samples from 1982-2020\n(Mastick et al. 2024)",
        xlab="Species", ylab="Nematodes per 100g",
        col=c("gray", "pink", "tomato", "#a7cdd6"),
        border=c("darkgray", "#f06790", "red", "#5894a3"))



# VISUALIZATION - Mastick WHG over time
ggplot(mast, aes(x = year, y = whg, color = salmon.species)) +
  geom_point() +
  geom_labelsmooth(aes(label = salmon.species), fill = "white",
                   method = "lm", formula = y ~ x,
                   size = 3, linewidth = 1, boxlinewidth = 0.4) +
  xlab ("Year") +
  xlim (1979, 2020) +
  theme_bw() + guides(color = 'none')

# VISUALIZATION - Combined sets of WHG over time
#wfish$year <- c(as.numeric("2024"))



ggplot(combo, aes(x = year, y = whg, color = species)) +
  geom_point() +
  geom_labelsmooth(aes(label = species), fill = "white",
                   method = "lm", formula = y ~ x,
                   size = 3, linewidth = 1, boxlinewidth = 0.4) +
  xlab ("Year") +
  xlim (1979, 2024) +
  theme_bw() + guides(color = 'none')

#########################################################


# VISUALIZATION - Distributions of WHG per data set
# Gardner 2024
ggplot (wfish, aes(y = SPECIES, x = WHG)) +
  stat_halfeye()

# Mastick 1979-2020
ggplot (mast, aes(y = salmon.species, x = whg)) +
  stat_halfeye()

# Combined 1979-2024
ggplot (combo, aes(y = species, x = whg)) +
  stat_halfeye()
# This tells us the data is not normally distributed
# and we need to use a zero-inflated gamma model.


# TEST - GLMM Model for Species Differences

SpeciesModel = glmmTMB(  
  whg ~ species + (1|year), 
  data = combo,
  ziformula = ~.,
  family=ziGamma("log"))
print(summary(SpeciesModel),show.residuals=TRUE)

SpAIC <- AIC(SpeciesModel)
SpAIC

# VISUALIZATION - Data sets on box plot together
ggplot(combo, aes(x=species, 
                  y=whg, 
                  fill=set,
                  #colour = 'red',
)) + 
  scale_fill_manual(values = c("tomato", "#a7cdd6")) +
  geom_boxplot() +
  ylab ("Nematodes per 100g") +
  xlab ("Species") 

### FIGURE X - Combined data set box plot

boxplot(combo$whg ~ combo$species,
        #subset=(threemast$salmon.species!="Coho"),
        main="Species Comparison of Combined Samples",
        xlab="Species", ylab="Nematodes per 100g",
        col=c("gray", "pink", "tomato", "#a7cdd6"),
        border=c("darkgray", "#f06790", "red", "#5894a3"))

##############################################

### RESULTS - Zero-inflated Gamma ANOVA

AnovaSpecies = glmmTMB(
  WHG ~ SPECIES,
  data = wfish,
  ziformula = ~.,
  family=ziGamma("log"),
)

anova_table_gammapart = glmmTMB:::Anova.glmmTMB(AnovaSpecies,type=3,component="cond")
anova_table_gammapart


anova_table_zipart = glmmTMB:::Anova.glmmTMB(AnovaSpecies,type=3,component="zi")
anova_table_zipart

emm_cond = emmeans(AnovaSpecies, ~SPECIES, component="cond")
tukey_cond = pairs(emm_cond, adjust = "tukey")
tukey_cond

emm_zi = emmeans(AnovaSpecies, ~SPECIES, component="zi")
tukey_zi = pairs(emm_zi, adjust = "tukey")
tukey_zi


# Species Chisq Value 0.01942 < 0.05 


# Significant diff between sockeye + pink
# Significant diff between sockeye + coho
# No other significant differences

############# MODEL LIST


##################################################
# 3 - MODEL COMPARISONS
##################################################


# Conway-Maxwell Poisson, compois(link = "log")
modConway = glmmTMB(  
  whg ~ species + (1|year), 
  data = combo,
  family=
    compois(link = "log"))
print(summary(SpeciesModel),show.residuals=TRUE)

# tweedie
modTwee = glmmTMB(  
  whg ~ species + (1|year), 
  data = combo,
  family= tweedie(link = "log"))

# ziGamma log
modZiGam = glmmTMB(  
  whg ~ species + (1|year), 
  data = combo,
  ziformula = ~.,
  family=ziGamma("log"))

# ziGamma inverse
modZiGamIn = glmmTMB(  
  whg ~ species + (1|year), 
  data = combo,
  ziformula = ~.,
  family=ziGamma("inverse"))


# best fit

models <- list(modTwee, modZiGam)
model_names <- c("Tweedie", "ziGamma Log")


residual_summaries <- map2_df(models, model_names, function(model, name) {
  sim <- simulateResiduals(fittedModel = model, plot = FALSE)
  data.frame(
    model = name,
    dispersion = testDispersion(sim)$p.value,
    uniformity = testUniformity(sim)$p.value,
    outliers = testOutliers(sim)$p.value,
    AIC = AIC(model)
  )
})

residual_summaries %>%
  arrange(AIC) %>%
  print()

modFinal <- modTwee


simulateResiduals(fittedModel = modFinal, plot = TRUE)
summary(modFinal)
plot(parameters(modFinal))


##################################################
# 3 - FINAL MODEL USED TO EVAL VARIABLES
##################################################


gModTwee = glmmTMB(  
  WHG ~ REGION, 
  data = wfish,
  family= tweedie(link = "log"))

print(summary(gModTwee),show.residuals=TRUE)

##################################################################################


cModAnovaSpecies = glmmTMB(
  whg ~ species + (1|year),
  data = comboModern,
  ziformula = ~.,
  family=ziGamma("log"),
)
print(summary(cModAnovaSpecies),show.residuals=TRUE)

anova_table_gammapart = glmmTMB:::Anova.glmmTMB(cModAnovaSpecies,type=3,component="cond")
anova_table_gammapart


anova_table_zipart = glmmTMB:::Anova.glmmTMB(cModAnovaSpecies,type=3,component="zi")
anova_table_zipart

emm_cond = emmeans(cModAnovaSpecies, ~species, component="cond")
tukey_cond = pairs(emm_cond, adjust = "tukey")
tukey_cond

emm_zi = emmeans(cModAnovaSpecies, ~species, component="zi")
tukey_zi = pairs(emm_zi, adjust = "tukey")
tukey_zi


##################################################################################

cHistAnovaSpecies = glmmTMB(
  whg ~ species,
  data = comboHist,
  ziformula = ~.,
  family=ziGamma("log"),
)
print(summary(cHistAnovaSpecies),show.residuals=TRUE)

anova_table_gammapart = glmmTMB:::Anova.glmmTMB(cHistAnovaSpecies,type=3,component="cond")
anova_table_gammapart


anova_table_zipart = glmmTMB:::Anova.glmmTMB(cHistAnovaSpecies,type=3,component="zi")
anova_table_zipart

emm_cond = emmeans(cHistAnovaSpecies, ~species, component="cond")
tukey_cond = pairs(emm_cond, adjust = "tukey")
tukey_cond

emm_zi = emmeans(cHistAnovaSpecies, ~species, component="zi")
tukey_zi = pairs(emm_zi, adjust = "tukey")
tukey_zi


###############################################

gModTwee = glmmTMB(  
  WHG ~ SPECIES, 
  data = wfish,
  family= tweedie(link = "log"))

print(summary(gModTwee),show.residuals=TRUE)



