setwd("D:/My_Documents/Dropbox/projects/Sushant/2015-08-05") # bergkamp-all

library('R.matlab')
library('arm')
library('lme4')
library('lattice')
library('languageR')
library('lmerTest')
library('DAAG')

# create df ---------------------------------------------------------------

f = readMat(Sys.glob('*trialinfo_fixedbaseline.mat'))
data = f$ALL.trialinfo
data = data[, , 1]

data$cue.ID = unlist(data$cue.ID, use.names=FALSE)

data = lapply(data, t)
data$cue.ID = t(data$cue.ID) # doesn't get done by lapply() for some reason

df <- as.data.frame(data,stringsAsFactors=FALSE)

ditch <- function(x) ifelse(is.nan(x), NA, x)
df[] <- lapply(df,ditch) # screws up factors because of ifelse()

df$cue.ID <- factor(df$cue.ID)

# inspect speed

dfs <- subset(df,block == 1 & (cue.ID == 'c1' | cue.ID == 'c5')) # only c1 and c5 for value block

boxplot(dfs$postCue.spd ~ dfs$cue.ID)

# construct some models
m1 = lm(dfs$postCue.spd ~ dfs$cue.ID)

mm_null = lmer(dfs$postCue.spd ~ (1|dfs$subject),REML=FALSE) # disable REML needed for likelihood ratio test
mm_spd = lmer(dfs$postCue.spd ~ (1|dfs$subject) + dfs$cue.ID,REML=FALSE)

anova(mm_null,mm_spd)

# with random slopes (worse fit..?)
mm_null = lmer(dfs$postCue.spd ~ (1|dfs$subject),REML=FALSE) # disable REML needed for likelihood ratio test
mm_spd = lmer(dfs$postCue.spd ~ (1 + dfs$cue.ID|dfs$subject),REML=FALSE)

anova(mm_null,mm_spd)

# better model
dfs$prevRew[is.na(dfs$.prevRew)] <- NULL;

mm_null = lmer(postCue.spd ~ (1|subject) + cue.time + xcoord + prevRew,data=dfs,REML=FALSE)
mm_spd = lmer(postCue.spd ~ (1|subject) + cue.time + xcoord + prevRew + cue.ID,data=dfs,REML=FALSE)

anova(mm_null,mm_spd)

# better model with random slopes
dfs$prevRew[is.na(dfs$.prevRew)] <- NULL;

mm_null = lmer(postCue.spd ~ (1|subject) + cue.time + xcoord + prevRew,data=dfs,REML=FALSE)
mm_spd = lmer(postCue.spd ~ cue.time + xcoord + prevRew + (cue.ID|subject),data=dfs,REML=FALSE)

anova(mm_null,mm_spd)
