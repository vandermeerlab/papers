setwd("D:/My_Documents/Dropbox/projects/Sushant/2015-08-05") # bergkamp-all

library('R.matlab')
library('arm')
library('lme4')
library('lattice')
library('languageR')
library('lmerTest')
library('DAAG')
library('DAAGxtras')

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

# pre vs post cue

# normalize to pre-baseline
dfb1 <- df

dfb1$preNP.g50 <- df$preNP.g50 / df$preNP.g50
dfb1$postNP.g50 <- df$postNP.g50 / df$preNP.g50


dfb1$preNP.g80 <- df$preNP.g80 / df$preNP.g80
dfb1$postNP.g80 <- df$postNP.g80 / df$preNP.g80

# whole session
dfb2 <- df

dfb2$preNP.g50 <- df$preNP.g50 / df$baseline.g50
dfb2$postNP.g50 <- df$postNP.g50 / df$baseline.g50

dfb2$preNP.g80 <- df$preNP.g80 / df$baseline.g80
dfb2$postNP.g80 <- df$postNP.g80 / df$baseline.g80

# inspect
summary(dfb2$preNP.g50)
summary(dfb2$preNP.g80)
summary(dfb2$postNP.g50)
summary(dfb2$postNP.g80)

boxplot(dfb2$preNP.g50,dfb2$postNP.g50)
boxplot(dfb2$preNP.g80,dfb2$postNP.g80)

# tests (replicates matlab, although not with same p-values..)
(wilcox.test(dfb2$preNP.g50,dfb2$postNP.g50,paired=TRUE))
(wilcox.test(dfb2$preNP.g80,dfb2$postNP.g80,paired=TRUE))

# for linear models, need to reformat data first
g50_pwr <- c(dfb2$preNP.g50,dfb2$postNP.g50)
g80_pwr <- c(dfb2$preNP.g80,dfb2$postNP.g80)

pp <- c(rep(0,1933),rep(1,1933))
pp <- factor(pp) # before/after cue, main variable of interest

subj <- c(df$subject,df$subject)
subj <- factor(subj)

cueTime <- c(df$cue.time,df$cue.time)
preSpd <- c(df$preCue.spd,df$preCue.spd)
postSpd <- c(df$postCue.spd,df$postCue.spd)
xcoord <- c(df$xcoord,df$xcoord)
prevRew <-c(df$prevRew,df$prevRew)

npGamma_b2 <- data.frame(pp,subj,g50_pwr,g80_pwr,cueTime,preSpd,postSpd,xcoord,prevRew)

# is the mixed model better than the one without?
ml <- lm(g50_pwr ~ 1 + pp,data=npGamma_b2)
m0 <- lmer(g50_pwr ~ (1|subj),data=npGamma_b2,REML=FALSE)
m1 <- lmer(g50_pwr ~ (1|subj) + pp,data=npGamma_b2,REML=FALSE)

anova(m1,ml)

# is there an effect of cue?
anova(m0,m1)

# coefficients
m1

# try some more complicated models
m0 <- lmer(g50_pwr ~ (1|subj) + cueTime + postSpd + xcoord + prevRew,data=npGamma_b2,REML=FALSE)
m1 <- lmer(g50_pwr ~ (1|subj) + cueTime + postSpd + xcoord + prevRew + pp,data=npGamma_b2,REML=FALSE)

anova(m0,m1)

m1

# --------------- gamma-80 -------------------------------------
ml <- lm(g80_pwr ~ 1 + pp,data=npGamma_b2)
m0 <- lmer(g80_pwr ~ (1|subj),data=npGamma_b2,REML=FALSE)
m1 <- lmer(g80_pwr ~ (1|subj) + pp,data=npGamma_b2,REML=FALSE)

anova(m1,ml)

anova(m0,m1)

m1

# try some more complicated models
m0 <- lmer(g80_pwr ~ (1|subj) + cueTime + postSpd + xcoord + prevRew,data=npGamma_b2,REML=FALSE)
m1 <- lmer(g80_pwr ~ (1|subj) + cueTime + postSpd + xcoord + prevRew + pp,data=npGamma_b2,REML=FALSE)

anova(m0,m1)

m1
