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

# NP 1p vd 5p


# NORMALIZATION -----------------------------------------------------------


# normalize to pre-baseline
dfb1 <- df

dfb1$preNP.g50 <- df$preNP.g50 / df$preNP.g50
dfb1$postNP.g50 <- df$postNP.g50 / df$preNP.g50


dfb1$preNP.g80 <- df$preNP.g80 / df$preNP.g80
dfb1$postNP.g80 <- df$postNP.g80 / df$preNP.g80

# whole session baseline
dfb2 <- df

dfb2$preNP.g50 <- df$preNP.g50 / df$baseline.g50
dfb2$postNP.g50 <- df$postNP.g50 / df$baseline.g50

dfb2$preNP.g80 <- df$preNP.g80 / df$baseline.g80
dfb2$postNP.g80 <- df$postNP.g80 / df$baseline.g80

# subset
dfb1s <- subset(dfb1,block == 1 & (cue.ID == 'c1' | cue.ID == 'c5'))
dfb2s <- subset(dfb2,block == 1 & (cue.ID == 'c1' | cue.ID == 'c5'))


# tests (replicates matlab, although not with same p-values..)
(wilcox.test(dfb1s$postNP.g80[which(df$cue.ID == "c1")],dfb1s$postNP.g80[which(df$cue.ID == "c5")]))
(wilcox.test(dfb1s$postNP.g50[which(df$cue.ID == "c1")],dfb1s$postNP.g50[which(df$cue.ID == "c5")]))
(wilcox.test(dfb2s$postNP.g80[which(df$cue.ID == "c1")],dfb2s$postNP.g80[which(df$cue.ID == "c5")]))
(wilcox.test(dfb2s$postNP.g50[which(df$cue.ID == "c1")],dfb2s$postNP.g50[which(df$cue.ID == "c5")]))


# G50 norm to preNP -------------------------------------------------------

# is the mixed model better than the one without?
ml <- lm(postNP.g50 ~ 1 + cue.ID,data=dfb1s)
m0 <- lmer(postNP.g50 ~ (1|subject),data=dfb1s,REML=FALSE)
m1 <- lmer(postNP.g50 ~ (1|subject) + cue.ID,data=dfb1s,REML=FALSE)

anova(m1,ml) # mixed model is better

# is there an effect of cue ID?
anova(m0,m1)

# try a more complicated model
m0 <- lmer(postNP.g50 ~ (1|subject) + prevRew + cue.time,data=dfb1s,REML=FALSE)
m1 <- lmer(postNP.g50 ~ (1|subject) + prevRew + cue.time + cue.ID,data=dfb1s,REML=FALSE)

anova(m0,m1)

# G50 norm to full session baseline  -------------------------------------------------------

# is the mixed model better than the one without?
ml <- lm(postNP.g50 ~ 1 + cue.ID,data=dfb2s)
m0 <- lmer(postNP.g50 ~ (1|subject),data=dfb2s,REML=FALSE)
m1 <- lmer(postNP.g50 ~ (1|subject) + cue.ID,data=dfb2s,REML=FALSE)

anova(m1,ml) # mixed model is better

# is there an effect of cue ID?
anova(m0,m1)

# try a more complicated model
m0 <- lmer(postNP.g50 ~ (1|subject) + prevRew + cue.time + preNP.g50,data=dfb2s,REML=FALSE)
m1 <- lmer(postNP.g50 ~ (1|subject) + prevRew + cue.time + preNP.g50 + cue.ID,data=dfb2s,REML=FALSE)

anova(m0,m1)

m0 <- lmer(postNP.g50 ~ (1|subject) + prevRew + cue.time,data=dfb2s,REML=FALSE)
m1 <- lmer(postNP.g50 ~ (1|subject) + prevRew + cue.time + cue.ID,data=dfb2s,REML=FALSE)

anova(m0,m1)

# G80 norm to preNP -------------------------------------------------------

# is the mixed model better than the one without?
ml <- lm(postNP.g80 ~ 1 + cue.ID,data=dfb1s)
m0 <- lmer(postNP.g80 ~ (1|subject),data=dfb1s,REML=FALSE)
m1 <- lmer(postNP.g80 ~ (1|subject) + cue.ID,data=dfb1s,REML=FALSE)

anova(m1,ml) # mixed model is better

# is there an effect of cue ID?
anova(m0,m1)

# try a more complicated model
m0 <- lmer(postNP.g80 ~ (1|subject) + prevRew + cue.time,data=dfb1s,REML=FALSE)
m1 <- lmer(postNP.g80 ~ (1|subject) + prevRew + cue.time + cue.ID,data=dfb1s,REML=FALSE)

anova(m0,m1)

# G80 norm to full session baseline  -------------------------------------------------------

# is the mixed model better than the one without?
ml <- lm(postNP.g80 ~ 1 + cue.ID,data=dfb2s)
m0 <- lmer(postNP.g80 ~ (1|subject),data=dfb2s,REML=FALSE)
m1 <- lmer(postNP.g80 ~ (1|subject) + cue.ID,data=dfb2s,REML=FALSE)

anova(m1,ml) # mixed model is better

# is there an effect of cue ID?
anova(m0,m1)

# try a more complicated model
m0 <- lmer(postNP.g80 ~ (1|subject) + prevRew + cue.time,data=dfb2s,REML=FALSE)
m1 <- lmer(postNP.g80 ~ (1|subject) + prevRew + cue.time + cue.ID,data=dfb2s,REML=FALSE)

anova(m0,m1)

m0 <- lmer(postNP.g80 ~ (1|subject) + prevRew + preNP.g80 + cue.time,data=dfb2s,REML=FALSE)
m1 <- lmer(postNP.g80 ~ (1|subject) + prevRew + preNP.g80 + cue.time + cue.ID,data=dfb2s,REML=FALSE)

anova(m0,m1)

# Verify we can still detect overall effect of cue... (preNP baseline)

g50_pwr <- c(dfb1s$preNP.g50,dfb1s$postNP.g50)
g80_pwr <- c(dfb1s$preNP.g80,dfb1s$postNP.g80)

pp <- c(rep(0,955),rep(1,955))
pp <- factor(pp) # before/after cue, main variable of interest

subj <- c(dfb1s$subject,dfb1s$subject)
subj <- factor(subj)

cueTime <- c(dfb1s$cue.time,dfb1s$cue.time)
preSpd <- c(dfb1s$preCue.spd,dfb1s$preCue.spd)
postSpd <- c(dfb1s$postCue.spd,dfb1s$postCue.spd)
xcoord <- c(dfb1s$xcoord,dfb1s$xcoord)
prevRew <-c(dfb1s$prevRew,dfb1s$prevRew)

npGamma_b1 <- data.frame(pp,subj,g50_pwr,g80_pwr,cueTime,preSpd,postSpd,xcoord,prevRew)

# is the mixed model better than the one without?
ml <- lm(g50_pwr ~ 1 + pp,data=npGamma_b1)
m0 <- lmer(g50_pwr ~ (1|subj),data=npGamma_b1,REML=FALSE)
m1 <- lmer(g50_pwr ~ (1|subj) + pp,data=npGamma_b1,REML=FALSE)

anova(m1,ml)

# is there an effect of cue?
anova(m0,m1)

# coefficients
m1

# try some more complicated models
m0 <- lmer(g50_pwr ~ (1|subj) + cueTime + postSpd + xcoord + prevRew,data=npGamma_b1,REML=FALSE)
m1 <- lmer(g50_pwr ~ (1|subj) + cueTime + postSpd + xcoord + prevRew + pp,data=npGamma_b1,REML=FALSE)

anova(m0,m1)

m1

## gamma-80

# is the mixed model better than the one without?
ml <- lm(g80_pwr ~ 1 + pp,data=npGamma_b1)
m0 <- lmer(g80_pwr ~ (1|subj),data=npGamma_b1,REML=FALSE)
m1 <- lmer(g80_pwr ~ (1|subj) + pp,data=npGamma_b1,REML=FALSE)

anova(m1,ml)

# is there an effect of cue?
anova(m0,m1)

# coefficients
m1

# try some more complicated models
m0 <- lmer(g80_pwr ~ (1|subj) + cueTime + postSpd + xcoord + prevRew,data=npGamma_b1,REML=FALSE)
m1 <- lmer(g80_pwr ~ (1|subj) + cueTime + postSpd + xcoord + prevRew + pp,data=npGamma_b1,REML=FALSE)

anova(m0,m1)

m1


# Verify we can still detect overall effect of cue... (full baseline)

g50_pwr <- c(dfb2s$preNP.g50,dfb2s$postNP.g50)
g80_pwr <- c(dfb2s$preNP.g80,dfb2s$postNP.g80)

pp <- c(rep(0,955),rep(1,955))
pp <- factor(pp) # before/after cue, main variable of interest

subj <- c(dfb1s$subject,dfb1s$subject)
subj <- factor(subj)

cueTime <- c(dfb1s$cue.time,dfb1s$cue.time)
preSpd <- c(dfb1s$preCue.spd,dfb1s$preCue.spd)
postSpd <- c(dfb1s$postCue.spd,dfb1s$postCue.spd)
xcoord <- c(dfb1s$xcoord,dfb1s$xcoord)
prevRew <-c(dfb1s$prevRew,dfb1s$prevRew)

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

## gamma-80

# is the mixed model better than the one without?
ml <- lm(g80_pwr ~ 1 + pp,data=npGamma_b2)
m0 <- lmer(g80_pwr ~ (1|subj),data=npGamma_b2,REML=FALSE)
m1 <- lmer(g80_pwr ~ (1|subj) + pp,data=npGamma_b2,REML=FALSE)

anova(m1,ml)

# is there an effect of cue?
anova(m0,m1)

# coefficients
m1

# try some more complicated models
m0 <- lmer(g80_pwr ~ (1|subj) + cueTime + postSpd + xcoord + prevRew,data=npGamma_b2,REML=FALSE)
m1 <- lmer(g80_pwr ~ (1|subj) + cueTime + postSpd + xcoord + prevRew + pp,data=npGamma_b2,REML=FALSE)

anova(m0,m1)

m1