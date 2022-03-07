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

# pre vs post cue

# normalize to pre-baseline
dfb1 <- df

dfb1$preCue.g50 <- df$preCue.g50 / df$preCue.g50
dfb1$postCue.g50 <- df$postCue.g50 / df$preCue.g50


dfb1$preCue.g80 <- df$preCue.g80 / df$preCue.g80
dfb1$postCue.g80 <- df$postCue.g80 / df$preCue.g80

# whole session
dfb2 <- df

dfb2$preCue.g50 <- df$preCue.g50 / df$baseline.g50
dfb2$postCue.g50 <- df$postCue.g50 / df$baseline.g50

dfb2$preCue.g80 <- df$preCue.g80 / df$baseline.g80
dfb2$postCue.g80 <- df$postCue.g80 / df$baseline.g80

# inspect
summary(dfb1$postCue.g50)
summary(dfb1$postCue.g80)

boxplot(dfb1$preCue.g50,dfb1$postCue.g50)
boxplot(dfb1$preCue.g80,dfb1$postCue.g80)

# tests (replicates matlab, although not with same p-values..)
(wilcox.test(dfb1$preCue.g50,dfb1$postCue.g50,paired=TRUE))
(wilcox.test(dfb1$preCue.g80,dfb1$postCue.g80,paired=TRUE))

# for linear models, need to reformat data first
g50_pwr <- c(dfb1$preCue.g50,dfb1$postCue.g50)
g80_pwr <- c(dfb1$preCue.g80,dfb1$postCue.g80)

pp <- c(rep(0,1933),rep(1,1933))
pp <- factor(pp) # before/after cue, main variable of interest

subj <- c(df$subject,df$subject)
subj <- factor(subj)

cueTime <- c(df$cue.time,df$cue.time)
preSpd <- c(df$preCue.spd,df$preCue.spd)
postSpd <- c(df$postCue.spd,df$postCue.spd)
xcoord <- c(df$xcoord,df$xcoord)
prevRew <-c(df$prevRew,df$prevRew)

cueGamma_b1 <- data.frame(pp,subj,g50_pwr,g80_pwr,cueTime,preSpd,postSpd,xcoord,prevRew)

# is the mixed model better than the one without?
ml <- lm(g50_pwr ~ 1 + pp,data=cueGamma_b1)
m0 <- lmer(g50_pwr ~ (1|subj),data=cueGamma_b1,REML=FALSE)
m1 <- lmer(g50_pwr ~ (1|subj) + pp,data=cueGamma_b1,REML=FALSE)

anova(m1,ml)

# is there an effect of cue?
anova(m0,m1)

# coefficients
m1

# try some more complicated models
m0 <- lmer(g50_pwr ~ (1|subj) + cueTime,data=cueGamma_b1,REML=FALSE,control=lmerControl(optimizer="bobyqa"))
m1 <- lmer(g50_pwr ~ (1|subj) + cueTime + pp,data=cueGamma_b1,REML=FALSE,control=lmerControl(optimizer="bobyqa"))

anova(m0,m1)

# --------------- gamma-80 -------------------------------------
ml <- lm(g80_pwr ~ 1 + pp,data=cueGamma_b1)
m0 <- lmer(g80_pwr ~ (1|subj),data=cueGamma_b1,REML=FALSE)
m1 <- lmer(g80_pwr ~ (1|subj) + pp,data=cueGamma_b1,REML=FALSE)

anova(m1,ml)

anova(m0,m1)

m1

# try some more complicated models
m0 <- lmer(g80_pwr ~ (1|subj) + cueTime,data=cueGamma_b1,REML=FALSE,control=lmerControl(optimizer="bobyqa"))
m1 <- lmer(g80_pwr ~ (1|subj) + cueTime + pp,data=cueGamma_b1,REML=FALSE,control=lmerControl(optimizer="bobyqa"))

anova(m0,m1)

# ------------------------------------------------------------------

# as above, but with session baseline

g50_pwr <- sqrt(c(dfb2$preCue.g50,dfb2$postCue.g50))
g80_pwr <- sqrt(c(dfb2$preCue.g80,dfb2$postCue.g80))

pp <- c(rep(0,1933),rep(1,1933))
pp <- factor(pp) # before/after cue, main variable of interest

subj <- c(df$subject,df$subject)
subj <- factor(subj)

cueTime <- c(df$cue.time,df$cue.time)
preSpd <- c(df$preCue.spd,df$preCue.spd)
postSpd <- c(df$postCue.spd,df$postCue.spd)
xcoord <- c(df$xcoord,df$xcoord)
prevRew <-c(df$prevRew,df$prevRew)

cueGamma_b1 <- data.frame(pp,subj,g50_pwr,g80_pwr,cueTime,preSpd,postSpd,xcoord,prevRew)

# inspect
summary(dfb2$preCue.g50)
summary(dfb2$postCue.g50)

summary(dfb2$preCue.g80)
summary(dfb2$postCue.g80)

boxplot(dfb1$preCue.g50,dfb1$postCue.g50)
boxplot(dfb1$preCue.g80,dfb1$postCue.g80)

# tests (replicates matlab, although not with same p-values..)
(wilcox.test(dfb2$preCue.g50,dfb2$postCue.g50,paired=TRUE))
(wilcox.test(dfb2$preCue.g80,dfb2$postCue.g80,paired=TRUE))

# is the mixed model better than the one without?
ml <- lm(g50_pwr ~ 1 + pp,data=cueGamma_b1)
m0 <- lmer(g50_pwr ~ (1|subj),data=cueGamma_b1,REML=FALSE)
m1 <- lmer(g50_pwr ~ (1|subj) + pp,data=cueGamma_b1,REML=FALSE)

anova(m1,ml)

# is there an effect of cue?
anova(m0,m1)

# coefficients
m1

# try some more complicated models
m0 <- lmer(g50_pwr ~ (1|subj) + cueTime,data=cueGamma_b1,REML=FALSE,control=lmerControl(optimizer="bobyqa"))
m1 <- lmer(g50_pwr ~ (1|subj) + cueTime + pp,data=cueGamma_b1,REML=FALSE,control=lmerControl(optimizer="bobyqa"))

anova(m0,m1)

# --------------- gamma-80 -------------------------------------
ml <- lm(g80_pwr ~ 1 + pp,data=cueGamma_b1)
m0 <- lmer(g80_pwr ~ (1|subj),data=cueGamma_b1,REML=FALSE)
m1 <- lmer(g80_pwr ~ (1|subj) + pp,data=cueGamma_b1,REML=FALSE)

anova(m1,ml)

anova(m0,m1)

m1

# try some more complicated models
m0 <- lmer(g80_pwr ~ (1|subj) + cueTime,data=cueGamma_b1,REML=FALSE,control=lmerControl(optimizer="bobyqa"))
m1 <- lmer(g80_pwr ~ (1|subj) + cueTime + pp,data=cueGamma_b1,REML=FALSE,control=lmerControl(optimizer="bobyqa"))

anova(m0,m1)
