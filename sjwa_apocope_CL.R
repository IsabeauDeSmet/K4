#Sjwa apocope

#=============================================
#---------------------------------------------
# 1. SETTING WORKING DIRECTORY & LOADING FILES
# --------------------------------------------
#=============================================
setwd("C:/Users/u0112465/OneDrive - KU Leuven/Onderzoek/K4/sjwa_apocope")

d <- read.csv("Dataset_sjwa_apocope_token_CL.csv")

library(dplyr)
library(reshape2)
library(lme4)
library(ModelMetrics)
library(MuMIn)

#=============================================
#---------------------------------------------
# 2. PREPPING TOKEN DATASET FOR ANALYSIS
# --------------------------------------------
#=============================================

#Random factor author -> less levels (too many levels with only very little attestations)
filter.infrequent <- function(words, threshold = 5, dummy = "OTHER") {
  # code from WBRS for recoding infrequent factor levels (default is <= 5
  # observations)
  return (relevel(
    as.factor(
      ifelse(words %in% levels(as.factor(words))[table(words) >= threshold],
             as.character(words), dummy)), dummy))
}
d$Author_filtered <- filter.infrequent(words = d$Author, threshold=5, dummy="OTHER")

d$Century <- as.numeric(d$Century)
d$cCentury <- scale(d$Century, center=TRUE, scale=TRUE)
d$clogFreq_meaning <- scale(d$logFreq_meaning, center=TRUE, scale=FALSE)
d$clogFreq_lemma <- scale(d$logFreq_lemma, center=TRUE, scale=FALSE)
d$clogFreq_lemma_century <- scale(d$logFreq_lemma_century, center=TRUE, scale=FALSE)
d$clag_variation <- scale(d$lag_variation, center=TRUE, scale=TRUE)

#=============================================
#---------------------------------------------
# 3. ANALYSIS TOKEN DATASET
# --------------------------------------------
#=============================================

#VIF-scores (Zuur et al. 2009)
corvif <- function(dataz) {
  dataz <- as.data.frame(dataz)
  
  #vif part
  form    <- formula(paste("fooy ~ ",paste(strsplit(names(dataz)," "),collapse=" + ")))
  dataz   <- data.frame(fooy=1 + rnorm(nrow(dataz)) ,dataz)
  lm_mod  <- lm(form,dataz)
  
  cat("\n\nVariance inflation factors\n\n")
  print(myvif(lm_mod))
}
#Support function for corvif. Will not be called by the user
myvif <- function(mod) {
  v <- vcov(mod)
  assign <- attributes(model.matrix(mod))$assign
  if (names(coefficients(mod)[1]) == "(Intercept)") {
    v <- v[-1, -1]
    assign <- assign[-1]
  } else warning("No intercept: vifs may not be sensible.")
  terms <- labels(terms(mod))
  n.terms <- length(terms)
  if (n.terms < 2) stop("The model contains fewer than 2 terms")
  if (length(assign) > dim(v)[1] ) {
    diag(tmp_cor)<-0
    if (any(tmp_cor==1.0)){
      return("Sample size is too small, 100% collinearity is present")
    } else {
      return("Sample size is too small")
    }
  }
  R <- cov2cor(v)
  detR <- det(R)
  result <- matrix(0, n.terms, 3)
  rownames(result) <- terms
  colnames(result) <- c("GVIF", "Df", "GVIF^(1/2Df)")
  for (term in 1:n.terms) {
    subs <- which(assign == term)
    result[term, 1] <- det(as.matrix(R[subs, subs])) * det(as.matrix(R[-subs, -subs])) / detR
    result[term, 2] <- length(subs)
  }
  if (all(result[, 2] == 1)) {
    result <- data.frame(GVIF=result[, 1])
  } else {
    result[, 3] <- result[, 1]^(1/(2 * result[, 2]))
  }
  invisible(result)
}



#Only one semantic factor per model


##Metaphor##
vars <- cbind(d$Variant, d$Metaphor, 
              d$cCentury, d$clogFreq_lemma, 
              d$clag_variation, d$Register, d$Lemma, d$Author_filtered, d$Rhyme)

corvif(vars)

fitA <- glmer(Variant ~ Metaphor + cCentury + clag_variation + clogFreq_lemma +
                Register + Poezie + (1+Metaphor|Lemma) + (1+cCentury|Lemma)  +
                (1+Register|Lemma) + (1+Rhyme|Lemma) + (1+Metaphor|Author_filtered),
              data=d, family=binomial)
#does not converge
#zero-correlation model
#transform categorical variables for zero correlation random slope
d$Metaphor_rs <- model.matrix(fitA)[,2]
d$Register_rs <- model.matrix(fitA)[,6]
d$Rhyme_rs <- model.matrix(fitA)[,7]

fitB <- glmer(Variant ~ Metaphor + cCentury + clag_variation + clogFreq_lemma +
                Register + Poezie + (0+Metaphor_rs|Lemma) + (0+cCentury|Lemma)  +
                (0+Register_rs|Lemma) + (0+Rhyme_rs|Lemma) + (0+Metafoor_rs|Author_filtered) +
                (1|Author_filtered) + (1|Lemma),
              data=d, family=binomial)
#does not converge
#drop random slope register by lemma
fitC <- glmer(Variant ~ Metaphor + cCentury + clag_variation + clogFreq_lemma +
                Register + Rhyme + (0+Metaphor_rs|Lemma) + (0+cCentury|Lemma) + (0+Rhyme_rs|Lemma) + (0+Metaphor_rs|Author_filtered) +
                (1|Author_filtered) + (1|Lemma),
              data=d, family=binomial)
anova(fitB, fitC)
#fitC better, yet doesn't converge
#drop random slope rhyme by lemma
fitD <- glmer(Variant ~ Metaphor + cCentury + clag_variation + clogFreq_lemma +
                Register + Poezie + (0+Metaphor_rs|Lemma) + (0+cCentury|Lemma) + (0+Metaphor_rs|Author_filtered) +
                (1|Author_filtered) + (1|Lemma),
              data=d, family=binomial)
anova(fitC, fitD)
#fitD better, yet doesn't converge
#drop random slope metaphor by author
fitE <- glmer(Variant ~ Metaphor + cCentury + clag_variation + clogFrequentie_lemma +
                Register + Rhyme + (0+Metaphor_rs|Lemma) + (0+cCentury|Lemma)  +
                (1|Author_filtered) + (1|Lemma),
              data=d, family=binomial)
anova(fitD, fitE)
#fitE better, yet doesn't converge
#drop random slope metaphor by lemma
fitF <- glmer(Variant ~ Metafoor + cCentury + clag_variation + clogFrequentie_lemma +
                Register + Poezie +  (0+cCentury|Lemma)  +
                (1|Author_filtered) + (1|Lemma),
              data=d, family=binomial)
anova(fitE, fitF)
#fitF not significantly worse, yet doesn't converge
#drop random slope century by lemma
fitG <- glmer(Variant ~ Metaphor + cCentury + clag_variation + clogFreq_lemma +
                Register + Rhyme  +
                (1|Author_filtered) + (1|Lemma),
              data=d, family=binomial)
anova(fitF, fitG)
#fitG worse
#drop random intercept author
fitH <- glmer(Variant ~ Metaphor + cCentury + clag_variation + clogFreq_lemma +
                Register + Rhyme +  (0+cCentury|Lemma) + (1|Lemma),
              data=d, family=binomial)
anova(fitF, fitH)
#fitF better
#add correlation again
fitI <- glmer(Variant ~ Metaphor + cCentury + clag_variation + clogFreq_lemma +
                Register + Rhyme +  (1+cCentury|Lemma)  +
                (1|Author_filtered),
              data=d, family=binomial)
anova(fitF, fitI)
#fitF better
#no convergence
#add bobyqa
fitJ <- glmer(Variant ~ Metaphor + cCentury + clag_variation + clogFreq_lemma +
                Register + Rhyme +  (0+cCentury|Lemma)  +
                (1|Author_filtered) + (1|Lemma),
              data=d, family=binomial, control=glmerControl(optimizer = "bobyqa"))
#converges
#=> no effect metaphor

#add region
fitK <- glmer(Variant ~ Metaphor + cCentury + clag_variation + clogFreq_lemma + Region +
                Register + Rhyme +  (0+cCentury|Lemma)  +
                (1|Author_filtered) + (1|Lemma),
              data=d, family=binomial, control=glmerControl(optimizer = "bobyqa"))

r.squaredGLMM(fitJ)
Preds <- predict(fitJ, type="response")
auc(d$Variant, Preds) 


##Original meaning##
vars <- cbind(d$Variant, d$Original_meaning, 
              d$cCentury, d$clogFreq_lemma, d$clag_variation, d$Register, d$Lemma, d$Author_filtered, d$Rhyme)

corvif(vars)

#full model
fit_A1 <- glmer(Variant ~ Original_meaning + cCentury + clag_variation + clogFreq_lemma +
              Register + Rhyme + (1+Original_meaning|Lemma) + (1+cCentury|Lemma)  +
              (1+Register|Lemma) + (1+Rhyme|Lemma) + (1+Original_meaning|Author_filtered),
            data=d, family=binomial)
#doesn't converge
#no correlation model
#transform original meaning:
d$Original_meaning_rs <- model.matrix(fit_A1)[,2]
fit_A <- glmer(Variant ~ Original_meaning + cCentury + clag_variation + clogFreq_lemma +
                 Register + Rhyme + (0+Original_meaning_rs|Lemma) + (0+cCentury|Lemma)  +
                 (0+Register_rs|Lemma) + (0+Rhyme_rs|Lemma) + (0+Original_meaning_rs|Author_filtered) +
                 (1|Author_filtered) + (1|Lemma),
               data=d, family=binomial)
#doesn't converge
#drop rs rhyme by lemma
fit_B <- glmer(Variant ~ Original_meaning + cCentury + clag_variation + clogFreq_lemma +
                 Register + Rhyme + (0+Original_meaning_rs|Lemma) + (0+cCentury|Lemma)  +
                 (0+Register_rs|Lemma) + (0+Original_meaning_rs|Author_filtered) +
                 (1|Author_filtered) + (1|Lemma),
               data=d, family=binomial)
anova(fit_A, fit_B)
#fit_B better, no convergence
#drop rs register by lemma
fit_C <- glmer(Variant ~ Original_meaning_rs + cCentury + clag_variation + clogFreq_lemma +
                 Register + Rhyme + (0+Original_meanings_rs|Lemma) + (0+cCentury|Lemma)  +
                 (0+Original_meaning_rs|Author_filtered) +
                 (1|Author_filtered) + (1|Lemma),
               data=d, family=binomial)
anova(fit_B, fit_C)
#fit_C better, no convergence
#drop rs original_meaning by lemma
fit_D <- glmer(Variant ~ Original_meaning_rs + cCentury + clag_variation + clogFreq_lemma +
                 Register + Rhyme +  (0+cCentury|Lemma)  +
                 (0+Original_meaning_rs|Author_filtered) +
                 (1|Author_filtered) + (1|Lemma),
               data=d, family=binomial)
anova(fit_C, fit_D)
#fit_C better
#drop rs original_meaning by author
fit_E <- glmer(Variant ~ Original_meaning_rs + cCentury + clag_variation + clogFreq_lemma +
                 Register + Rhyme + (0+Original_meaning_rs|Lemma) + (0+cCentury|Lemma)  +
                 (1|Author_filtered) + (1|Lemma),
               data=d, family=binomial)
anova(fit_C, fit_E)
#fit_C beter
#add correlation parameters again, one by one
#correlation rs original meaning by lemma
fit_F <- glmer(Variant ~ Original_meaning + cCentury + clag_variation + clogFreq_lemma +
                 Register + Rhyme + (1+Original_meaning_rs|Lemma) + (0+cCentury|Lemma)  +
                 (0+Original_meaning_rs|Author_filtered) +
                 (1|Author_filtered),
               data=d, family=binomial)
anova(fit_C, fit_F)
#fit_C better
#correlation parameter century by lemma
fit_H <- glmer(Variant ~ Original_meaning + cCentury + clag_variation + clogFreq_lemma +
                 Register + Rhyme + (0+Original_meaning_rs|Lemma) + (1+cCentury|Lemma)  +
                 (0+Original_meaning_rs|Author_filtered) +
                 (1|Author_filtered),
               data=d, family=binomial)
anova(fit_C, fit_H)
#fit_H better
#correlation parameter original_meaning by author
fit_I <- glmer(Variant ~ Original_meaning + cCentury + clag_variation + clogFreq_lemma +
                 Register + Rhyme + (0+Original_meaning_rs|Lemma) + (1+cCentury|Lemma)  +
                 (1+Original_meaning_rs|Author_filtered),
               data=d, family=binomial)
anova(fit_H, fit_I)
#fit_I better
#still no convergence
#add bobyqa
fit_K <- glmer(Variant ~ Original_meaning + cCentury + clag_variation + clogFreq_lemma +
                 Register + Rhyme + (0+Original_meaning_rs|Lemma) + (1+cCentury|Lemma)  +
                 (1+Original_meaning_rs|Author_filtered),
               data=d, family=binomial, control=glmerControl(optimizer = "bobyqa"))
#still no convergence
#drop both correlated random slopes, even though they make the model significantly better to allow for convergence
fit_L <- glmer(Variant ~ Original_meaning + cCentury + clag_variation + clogFreq_lemma +
                 Register + Rhyme + (0+Original_meaning_rs|Lemma) + (0+cCentury|Lemma)  +
                 (0+Original_meaning_rs|Author_filtered) + (1|Lemma) +
                 (1|Author_filtered),
               data=d, family=binomial, control=glmerControl(optimizer = "bobyqa"))
#model converges, no large differences due to dropping correlation parameters
#-> no effect original_meaning

#add region
fit_M <- glmer(Variant ~ Original_meaning + cCentury + clag_variation + clogFreq_lemma + Region +
                 Register + Rhyme + (0+Original_meaning_rs|Lemma) + (0+cCentury|Lemma)  +
                 (0+Original_meaning_rs|Author_filtered) + (1|Lemma) +
                 (1|Author_filtered),
               data=d, family=binomial, control=glmerControl(optimizer = "bobyqa"))

#R-squared and C-value
r.squaredGLMM(fit_L)
Preds <- predict(fit_L, type="response")
auc(d$Variant, Preds) 


##Meaning frequency##
vars <- cbind(d$Variant,
              d$cCentury, d$clogFreq_lemma, 
              d$Freq_meaning_cat, d$clag_variation, d$Register, d$Lemma, d$Author_filtered, d$Rhyme)

corvif(vars)
#full model
model0 <- glmer(Variant ~ Freq_meaning_cat + cCentury + clag_variation + clogFreq_lemma +
                  Register + Rhyme + (1+Freq_meaning_cat_rs|Lemma) + (1+cCentury|Lemma)  +
                  (1+Register_rs|Lemma) + (1+Rhyme_rs|Lemma) + (1+Freq_meaning_cat_rs|Author_filtered) +
                  (1|Author_filtered) + (1|Lemma),
                data=d, family=binomial)
#doesn't converge
#no correlation model
#transform freq meaning:
d$Freq_meaning_cat_rs <- model.matrix(model0)[,2]
model1 <- glmer(Variant ~ Freq_meaning_cat + cCentury + clag_variation + clogFreq_lemma +
                  Register + Rhyme + (0+Freq_meaning_cat_rs|Lemma) + (0+cCentury|Lemma)  +
                  (0+Register_rs|Lemma) + (0+Rhyme_rs|Lemma) + (0+Freq_meaning_cat_rs|Author_filtered) +
                  (1|Author_filtered) + (1|Lemma),
                data=d, family=binomial)
#doesn't converge
#drop rhyme by lemma
model2 <- glmer(Variant ~ Freq_meaning_cat + cCentury + clag_variation + clogFreq_lemma +
                  Register + Rhyme + (0+Freq_meaning_cat_rs|Lemma) + (0+cCentury|Lemma)  +
                  (0+Register_rs|Lemma) + (0+Freq_meaning_cat_rs|Author_filtered) +
                  (1|Author_filtered) + (1|Lemma),
                data=d, family=binomial)
anova(model1, model2)
#model2 better, no convergence
#drop rs register by lemma
model3 <- glmer(Variant ~ Freq_meaning_cat + cCentury + clag_variation + clogFreq_lemma +
                  Register + Rhyme + (0+Freq_meaning_cat_rs|Lemma) + (0+cCentury|Lemma) + (0+Freq_meaning_cat_rs|Author_filtered) +
                  (1|Author_filtered) + (1|Lemma),
                data=d, family=binomial)
anova(model2, model3)
#model3 better, no convergence
#frequentie|lemma
model4 <- glmer(Variant ~ Freq_meaning_cat + cCentury + clag_variation + clogFreq_lemma +
                  Register + Rhyme + (0+cCentury|Lemma) + (0+Freq_meaning_cat_rs|Author_filtered) +
                  (1|Author_filtered) + (1|Lemma),
                data=d, family=binomial)
anova(model3, model4)
#model3 better, no convergence
#drop freq by author
model5 <- glmer(Variant ~ Freq_meaning_cat + cCentury + clag_variation + clogFreq_lemma +
                  Register + Rhyme + (0+Freq_meaning_cat_rs|Lemma) + (0+cCentury|Lemma)  +
                  (1|Author_filtered) + (1|Lemma),
                data=d, family=binomial)
anova(model3, model5)
#model5 is not significantly worse, still no convergence
#add correlation parameters again
#correlation freq_meaning by lemma
model6 <- glmer(Variant ~ Freq_meaning_cat + cCentury + clag_variation + clogFreq_lemma +
                  Register + Poezie + (1+Freq_meaning_cat_rs|Lemma) + (0+cCentury|Lemma)  +
                  (1|Author_filtered),
                data=d, family=binomial)
anova(model5, model6)
#model6 better
model7 <- glmer(Variant ~ Freq_meaning_cat + cCentury + clag_variation + clogFreq_lemma +
                  Register + Rhyme + (1+Freq_meaning_cat_rs|Lemma) + (1+cCentury|Lemma)  +
                  (1|Author_filtered),
                data=d, family=binomial)
anova(model6, model7)
#model7 worse
#add bobyqa
model8 <- glmer(Variant ~ Freq_meaning_cat + cCentury + clag_variation + clogFreq_lemma +
                  Register + Rhyme + (1+Freq_meaning_cat_rs|Lemma) + (0+cCentury|Lemma)  +
                  (1|Author_filtered),
                data=d, family=binomial, control=glmerControl(optimizer = "bobyqa"))

#=>no effect frequency meaning

#add region
model9 <- glmer(Variant ~ Freq_meaning_cat + cCentury + clag_variation + clogFreq_lemma + Region +
                  Register + Rhyme + (1+Freq_meaning_cat_rs|Lemma) + (0+cCentury|Lemma)  +
                  (1|Author_filtered),
                data=d, family=binomial, control=glmerControl(optimizer = "bobyqa"))

#R-squared and C-value
r.squaredGLMM(model8)
Preds <- predict(model8, type="response")
auc(d$Variant, Preds) 

##Concreteness##
vars <- cbind(d$Variant, d$Concreteness, 
              d$cCentury, d$clogFreq_lemma, d$clag_variation, d$Register, d$Lemma, d$Author_filtered, d$Rhyme)

corvif(vars)
#full model
modela0 <- glmer(Variant ~ Concreteness + cCentury + clag_variation + clogFreq_lemma +
                  Register + Rhyme + (1+Concreteness|Lemma) + (1+cCentury|Lemma)  +
                  (1+Register|Lemma) + (1+Rhyme|Lemma) + (1+Concreteness|Author_filtered) +
                  (1|Author_filtered) + (1|Lemma),
                data=d, family=binomial)
#doesn't converge
#no correlation model
#transform freq meaning:
d$Concreteness_rs <- model.matrix(modela0)[,2]
modela <- glmer(Variant ~ Concreteness + cCentury + clag_variation + clogFreq_lemma +
                  Register + Rhyme + (0+Concreteness_rs|Lemma) + (0+cCentury|Lemma)  +
                  (0+Register_rs|Lemma) + (0+Rhyme_rs|Lemma) + (0+Concreteness_rs|Author_filtered) +
                  (1|Author_filtered) + (1|Lemma),
                data=d, family=binomial)
#drop register by lemma
modelb <- glmer(Variant ~ Concreteness + cCentury + clag_variation + clogFreq_lemma +
                  Register + Rhyme + (0+Concreteness_rs|Lemma) + (0+cCentury|Lemma)  + (0+Rhyme_rs|Lemma) + (0+Concreteness_rs|Author_filtered) +
                  (1|Author_filtered) + (1|Lemma),
                data=d, family=binomial)
anova(modela, modelb)
#modelb better, no convergence
#drop rhyme by lemma
modelc <- glmer(Variant ~ Concreteness + cCentury + clag_variation + clogFreq_lemma +
                  Register + Rhyme + (0+Concreteness_rs|Lemma) + (0+cCentury|Lemma) + (0+Concreteness_rs|Author_filtered) +
                  (1|Author_filtered) + (1|Lemma),
                data=d, family=binomial)
anova(modelb, modelc)
#modelc better, no convergence
#drop concreteness by author
modeld <- glmer(Variant ~ Concreteness + cCentury + clag_variation + clogFreq_lemma +
                  Register + Rhyme + (0+Concreteness_rs|Lemma) + (0+cCentury|Lemma)+
                  (1|Author_filtered) + (1|Lemma),
                data=d, family=binomial)
anova(modelc, modeld)
#modelc better, no convergence
#drop concreteness by lemma
modele <- glmer(Variant ~ Concreteness + cCentury + clag_variation + clogFreq_lemma +
                  Register + Rhyme  + (0+cCentury|Lemma) + (0+Concreteness_rs|Author_filtered) +
                  (1|Author_filtered) + (1|Lemma),
                data=d, family=binomial)
anova(modelc, modele)
#modelc better
#drop century by lemma
modelf <- glmer(Variant ~ Concreteness + cCentury + clag_variation + clogFreq_lemma +
                  Register + Rhyme + (0+Concreteness_rs|Lemma) + (0+Concreteness_rs|Author_filtered) +
                  (1|Author_filtered) + (1|Lemma),
                data=d, family=binomial)
anova(modelc, modelf)
#modelc better
#add correlated random slopes again
#correlated rs concreteness by lemma
modelg <- glmer(Variant ~ Concreteness + cCentury + clag_variation + clogFreq_lemma +
                  Register + Rhyme + (1+Concreteness_rs|Lemma) + (0+cCentury|Lemma) + (0+Concreteness_rs|Author_filtered) +
                  (1|Author_filtered),
                data=d, family=binomial)
anova(modelc, modelg)
#modelc better
#correlated rs century by lemma
modelh <- glmer(Variant ~ Concreteness + cCentury + clag_variation + clogFreq_lemma +
                  Register + Rhyme + (0+Concreteness_rs|Lemma) + (1+cCentury|Lemma) + (0+Concreteness_rs|Author_filtered) +
                  (1|Author_filtered),
                data=d, family=binomial)
anova(modelc, modelh)
#modelc better
#correlated rs concreteness by author
modeli <- glmer(Variant ~ Concreteness + cCentury + clag_variation + clogFreq_lemma +
                  Register + Rhyme + (0+Concreteness_rs|Lemma) + (0+cCentury|Lemma) + (1+Concreteness_rs|Author_filtered)+ (1|Lemma),
                data=d, family=binomial)
anova(modeli, modelc)
#modelc better
#no convergence
#add bobyqa
modelj <- glmer(Variant ~ Concreteness + cCentury + clag_variation + clogFreq_lemma +
                  Register + Rhyme + (0+Concreteness_rs|Lemma) + (0+cCentury|Lemma) + (0+Concreteness_rs|Author_filtered) +
                  (1|Author_filtered) + (1|Lemma),
                data=d, family=binomial, control=glmerControl(optimizer = "bobyqa"))
#=> no effect concreteness
#add region
modelj <- glmer(Variant ~ Concreteness + cCentury + clag_variation + clogFreq_lemma + Region +
                  Register + Rhyme + (0+Concreteness_rs|Lemma) + (0+cCentury|Lemma) + (0+Concreteness_rs|Author_filtered) +
                  (1|Author_filtered) + (1|Lemma),
                data=d, family=binomial, control=glmerControl(optimizer = "bobyqa"))

#R-squared and C-value
r.squaredGLMM(modelj)
Preds <- predict(modelj, type="response")
auc(d$Variant, Preds) 

#=============================================
#---------------------------------------------
# 4. PREPPING TYPE DATASET FOR ANALYSIS
# --------------------------------------------
#=============================================
setwd("C:/Users/u0112465/OneDrive - KU Leuven/Onderzoek/K4/sjwa_apocope")

df <- read.csv("Dataset_sjwa_apocope_type_CL.csv")


# Reduce levels for random factor author, omdat deze heel erg 'skewed' is, because it's skewed
filter.infrequent <- function(words, threshold = 5, dummy = "OTHER") {
  # code from WBRS for recoding infrequent factor levels (default is <= 5
  # observations)
  return (relevel(
    as.factor(
      ifelse(words %in% levels(as.factor(words))[table(words) >= threshold],
             as.character(words), dummy)), dummy))
}
df$Author_filtered <- filter.infrequent(words = df$Author, threshold=5, dummy="OTHER")

# Reduce levels voor random factor author, because it's skewed
filter.infrequent <- function(words, threshold = 5, dummy = "OTHER") {
  # code from WBRS for recoding infrequent factor levels (default is <= 5
  # observations)
  return (relevel(
    as.factor(
      ifelse(words %in% levels(as.factor(words))[table(words) >= threshold],
             as.character(words), dummy)), dummy))
}
df$Source_filtered <- filter.infrequent(words = df$Source, threshold=5, dummy="OTHER")


df$Century <- as.numeric(df$Century)
df$cCentury <- scale(df$Century, center=TRUE, scale=TRUE)
df$logFrequency <- log10(df$Frequency + 1)
df$clogFrequency <- scale(df$logFrequency, center=TRUE, scale=FALSE)
df$cAverage_AOA <- scale(df$Average_AOA, center=TRUE, scale=TRUE)
df$cConcreteness <- scale(df$Concreteness, center=TRUE, scale=TRUE)
df$cLength <- scale(df$Length, center=TRUE, scale=TRUE)
df$cLag_variation <- scale(df$lag_variation, center=TRUE, scale=TRUE)
df$cPolysemy <- scale(df$Polysemy, center=TRUE, scale=TRUE)

#=============================================
#---------------------------------------------
# 4. ANALYSIS TYPE LEVEL
# --------------------------------------------
#=============================================




##Age of acquisition##
vars <- cbind(df$Variant, df$cCentury, df$clogFrequency, df$cAverage_AOA,
              df$clag_variation, df$Register,
              df$Lemma, df$Author_filtered)

corvif(vars)
#full model
m0 <- glmer(Variant ~ cCentury + clogFrequency + cAverage_AOA  +  cLag_variation 
            + Register + Rhyme + (1+cCentury|Lemma) 
            + (1+Register|Lemma) + (1+Rhyme|Lemma) + (1+cAverage_AOA|Author_filtered)
            , data=df, family=binomial)
#no convergence
#model without correlation parameters
#transform register & rhyme
df$Register_rs <- model.matrix(m0)[,6]
df$Rhyme_rs <- model.matrix(m0)[,7]
#no correlation model
m1 <- glmer(Variant ~ cCentury + clogFrequency + cAverage_AOA  +  cLag_variation 
            + Register + Rhyme + (0+cCentury|Lemma) 
            + (0+Register_rs|Lemma) + (0+Rhyme_rs|Lemma) + (0+cAverage_AOA|Author_filtered)  +
              (1|Lemma) + (1|Author_filtered), data=df, family=binomial)
#no convergence
#drop rhyme by lemma
m2 <- glmer(Variant ~ cCentury + clogFrequency + cAverage_AOA  +  cLag_variation 
            + Register + Rhyme + (0+cCentury|Lemma) 
            + (0+Register_rs|Lemma) + (0+cAverage_AOA|Author_filtered)  +
              (1|Lemma) + (1|Author_filtered), data=df, family=binomial)
anova(m1, m2)
#m1 better
#drop AoA by Author
m3 <- glmer(Variant ~ cCentury + clogFrequency + cAverage_AOA  +  cLag_variation 
            + Register + Rhyme + (0+cCentury|Lemma) 
            + (0+Register_rs|Lemma) + (0+Rhyme_rs|Lemma)   +
              (1|Lemma) + (1|Author_filtered), data=df, family=binomial)
anova(m1, m3)
#m1 better
#drop register by lemma
m4 <- glmer(Variant ~ cCentury + clogFrequency + cAverage_AOA  +  cLag_variation 
            + Register + Rhyme + (0+cCentury|Lemma) 
            + (0+Rhyme_rs|Lemma) + (0+cAverage_AOA|Author_filtered)  +
              (1|Lemma) + (1|Author_filtered), data=df, family=binomial)
anova(m1, m4)
#m1 better
#add correlation parameters again
#correlation for century by lemma
m5 <- glmer(Variant ~ cCentury + clogFrequency + cAverage_AOA  +  cLag_variation 
            + Register + Rhyme + (1+cCentury|Lemma) 
            + (0+Register_rs|Lemma) + (0+Rhyme_rs|Lemma) + (0+cAverage_AOA|Author_filtered)  +
              (1|Author_filtered), data=df, family=binomial)
anova(m1, m5)
#m5 better
#correlation for register by lemma
m6 <- glmer(Variant ~ cCentury + clogFrequency + cAverage_AOA  +  cLag_variation 
            + Register + Rhyme + (1+cCentury|Lemma) 
            + (1+Register_rs|Lemma) + (0+Rhyme_rs|Lemma) + (0+cAverage_AOA|Author_filtered)  +
              (1|Author_filtered), data=df, family=binomial)
anova(m5, m6)
#m5 better
#correlation for rhyme by lemma
m7 <- glmer(Variant ~ cCentury + clogFrequency + cAverage_AOA  +  cLag_variation 
            + Register + Rhyme + (1+cCentury|Lemma) 
            + (0+Register_rs|Lemma) + (1+Rhyme_rs|Lemma) + (0+cAverage_AOA|Author_filtered)  +
              (1|Author_filtered), data=df, family=binomial)
anova(m5, m7)
#m5 better
#correlation AoA by author
m8 <- glmer(Variant ~ cCentury + clogFrequency + cAverage_AOA  +  cLag_variation 
            + Register + Rhyme + (1+cCentury|Lemma) 
            + (0+Register_rs|Lemma) + (0+Rhyme_rs|Lemma) + (1+cAverage_AOA|Author_filtered) , data=df, family=binomial)
anova(m5, m8)
#m5 better
#no convergence
#add bobyqa 
m5a <- glmer(Variant ~ cCentury + clogFrequency + cAverage_AOA  +  cLag_variation 
             + Register + Rhyme + (1+cCentury|Lemma) 
             + (0+Register_rs|Lemma) + (0+Rhyme_rs|Lemma) + (0+cAverage_AOA|Author_filtered)  +
               (1|Author_filtered), data=df, family=binomial,
             control=glmerControl(optimizer = "bobyqa"))
#add region
m9 <- glmer(Variant ~ cCentury + clogFrequency + cAverage_AOA  +  cLag_variation + Region 
            + Register + Rhyme + (1+cCentury|Lemma) 
            + (0+Register_rs|Lemma) + (0+Rhyme_rs|Lemma) + (0+cAverage_AOA|Author_filtered)  +
              (1|Author_filtered), data=df, family=binomial)

#->no effect of AOA
#R-squared and C-value
r.squaredGLMM(m5a)
Preds <- predict(m5a, type="response")
auc(df$Variant, Preds) 



##Concreteness##
vars <- cbind(df$Variant, df$cCentury, df$clogFrequency, df$cConcreteness,
              df$clag_variation, df$Register,
              df$Lemma, df$Author_filtered)

corvif(vars)
#full model
m0 <- glmer(Variant ~ cCentury + clogFrequency + cConcreteness  +  cLag_variation 
            + Register + Rhyme + (1+cCentury|Lemma) 
            + (1+Register|Lemma) + (1+Rhyme|Lemma) + (1+cConcreteness|Author_filtered)
            , data=df, family=binomial)
#no convergence
#no correlation model
m_1 <- glmer(Variant ~ cCentury + clogFrequency + cConcreteness  +  cLag_variation 
             + Register + Rhyme + (0+cCentury|Lemma) 
             + (0+Register_rs|Lemma) + (0+Rhyme_rs|Lemma) + (0+cConcreteness|Author_filtered)  +
               (1|Lemma) + (1|Author_filtered), data=df, family=binomial)
#drop rhyme by lemma
m_2 <- glmer(Variant ~ cCentury + clogFrequency + cConcreteness  +  cLag_variation 
             + Register + Poezie + (0+cCentury|Lemma) 
             + (0+Register_rs|Lemma) + (0+cConcreteness|Author_filtered)  +
               (1|Lemma) + (1|Author_filtered), data=df, family=binomial)
anova(m_1, m_2)
#m_1 better
#drop concreteness by author
m_3 <- glmer(Variant ~ cCentury + clogFrequency + cConcreteness  +  cLag_variation 
             + Register + Rhyme + (0+cCentury|Lemma) 
             + (0+Register_rs|Lemma) + (0+Rhyme_rs|Lemma) + 
               (1|Lemma) + (1|Author_filtered), data=df, family=binomial)
anova(m_1, m_3)
#m_1 better
#drop register by lemma
m_4 <- glmer(Variant ~ cCentury + clogFrequency + cConcreteness  +  cLag_variation 
             + Register + Rhyme + (0+cCentury|Lemma) 
             + (0+Rhyme_rs|Lemma) + (0+cConcreteness|Author_filtered)  +
               (1|Lemma) + (1|Author_filtered), data=df, family=binomial)
anova(m_1, m_4)
#m_1 better 
#add correlation parameters again
#correlation century by lemma
m_5 <- glmer(Variant ~ cCentury + clogFrequency + cConcreteness  +  cLag_variation 
             + Register + Rhyme + (1+cCentury|Lemma) 
             + (0+Register_rs|Lemma) + (0+Rhyme_rs|Lemma) + (0+cConcreteness|Author_filtered)  +
               (1|Author_filtered), data=df, family=binomial)
anova(m_1, m_5)
#m_5 not significantly better
#correlation register by lemma
m_6 <- glmer(Variant ~ cCentury + clogFrequency + cConcreteness  +  cLag_variation 
             + Register + Rhyme + (0+cCentury|Lemma) 
             + (1+Register_rs|Lemma) + (0+Rhyme_rs|Lemma) + (0+cConcreteness|Author_filtered)  +
               (1|Author_filtered), data=df, family=binomial)
anova(m_1, m_6)
#m_1 better
#correlation rhyme by lemma
m_7 <- glmer(Variant ~ cCentury + clogFrequency + cConcreteness  +  cLag_variation 
             + Register + Rhyme + (0+cCentury|Lemma) 
             + (0+Register_rs|Lemma) + (1+Rhyme_rs|Lemma) + (0+cConcreteness|Author_filtered)  +
               (1|Author_filtered), data=df, family=binomial)
anova(m_1, m_7)
#m_1 better
#correlation concreteness by author
m_8 <- glmer(Variant ~ cCentury + clogFrequency + cConcreteness  +  cLag_variation 
             + Register + Rhyme + (0+cCentury|Lemma) 
             + (0+Register_rs|Lemma) + (0+Rhyme_rs|Lemma) + (1+cConcreteness|Author_filtered)  +
               (1|Lemma), data=df, family=binomial)
anova(m_1, m_8)
#m_1 better
#no convergence
#add bobyqa 
m_1a <- glmer(Variant ~ cCentury + clogFrequency + cConcreteness  +  cLag_variation 
              + Register + Rhyme + (0+cCentury|Lemma) 
              + (0+Register_rs|Lemma) + (0+Rhyme_rs|Lemma) + (0+cConcreteness|Author_filtered)  +
                (1|Lemma) + (1|Author_filtered), data=df, family=binomial
              ,control=glmerControl(optimizer = "bobyqa"))
#=> effect of concreteness, yet nog significant
#add region
m_9 <- glmer(Variant ~ cCentury + clogFrequency + cConcreteness  +  cLag_variation 
             + Register + Poezie + Region + (0+cCentury|Lemma) 
             + (0+Register_rs|Lemma) + (0+Rhyme_rs|Lemma) + (0+cConcreteness|Author_filtered)  +
               (1|Lemma) + (1|Author_filtered), data=df, family=binomial)
#R-squared and C-value
r.squaredGLMM(m_1a)
Preds <- predict(m_1a, type="response")
auc(df$Variant, Preds) 

##Polysemy##
vars <- cbind(df$Variant, df$cCentury, df$clogFrequency, df$cPolysemy,
              df$clag_variation, df$Register,
              df$Lemma, df$Author_filtered)

corvif(vars)
#full model
a0 <- glmer(Variant ~ cCentury + clogFrequency + cPolysemy  +  cLag_variation 
            + Register + Rhyme + (1+cCentury|Lemma) 
            + (1+Register|Lemma) + (1+Rhyme|Lemma) + (1+cPolysemy|Author_filtered)
            , data=df, family=binomial)
#no convergence
#no correlation model
a1 <- glmer(Variant ~ cCentury + clogFrequency + cPolysemy  +  cLag_variation 
            + Register + Rhyme + (0+cCentury|Lemma) 
            + (0+Register_rs|Lemma) + (0+Rhyme_rs|Lemma) + (0+cPolysemy|Author_filtered)  +
              (1|Lemma) + (1|Author_filtered), data=df, family=binomial)
#no convergence
#drop polysemy by author
a2 <- glmer(Variant ~ cCentury + clogFrequency + cPolysemy  +  cLag_variation 
            + Register + Rhyme + (0+cCentury|Lemma) 
            + (0+Register_rs|Lemma) + (0+Rhyme_rs|Lemma)   +
              (1|Lemma) + (1|Author_filtered), data=df, family=binomial)
anova(a1, a2)
#a1 better
#drop rhyme by lemma
a3 <- glmer(Variant ~ cCentury + clogFrequency + cPolysemy  +  cLag_variation 
            + Register + Rhyme + (0+cCentury|Lemma) 
            + (0+Register_rs|Lemma) + (0+cPolysemy|Author_filtered)  +
              (1|Lemma) + (1|Author_filtered), data=df, family=binomial)
anova(a1, a3)
#a1 better
#drop register by lemma
a4 <- glmer(Variant ~ cCentury + clogFrequency + cPolysemy  +  cLag_variation 
            + Register + Rhyme + (0+cCentury|Lemma) 
            + (0+Rhyme_rs|Lemma) + (0+cPolysemy|Author_filtered)  +
              (1|Lemma) + (1|Author_filtered), data=df, family=binomial)
anova(a1, a4)
#a1 better
#add correlation again
#correlation century by lemma
a5 <- glmer(Variant ~ cCentury + clogFrequency + cPolysemy  +  cLag_variation 
            + Register + Rhyme + (1+cCentury|Lemma) 
            + (0+Register_rs|Lemma) + (0+Rhyme_rs|Lemma) + (0+cPolysemy|Author_filtered)  +
              (1|Author_filtered), data=df, family=binomial)
anova(a1, a5)
#a5 better
#correlation register by lemma
a6 <- glmer(Variant ~ cCentury + clogFrequency + cPolysemy  +  cLag_variation 
            + Register + Rhyme + (1+cCentury|Lemma) 
            + (1+Register_rs|Lemma) + (0+Rhyme_rs|Lemma) + (0+cPolysemy|Author_filtered)  +
              (1|Author_filtered), data=df, family=binomial)
anova(a5, a6)
#a5 better
#correlation rhyme by lemma
a7 <- glmer(Variant ~ cCentury + clogFrequency + cPolysemy  +  cLag_variation 
            + Register + Rhyme + (1+cCentury|Lemma) 
            + (0+Register_rs|Lemma) + (1+Rhyme_rs|Lemma) + (0+cPolysemy|Author_filtered)  +
              (1|Author_filtered), data=df, family=binomial)
anova(a5, a7)
#a5 better
#no convergence
#add bobyqa
a5a <- glmer(Variant ~ cCentury + clogFrequency + cPolysemy  +  cLag_variation 
             + Register + Rhyme + (1+cCentury|Lemma) 
             + (0+Register_rs|Lemma) + (0+Rhyme_rs|Lemma) + (0+cPolysemy|Author_filtered)  +
               (1|Author_filtered), data=df, family=binomial, control=glmerControl(optimizer = "bobyqa"))
#add region (no effect)
a8 <- glmer(Variant ~ cCentury + clogFrequency + cPolysemy  +  cLag_variation + Regio
            + Register + Rhyme + (1+cCentury|Lemma) 
            + (0+Register_rs|Lemma) + (0+Rhyme_rs|Lemma) + (0+cPolysemy|Author_filtered)  +
              (1|Author_filtered), data=df, family=binomial)
# => no effect polysemy

#R-squared and C-value
r.squaredGLMM(a5a)
Preds <- predict(a5a, type="response")
auc(df$Variant, Preds) 