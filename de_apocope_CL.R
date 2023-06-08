#de apocope#

#=============================================
#---------------------------------------------
# 1. SETTING WORKING DIRECTORY & LOADING FILES
# --------------------------------------------
#=============================================
setwd("C:/Users/u0112465/OneDrive - KU Leuven/Onderzoek/K4/de_apocope")

d <- read.csv("Dataset_de_apocope_token_CL.csv")

library(dplyr)
library(reshape2)
library(lme4)
library(ModelMetrics)
library(MuMIn)
#=============================================
#---------------------------------------------
# 2. CLEANING DATASET
# --------------------------------------------
#=============================================

#center & scale
d$Century <- as.numeric(d$Century)
d$cCentury <- scale(d$Century, center=TRUE, scale=TRUE)
d$clogFreq_meaning <- scale(d$logFreq_meaning, center=TRUE, scale=FALSE)
d$clogFreq_lemma <- scale(d$logFreq_lemma, center=TRUE, scale=FALSE)
d$clag_variation <- scale(d$lag_variation, center=TRUE, scale=TRUE)



# Author variable is skewed -> reduce levels
filter.infrequent <- function(words, threshold = 5, dummy = "OTHER") {
  # code from WBRS for recoding infrequent factor levels (default is <= 5
  # observations)
  return (relevel(
    as.factor(
      ifelse(words %in% levels(as.factor(words))[table(words) >= threshold],
             as.character(words), dummy)), dummy))
}
d$Author_filtered <- filter.infrequent(words = d$Author, threshold=5, dummy="OTHER")

# Source variable is skewed -> reduce levels
filter.infrequent <- function(words, threshold = 5, dummy = "OTHER") {
  # code from WBRS for recoding infrequent factor levels (default is <= 5
  # observations)
  return (relevel(
    as.factor(
      ifelse(words %in% levels(as.factor(words))[table(words) >= threshold],
             as.character(words), dummy)), dummy))
}
d$Source_filtered <- filter.infrequent(words = d$Source, threshold=5, dummy="OTHER")

#clean up outcome factor -> make binary
d$Variant_bin <- d$Variant
levels(d$Variant_bin) <- c("de", "de", "d_syncope", "d_syncope")

#=============================================
#---------------------------------------------
# 3. ANALYSIS TOKEN LEVEL
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




##METAPHOR##
vars <- cbind(d$Variant_bin, d$Metaphor, d$logFreq_lemma, d$cCentury,
              d$clag_variation, d$Register, d$Lemma, d$Author_filtered, d$Rhyme)

corvif(vars)
#full model
fit_0 <- glmer(Variant_bin ~ cCentury  + Metaphor + Rhyme + Register +  clag_variation + clogFreq_lemma +
                 (1+cCentury|Lemma) + (1+Metaphor|Lemma) + (1+Rhyme|Lemma) + (1+Register|Lemma)
               + (1+Metaphor|Author_filtered),
               data=d, family=binomial)
#no convergence
#drop correlation parameters
#transform factors
d$Metaphor_rs <- model.matrix(fit_0)[,3]
d$Rhyme_rs <- model.matrix(fit_0)[,4]
d$Register_rs <- model.matrix(fit_0)[,5]

fit_a <- glmer(Variant_bin ~ cCentury  + Metaphor + Rhyme + Register +  clag_variation + clogFreq_lemma +
                 (0+cCentury|Lemma) + (0+Metaphor_rs|Lemma) + (0+Rhyme_rs|Lemma) + (0+Register_rs|Lemma)
               + (0+Metaphor_rs|Author_filtered)  +
                 (1|Author_filtered) + (1|Lemma),
               data=d, family=binomial)
#no convergence
#drop metaphor by lemma
fit_b <- glmer(Variant_bin ~ cCentury  + Metaphor + Rhyme + Register +  clag_variation + clogFreq_lemma +
                 (0+cCentury|Lemma)  + (0+Rhyme_rs|Lemma) + (0+Register_rs|Lemma)
               + (0+Metaphor_rs|Author_filtered)  +
                 (1|Author_filtered) + (1|Lemma),
               data=d, family=binomial)
anova(fit_a, fit_b)
#fit_b better, no convergence
#drop rhyme by lemma 
fit_c <- glmer(Variant_bin ~ cCentury  + Metaphor + Rhyme + Register +  clag_variation + clogFreq_lemma +
                 (0+cCentury|Lemma)  + (0+Register_rs|Lemma)
               + (0+Metaphor_rs|Author_filtered)  +
                 (1|Author_filtered) + (1|Lemma),
               data=d, family=binomial)
anova(fit_b, fit_c)
#fit_c better, no convergence
#drop metaphor by author
fit_d <- glmer(Variant_bin ~ cCentury  + Metaphor + Rhyme + Register +  clag_variation + clogFreq_lemma +
                  (0+cCentury|Lemma)  + (0+Register_rs|Lemma)
               + (1|Author_filtered) + (1|Lemma),
               data=d, family=binomial)
anova(fit_c, fit_d)
#fit_d better, no convergence
#drop register by lemma 
fit_e <- glmer(Variant_bin ~ cCentury  + Metaphor + Rhyme + Register +  clag_variation + clogFreq_lemma +
                 (0+cCentury|Lemma)  + (1|Author_filtered) + (1|Lemma),
               data=d, family=binomial)
anova(fit_d, fit_e)
#fit_e better, no convergence
#drop century by lemma
fit_f <- glmer(Variant_bin ~ cCentury  + Metaphor + Rhyme + Register +  clag_variation + clogFreq_lemma +
               + (1|Author_filtered) + (1|Lemma),
               data=d, family=binomial)
anova(fit_e, fit_f)
#fit_e better
#drop random effect author
fit_g <- glmer(Variant_bin ~ cCentury  + Metaphor + Rhyme + Register +  clag_variation + clogFreq_lemma +
                 (0+cCentury|Lemma)   + (1|Lemma),
               data=d, family=binomial)
anova(fit_e, fit_g)
#fit_e better
#add correlation parameters again
#correlation century by lemma
fit_h <- glmer(Variant_bin ~ cCentury  + Metaphor + Rhyme + Register +  clag_variation + clogFreq_lemma +
                 (1+cCentury|Lemma)  
               +  (1|Author_filtered) ,
               data=d, family=binomial)
anova(fit_e, fit_h)
#fit_e better
#add bobyqa
fit_i <- glmer(Variant_bin ~ cCentury  + Metaphor + Rhyme + Register +  clag_variation + clogFreq_lemma +
                 (0+cCentury|Lemma)  
               +  (1|Author_filtered) + (1|Lemma),
               data=d, family=binomial, control=glmerControl(optimizer = "bobyqa"))
#still no convergence
#drop century by lemma
fit_j <- glmer(Variant_bin ~ cCentury  + Metaphor + Rhyme + Register +  clag_variation + clogFreq_lemma 
               +  (1|Author_filtered) + (1|Lemma),
               data=d, family=binomial, control=glmerControl(optimizer = "bobyqa"))
#=> small effect metaphor (but not with Bonferroni correction)

#add region
fit_k <- glmer(Variant_bin ~ cCentury  + Metaphor + Rhyme + Register +  clag_variation + Region  
               +  (1|Author_filtered) + (1|Lemma),
               data=d, family=binomial, control=glmerControl(optimizer = "bobyqa"))
#R-squared and C-value
r.squaredGLMM(fit_j)
Preds <- predict(fit_j, type="response")
auc(d$Variant_bin, Preds) 


##Freq meaning##
vars <- cbind(d$Variant_bin, d$Freq_meaning_cat, d$logFreq_lemma, d$cCentury,
              d$clag_variation, d$Register, d$Lemma, d$Author_filtered, d$Rhyme)

corvif(vars)
#full model
fit0 <- glmer(Variant_bin ~ cCentury  + Freq_meaning_cat + Rhyme + Register +  clag_variation +
                 (1+cCentury|Lemma) + (1+Metaphor|Lemma) + (1+Rhyme|Lemma) + (1+Register|Lemma)
               + (1+Freq_meaning_cat|Author_filtered),
               data=d, family=binomial)
#no convergence
#drop correlation parameters
#transform factor
d$Freq_meaning_cat_rs <- model.matrix(fit0)[,3]

fita <- glmer(Variant_bin ~ cCentury  + Freq_meaning_cat + Rhyme + Register +  clag_variation + clogFreq_lemma +
                (0+cCentury|Lemma) + (0+Freq_meaning_cat_rs|Lemma) + (0+Rhyme_rs|Lemma) + (0+Register_rs|Lemma)
              + (0+Freq_meaning_cat_rs|Author_filtered)  +
                (1|Author_filtered) + (1|Lemma),
              data=d, family=binomial)
#no convergence
#drop rhyme by lemma
fitb <- glmer(Variant_bin ~ cCentury  + Freq_meaning_cat + Rhyme + Register +  clag_variation + clogFreq_lemma +
                (0+cCentury|Lemma) + (0+Freq_meaning_cat_rs|Lemma) + (0+Register_rs|Lemma)
              + (0+Freq_meaning_cat_rs|Author_filtered)  +
                (1|Author_filtered) + (1|Lemma),
              data=d, family=binomial)
anova(fita, fitb)
#fitb better, no convergence
#drop register by lemma 
fitc <- glmer(Variant_bin ~ cCentury  + Freq_meaning_cat + Rhyme + Register +  clag_variation + clogFreq_lemma +
                (0+cCentury|Lemma) + (0+Freq_meaning_cat_rs|Lemma) 
              + (0+Freq_meaning_cat_rs|Author_filtered)  +
                (1|Author_filtered) + (1|Lemma),
              data=d, family=binomial)
anova(fitb, fitc)
#fitc better, no convergence
#drop century by lemma
fitd <- glmer(Variant_bin ~ cCentury  + Freq_meaning_cat + Rhyme + Register +  clag_variation + clogFreq_lemma +
                (0+Freq_meaning_cat_rs|Lemma) 
              + (0+Freq_meaning_cat_rs|Author_filtered)  +
                (1|Author_filtered) + (1|Lemma),
              data=d, family=binomial)
anova(fitc, fitd)
#fitc better
#drop freq meaning by lemma
fite <- glmer(Variant_bin ~ cCentury  + Freq_meaning_cat + Rhyme + Register +  clag_variation + clogFreq_lemma +
              + (0+Freq_meaning_cat_rs|Author_filtered)  +
                (1|Author_filtered) + (1|Lemma),
              data=d, family=binomial)
anova(fitc, fite)
#fit c better
#drop freq meaning by author
fitf <- glmer(Variant_bin ~ cCentury  + Freq_meaning_cat + Rhyme + Register +  clag_variation + clogFreq_lemma +
                (0+cCentury|Lemma) + (0+Freq_meaning_cat_rs|Lemma) 
              + (1|Author_filtered) + (1|Lemma),
              data=d, family=binomial)
anova(fitc, fitf)
#fitc better
#add correlated random slopes again
#cor ran slope century by lemma
fitg <- glmer(Variant_bin ~ cCentury  + Freq_meaning_cat + Rhyme + Register +  clag_variation + clogFreq_lemma +
                (1+cCentury|Lemma) + (0+Freq_meaning_cat_rs|Lemma) 
              + (0+Freq_meaning_cat_rs|Author_filtered)  +
                (1|Author_filtered),
              data=d, family=binomial)
anova(fitc, fitg)
#fitg not significantly better
#cor ran slope freq meaning by lemma
fith <- glmer(Variant_bin ~ cCentury  + Freq_meaning_cat + Rhyme + Register +  clag_variation + clogFreq_lemma +
                (0+cCentury|Lemma) + (1+Freq_meaning_cat_rs|Lemma) 
              + (0+Freq_meaning_cat_rs|Author_filtered)  +
                (1|Author_filtered),
              data=d, family=binomial)
anova(fitc, fith)
#fitc better
#cor random slope freq meaning by source
fiti <- glmer(Variant_bin ~ cCentury  + Freq_meaning_cat + Rhyme + Register +  clag_variation + clogFreq_lemma +
                (0+cCentury|Lemma) + (0+Freq_meaning_cat_rs|Lemma) 
              + (1+Freq_meaning_cat_rs|Author_filtered) + (1|Lemma),
              data=d, family=binomial)
anova(fitc, fiti)
#fitc better
#add bobyqa
fitj <- glmer(Variant_bin ~ cCentury  + Freq_meaning_cat + Rhyme + Register +  clag_variation + clogFreq_lemma +
                (0+cCentury|Lemma) + (0+Freq_meaning_cat_rs|Lemma) 
              + (0+Freq_meaning_cat_rs|Author_filtered)  +
                (1|Author_filtered) + (1|Lemma),
              data=d, family=binomial, control=glmerControl(optimizer = "bobyqa"))
#still no convergence
#delete random slope that makes the least difference
fitk <- glmer(Variant_bin ~ cCentury  + Freq_meaning_cat + Rhyme + Register +  clag_variation + clogFreq_lemma +
                (0+Freq_meaning_cat_rs|Lemma) 
              + (0+Freq_meaning_cat_rs|Author_filtered)  +
                (1|Author_filtered) + (1|Lemma),
              data=d, family=binomial, control=glmerControl(optimizer = "bobyqa"))
#still no convergence
#delete second random slope
fitl <- glmer(Variant_bin ~ cCentury  + Freq_meaning_cat + Rhyme + Register +  clag_variation + clogFreq_lemma
              + (0+Freq_meaning_cat_rs|Author_filtered)  +
                (1|Author_filtered) + (1|Lemma),
              data=d, family=binomial, control=glmerControl(optimizer = "bobyqa"))
#delete third random slope
fitm <- glmer(Variant_bin ~ cCentury  + Freq_meaning_cat + Rhyme + Register +  clag_variation + clogFreq_lemma
               +
                (1|Author_filtered) + (1|Lemma),
              data=d, family=binomial, control=glmerControl(optimizer = "bobyqa"))
#=>no effect freq meaning

#add region
fitn <- glmer(Variant_bin ~ cCentury  + Freq_meaning_cat + Rhyme + Register +  clag_variation + clogFreq_lemma
              + Region +
                (1|Author_filtered) + (1|Lemma),
              data=d, family=binomial, control=glmerControl(optimizer = "bobyqa"))

#R-squared and C-value
r.squaredGLMM(fitm)
Preds <- predict(fitm, type="response")
auc(d$Variant_bin, Preds) 

##original meaning##
vars <- cbind(d$Variant_bin, d$Original_meaning, d$logFreq_lemma, d$cCentury,
              d$clag_variation, d$Register, d$Lemma, d$Author_filtered, d$Rhyme)

corvif(vars)


#full model
model0 <- glmer(Variant_bin ~ cCentury  + Original_meaning + Rhyme + Register +  clag_variation +
                (1+cCentury|Lemma) + (1+Original_meaning|Lemma) + (1+Rhyme|Lemma) + (1+Register|Lemma)
              + (1+Original_meaning|Author_filtered),
              data=d, family=binomial)
#no convergence
#drop correlation parameters
#transform factor
d$Original_meaning_rs <- model.matrix(model0)[,3]

model1 <- glmer(Variant_bin ~ cCentury  + Original_meaning + Rhyme + Register +  clag_variation + clogFreq_lemma +
                  (0+cCentury|Lemma) + (0+Original_meaning_rs|Lemma) + (0+Rhyme_rs|Lemma) + (0+Register_rs|Lemma)
                + (0+Original_meaning_rs|Author_filtered)  +
                  (1|Author_filtered) + (1|Lemma),
                data=d, family=binomial)
#no convergence
#drop rhyme by lemma
model2 <- glmer(Variant_bin ~ cCentury  + Original_meaning + Rhyme + Register +  clag_variation + clogFreq_lemma +
                  (0+cCentury|Lemma) + (0+Original_meaning_rs|Lemma)  + (0+Register_rs|Lemma)
                + (0+Original_meaning_rs|Author_filtered)  +
                  (1|Author_filtered) + (1|Lemma),
                data=d, family=binomial)
anova(model1, model2)
#model 2 better, no convergence
#drop register by lemma
model3 <- glmer(Variant_bin ~ cCentury  + Original_meaning + Rhyme + Register +  clag_variation + clogFreq_lemma +
                  (0+cCentury|Lemma) + (0+Original_meaning_rs|Lemma) 
                + (0+Original_meaning_rs|Author_filtered)  +
                  (1|Author_filtered) + (1|Lemma),
                data=d, family=binomial)
anova(model2, model3)
#model3 beter, no convergence
#drop century by lemma
model4 <- glmer(Variant_bin ~ cCentury  + Original_meaning + Rhyme + Register +  clag_variation + clogFreq_lemma +
                  + (0+Original_meaning_rs|Lemma) 
                + (0+Original_meaning_rs|Author_filtered)  +
                  (1|Author_filtered) + (1|Lemma),
                data=d, family=binomial)
anova(model3, model4)
#model4 better
#drop original meaning by lemma
model5 <- glmer(Variant_bin ~ cCentury  + Original_meaning + Rhyme + Register +  clag_variation + clogFreq_lemma +
                + (0+Original_meaning_rs|Author_filtered)  +
                  (1|Author_filtered) + (1|Lemma),
                data=d, family=binomial)
anova(model4, model5)
#model4 better
#drop original meaning by author
model6 <- glmer(Variant_bin ~ cCentury  + Original_meaning + Rhyme + Register +  clag_variation + clogFreq_lemma +
                  (0+Original_meaning_rs|Lemma) 
                + (1|Author_filtered) + (1|Lemma),
                data=d, family=binomial)
anova(model4, model6)
#model4 better
#add cor random slopes again
# correlation original meaning by lemma
model7 <- glmer(Variant_bin ~ cCentury  + Original_meaning + Rhyme + Register +  clag_variation + clogFreq_lemma +
                + (1+Original_meaning_rs|Lemma) 
                + (0+Original_meaning_rs|Author_filtered)  +
                  (1|Author_filtered),
                data=d, family=binomial)
anova(model4, model7)
#model4 better
#correlation original meaning by author
model8 <- glmer(Variant_bin ~ cCentury  + Original_meaning + Rhyme + Register +  clag_variation + clogFreq_lemma +
                  (0+Original_meaning_rs|Lemma) 
                + (1+Original_meaning_rs|Author_filtered)  +
                  (1|Lemma),
                data=d, family=binomial)
anova(model4, model8)
#model8 better
#add bobyqa
model10 <- glmer(Variant_bin ~ cCentury  + Original_meaning + Rhyme + Register +  clag_variation + clogFreq_lemma +
                 + (0+Original_meaning_rs|Lemma) 
                 + (1+Original_meaning_rs|Author_filtered)  +
                   (1|Lemma),
                 data=d, family=binomial, control=glmerControl(optimizer = "bobyqa"))
#still no convergence
#drop correlation again
model11 <- glmer(Variant_bin ~ cCentury  + Original_meaning + Rhyme + Register +  clag_variation + clogFreq_lemma +
                   (0+Original_meaning_rs|Lemma) 
                 + (0+Original_meaning_rs|Author_filtered)  + (1|Author_filtered) +
                   (1|Lemma),
                 data=d, family=binomial, control=glmerControl(optimizer = "bobyqa"))
#still no convergence
#try dropping random slopes
model12 <- glmer(Variant_bin ~ cCentury  + Original_meaning + Rhyme + Register +  clag_variation + clogFreq_lemma +
                   (0+Original_meaning_rs|Author_filtered) 
                 + (1|Author_filtered) +
                   (1|Lemma),
                 data=d, family=binomial, control=glmerControl(optimizer = "bobyqa"))
#no convergence
model13 <- glmer(Variant_bin ~ cCentury  + Original_meaning + Rhyme + Register +  clag_variation  + clogFreq_lemma +
                 +  (1|Author_filtered) +
                   (1|Lemma),
                 data=d, family=binomial, control=glmerControl(optimizer = "bobyqa"))
#no convergence
model14 <- glmer(Variant_bin ~ cCentury  + Original_meaning + Rhyme + Register +  clag_variation  + clogFreq_lemma +
                   +  (1|Author_filtered),
                 data=d, family=binomial, control=glmerControl(optimizer = "bobyqa"))
#converges
#=> no effect original meaning

#add region
model16 <- glmer(Variant_bin ~ cCentury  + Original_meaning + Rhyme + Register +  clag_variation  + Region + clogFreq_lemma +
                   +  (1|Author_filtered),
                 data=d, family=binomial, control=glmerControl(optimizer = "bobyqa"))

#R-squared and C-value
r.squaredGLMM(model14)
Preds <- predict(model14, type="response")
auc(d$Variant_bin, Preds) 

##Concreteness##
vars <- cbind(d$Variant_bin, d$Concreteness, d$logFreq_lemma, d$cCentury,
              d$clag_variation, d$Register, d$Lemma, d$Author_filtered, d$Rhyme)

corvif(vars)

#full model
model_0 <- glmer(Variant_bin ~ cCentury  + Concreteness + Rhyme + Register +  clag_variation +
                  (1+cCentury|Lemma) + (1+Concreteness|Lemma) + (1+Rhyme|Lemma) + (1+Register|Lemma)
                + (1+Concreteness|Author_filtered),
                data=d, family=binomial)
#no convergence
#drop correlation parameters
#transform factor
d$Concreteness_rs <- model.matrix(model_0)[,3]
model_1 <- glmer(Variant_bin ~ cCentury  + Concreteness + Rhyme + Register +  clag_variation + clogFreq_lemma +
                   (0+cCentury|Lemma) + (0+Concreteness_rs|Lemma) + (0+Rhyme_rs|Lemma) + (0+Register_rs|Lemma)
                 + (0+Concreteness_rs|Author_filtered)  +
                   (1|Author_filtered) + (1|Lemma),
                 data=d, family=binomial)
#no convergence
#drop rhyme by lemma
model_2 <- glmer(Variant_bin ~ cCentury  + Concreteness + Rhyme + Register +  clag_variation + clogFreq_lemma +
                   (0+cCentury|Lemma) + (0+Concreteness_rs|Lemma)  + (0+Register_rs|Lemma)
                 + (0+Concreteness_rs|Author_filtered)  +
                   (1|Author_filtered) + (1|Lemma),
                 data=d, family=binomial)
anova(model_1, model_2)
#model_2 better
#drop register by lemma
model_3 <- glmer(Variant_bin ~ cCentury  + Concreteness + Rhyme + Register +  clag_variation + clogFreq_lemma +
                   (0+cCentury|Lemma) + (0+Concreteness_rs|Lemma) 
                 + (0+Concreteness_rs|Author_filtered)  +
                   (1|Author_filtered) + (1|Lemma),
                 data=d, family=binomial)
anova(model_2, model_3)
#model_3 better
#drop concreteness by lemma
model_4 <- glmer(Variant_bin ~ cCentury  + Concreteness + Rhyme + Register +  clag_variation + clogFreq_lemma +
                   (0+cCentury|Lemma) 
                 + (0+Concreteness_rs|Author_filtered)  +
                   (1|Author_filtered) + (1|Lemma),
                 data=d, family=binomial)
anova(model3, model4)
#model_4 better
#drop century by lemma
model_5 <- glmer(Variant_bin ~ cCentury  + Concreteness + Rhyme + Register +  clag_variation + clogFreq_lemma +
                 + (0+Concreteness_rs|Author_filtered)  +
                   (1|Author_filtered) + (1|Lemma),
                 data=d, family=binomial)
anova(model_4, model_5)
#model_4 better
#drop concreteness by author
model_6 <- glmer(Variant_bin ~ cCentury  + Concreteness + Rhyme + Register +  clag_variation + clogFreq_lemma +
                   (0+cCentury|Lemma)   +
                   (1|Author_filtered) + (1|Lemma),
                 data=d, family=binomial)
anova(model_4, model_6)
#model_4 better
#add correlation parameters again
#correlation century by lemma
model_7 <- glmer(Variant_bin ~ cCentury  + Concreteness + Rhyme + Register +  clag_variation + clogFreq_lemma +
                   (1+cCentury|Lemma) 
                 + (0+Concreteness_rs|Author_filtered)  +
                   (1|Author_filtered),
                 data=d, family=binomial)
anova(model_4, model_7)
#model_4 better
#correlation concreteness by author
model_8 <- glmer(Variant_bin ~ cCentury  + Concreteness + Rhyme + Register +  clag_variation + clogFreq_lemma +
                   (0+cCentury|Lemma) 
                 + (1+Concreteness_rs|Author_filtered) + (1|Lemma),
                 data=d, family=binomial)
anova(model_4, model_8)
#model_8 better
#add bobyqa
model_9 <- glmer(Variant_bin ~ cCentury  + Concreteness + Rhyme + Register +  clag_variation + clogFreq_lemma +
                   (0+cCentury|Lemma) 
                 + (1+Concreteness_rs|Author_filtered) + (1|Lemma),
                 data=d, family=binomial, control=glmerControl(optimizer = "bobyqa"))
#still no convergence
#delete correlation again
model_10 <- glmer(Variant_bin ~ cCentury  + Concreteness + Rhyme + Register +  clag_variation + clogFreq_lemma +
                   (0+cCentury|Lemma) 
                 + (0+Concreteness_rs|Author_filtered) + (1|Lemma) + (1|Author_filtered),
                 data=d, family=binomial, control=glmerControl(optimizer = "bobyqa"))
#=> no effect concreteness

#add region
model_11 <- glmer(Variant_bin ~ cCentury  + Region + Concreteness + Rhyme + Register +  clag_variation + clogFreq_lemma +
                   (0+cCentury|Lemma) 
                 + (0+Concreteness_rs|Author_filtered) + (1|Lemma) + (1|Author_filtered),
                 data=d, family=binomial, control=glmerControl(optimizer = "bobyqa"))
#R-squared and C-value
r.squaredGLMM(model_10)
Preds <- predict(model_10, type="response")
auc(d$Variant_bin, Preds) 

#=============================================
#---------------------------------------------
# 4. PREPPING DATASET TYPE LEVEL
# --------------------------------------------
#=============================================
#center and scale
df$Century <- as.numeric(df$Century)
df$cCentury <- scale(df$Century, center=TRUE, scale=TRUE)
df$logFrequency <- log10(df$Frequency + 1)
df$clogFrequency <- scale(df$logFrequency, center=TRUE, scale=FALSE)
df$cAverage_AOA <- scale(df$Average_AOA, center=TRUE, scale=TRUE)
df$cConcreteness <- scale(df$Concreteness, center=TRUE, scale=TRUE)
df$cLag_variation <- scale(df$lag_variation, center=TRUE, scale=TRUE)
df$cPolysemy <- scale(df$Polysemy, center=TRUE, scale=TRUE)

# Fix skew in author
filter.infrequent <- function(words, threshold = 5, dummy = "OTHER") {
  # code from WBRS for recoding infrequent factor levels (default is <= 5
  # observations)
  return (relevel(
    as.factor(
      ifelse(words %in% levels(as.factor(words))[table(words) >= threshold],
             as.character(words), dummy)), dummy))
}
df$Author_filtered <- filter.infrequent(words = df$Author, threshold=5, dummy="OTHER")

# Fix skew in source
filter.infrequent <- function(words, threshold = 5, dummy = "OTHER") {
  # code from WBRS for recoding infrequent factor levels (default is <= 5
  # observations)
  return (relevel(
    as.factor(
      ifelse(words %in% levels(as.factor(words))[table(words) >= threshold],
             as.character(words), dummy)), dummy))
}
df$Source_filtered <- filter.infrequent(words = df$Source, threshold=5, dummy="OTHER")

#=============================================
#---------------------------------------------
# 4. ANALYSIS TYPE LEVEL
# --------------------------------------------
#=============================================

##Age of acquisition##
vars <- cbind(df$Variant_bin, df$clogFrequency, df$cAverage_AOA,
              df$cCentury, 
              df$Register, df$Lemma, df$Source_filtered, df$Rhyme, df$cLag_variation)
corvif(vars)
#source_filtered te hoge vif's -> opgelost met author ipv source

#full model
m1 <- glmer(Variant ~ clogFrequency + cAverage_AOA + cCentury + Register + Rhyme
            + cLag_variation + (1+cCentury|Lemma) + (1+Register|Lemma) + (1+Rhyme|Lemma)
            + (1+cAverage_AOA|Author_filtered), data=df, family=binomial)
#no convergence
#no correlation model
#transform register & rhyme
df$Register_rs <- model.matrix(m1)[,5]
df$Rhyme_rs <- model.matrix(m1)[,6]

m2 <- glmer(Variant_bin ~ clogFrequency + cAverage_AOA + cCentury + Register + Rhyme
            + cLag_variation + (0+cCentury|Lemma) + (0+Register_rs|Lemma) + (0+Rhyme_rs|Lemma)
            + (0+cAverage_AOA|Author_filtered) + (1|Lemma) + (1|Author_filtered),
            family=binomial, data=df)
#no convergence
#drop aoa by author
m3 <- glmer(Variant_bin ~ clogFrequency + cAverage_AOA + cCentury + Register + Rhyme
            + cLag_variation + (0+cCentury|Lemma) + (0+Register_rs|Lemma) + (0+Rhyme_rs|Lemma)
            + (1|Lemma) + (1|Author_filtered),
            family=binomial, data=df)
anova(m2, m3)
#m2 better
#drop century by lemma
m4 <- glmer(Variant_bin ~ clogFrequency + cAverage_AOA + cCentury + Register + Rhyme
            + cLag_variation  + (0+Register_rs|Lemma) + (0+Rhyme_rs|Lemma)
            + (0+cAverage_AOA|Author_filtered) + (1|Lemma) + (1|Author_filtered),
            family=binomial, data=df)
anova(m2, m4)
#m2 better
#drop register by lemma
m5 <- glmer(Variant_bin ~ clogFrequency + cAverage_AOA + cCentury + Register + Rhyme
            + cLag_variation + (0+cCentury|Lemma) + (0+Rhyme_rs|Lemma)
            + (0+cAverage_AOA|Author_filtered) + (1|Lemma) + (1|Author_filtered),
            family=binomial, data=df)
anova(m2, m5)
#m5 not significantly worse, but simpler, yet no convergence
#add correlation random slopes again
#correlation century by lemma
m6 <- glmer(Variant_bin ~ clogFrequency + cAverage_AOA + cCentury + Register + Rhyme
            + cLag_variation + (1+cCentury|Lemma) + (0+Rhyme_rs|Lemma)
            + (0+cAverage_AOA|Author_filtered)+ (1|Author_filtered),
            family=binomial, data=df)
anova(m5, m6)
#m5 better
#correlation rhyme by lemma
m7 <- glmer(Variant_bin ~ clogFrequency + cAverage_AOA + cCentury + Register + Rhyme
            + cLag_variation + (0+cCentury|Lemma) + (1+Rhyme_rs|Lemma)
            + (0+cAverage_AOA|Author_filtered)+ (1|Author_filtered),
            family=binomial, data=df)
anova(m5, m7)
#m5 better
m8 <- glmer(Variant_bin ~ clogFrequency + cAverage_AOA + cCentury + Register + Rhyme
            + cLag_variation + (0+cCentury|Lemma) + (0+Rhyme_rs|Lemma)
            + (1+cAverage_AOA|Author_filtered) + (1|Lemma),
            family=binomial, data=df)
anova(m5, m8)
#m5 better
#add bobyqa
m5a <- glmer(Variant_bin ~ clogFrequency + cAverage_AOA + cCentury + Register + Rhyme
             + cLag_variation + (0+cCentury|Lemma) + (0+Rhyme_rs|Lemma)
             + (0+cAverage_AOA|Author_filtered) + (1|Lemma) + (1|Author_filtered),
             family=binomial, data=df, control=glmerControl(optimizer = "bobyqa"))
#-> effect of age of acquisition, yet not when p-value is adjusted with Bonferroni correction

#add region
m9 <- glmer(Variant_bin ~ clogFrequency + cAverage_AOA + cCentury + Register + Rhyme + Region +
             + cLag_variation + (0+cCentury|Lemma) + (0+Rhyme_rs|Lemma)
             + (0+cAverage_AOA|Author_filtered) + (1|Lemma) + (1|Author_filtered),
             family=binomial, data=df, control=glmerControl(optimizer = "bobyqa"))

#R-squared and C-value
r.squaredGLMM(m5a)
Preds <- predict(m5a, type="response")
auc(df$Variant_bin, Preds) 

#CONCRETENESS
vars <- cbind(df$Variant_bin, df$clogFrequency, df$cConcreteness,
              df$cCentury, 
              df$Register, df$Lemma, df$Author_filtered, df$Rhyme, df$cLag_variation)
corvif(vars)
#ok

#full model
m_1 <- glmer(Variant ~ clogFrequency + cConcreteness + cCentury + Register + Rhyme
            + cLag_variation + (1+cCentury|Lemma) + (1+Register|Lemma) + (1+Rhyme|Lemma)
            + (1+cConcreteness|Author_filtered), data=df, family=binomial)
#no convergence
#no correlation model
m_2 <- glmer(Variant_bin ~ clogFrequency + cConcreteness + cCentury + Register + Rhyme
             + cLag_variation + (0+cCentury|Lemma) + (0+Register_rs|Lemma) + (0+Rhyme_rs|Lemma)
             + (0+cConcreteness|Author_filtered) + (1|Lemma) + (1|Author_filtered),
             data=df, family=binomial)
#no convergence
#drop register by lemma
m_3 <- glmer(Variant_bin ~ clogFrequency + cConcreteness + cCentury + Register + Rhyme
             + cLag_variation + (0+cCentury|Lemma)  + (0+Rhyme_rs|Lemma)
             + (0+cConcreteness|Author_filtered) + (1|Lemma) + (1|Author_filtered),
             data=df, family=binomial)
anova(m_2, m_3)
#m_3 better, no convergence
#drop concreteness by author
m_4 <- glmer(Variant_bin ~ clogFrequency + cConcreteness + cCentury + Register + Rhyme
             + cLag_variation + (0+cCentury|Lemma)  + (0+Rhyme_rs|Lemma)
             + (1|Lemma) + (1|Author_filtered),
             data=df, family=binomial)
anova(m_3, m_4)
#m_3 better
#drop century by lemma
m_5 <- glmer(Variant_bin ~ clogFrequency + cConcreteness + cCentury + Register + Rhyme
             + cLag_variation + (0+Rhyme_rs|Lemma)
             + (0+cConcreteness|Author_filtered) + (1|Lemma) + (1|Author_filtered),
             data=df, family=binomial)
anova(m_3, m_5)
#m_3 better
#add correlation again
#correlation century by lemma
m_6 <- glmer(Variant_bin ~ clogFrequency + cConcreteness + cCentury + Register + Rhyme
             + cLag_variation + (1+cCentury|Lemma)  + (0+Rhyme_rs|Lemma)
             + (0+cConcreteness|Author_filtered) + (1|Author_filtered),
             data=df, family=binomial)
anova(m_3, m_6)
#m_3 better
#correlation rhyme by lemma
m_7 <- glmer(Variant_bin ~ clogFrequency + cConcreteness + cCentury + Register + Rhyme
             + cLag_variation + (0+cCentury|Lemma)  + (1+Rhyme_rs|Lemma)
             + (0+cConcreteness|Author_filtered) + (1|Author_filtered),
             data=df, family=binomial)
anova(m_3, m_7)
#m_3 better
#correlation concreteness by author
m_8 <- glmer(Variant_bin ~ clogFrequency + cConcreteness + cCentury + Register + Rhyme
             + cLag_variation + (0+cCentury|Lemma)  + (0+Rhyme_rs|Lemma)
             + (1+cConcreteness|Author_filtered) + (1|Lemma),
             data=df, family=binomial)
anova(m_3, m_8)
#m_3 better
#add bobyqa
m_3a <- glmer(Variant_bin ~ clogFrequency + cConcreteness + cCentury + Register + Rhyme
              + cLag_variation + (0+cCentury|Lemma)  + (0+Rhyme_rs|Lemma)
              + (0+cConcreteness|Author_filtered) + (1|Lemma) + (1|Author_filtered),
              data=df, family=binomial, control=glmerControl(optimizer = "bobyqa"))
#still no convergence
#drop random slope with least impact, i.e. concreteness
m_3b <- glmer(Variant_bin ~ clogFrequency + cConcreteness + cCentury + Register + Rhyme
              + cLag_variation + (0+cCentury|Lemma)  + (0+Rhyme_rs|Lemma)
              + (1|Lemma) + (1|Author_filtered),
              data=df, family=binomial, control=glmerControl(optimizer = "bobyqa"))
#=> effect concreteness

#add region
m_3c <- glmer(Variant_bin ~ clogFrequency + cConcreteness + cCentury + Register + Rhyme + Region +
              + cLag_variation + (0+cCentury|Lemma)  + (0+Rhyme_rs|Lemma)
              + (1|Lemma) + (1|Author_filtered),
              data=df, family=binomial, control=glmerControl(optimizer = "bobyqa"))

#R-squared and C-value
r.squaredGLMM(m_3b)
Preds <- predict(m_3b, type="response")
auc(df$Variant_bin, Preds) 

##polysemy##
vars <- cbind(df$Variant_bin, df$clogFrequency, df$cPolysemy,
              df$cCentury, 
              df$Register, df$Lemma, df$Author_filtered, df$Rhyme, df$cLag_variation)
corvif(vars)
#ok
#full model
f0 <- glmer(Variant ~ clogFrequency + cPolysemy + cCentury + Register + Rhyme
             + cLag_variation + (1+cCentury|Lemma) + (1+Register|Lemma) + (1+Rhyme|Lemma)
             + (1+cPolysemy|Author_filtered), data=df, family=binomial)
#no convergence
#no correlation model
f1 <- glmer(Variant_bin ~ clogFrequency + cPolysemy + cCentury + Register + Rhyme
            + cLag_variation + (0+cCentury|Lemma) + (0+Register_rs|Lemma) + (0+Rhyme_rs|Lemma)
            + (0+cPolysemy|Author_filtered) + (1|Lemma) + (1|Author_filtered),
            data=df, family=binomial)
#no convergence
#drop polysemy by author
f2 <- glmer(Variant_bin ~ clogFrequency + cPolysemy + cCentury + Register + Rhyme
            + cLag_variation + (0+cCentury|Lemma) + (0+Register_rs|Lemma) + (0+Rhyme_rs|Lemma)
            + (1|Lemma) + (1|Author_filtered),
            data=df, family=binomial)
anova(f1, f2)
#f1 better
#drop register by lemma
f3 <- glmer(Variant_bin ~ clogFrequency + cPolysemy + cCentury + Register + Rhyme
            + cLag_variation + (0+cCentury|Lemma)  + (0+Rhyme_rs|Lemma)
            + (0+cPolysemy|Author_filtered) + (1|Lemma) + (1|Author_filtered),
            data=df, family=binomial)
anova(f1, f3)
#f1 better
#add correlated random slopes again
#correlation century by lemma
f4 <- glmer(Variant_bin ~ clogFrequency + cPolysemy + cCentury + Register + Rhyme
            + cLag_variation + (1+cCentury|Lemma) + (0+Register_rs|Lemma) + (0+Rhyme_rs|Lemma)
            + (0+cPolysemy|Author_filtered) + (1|Author_filtered),
            data=df, family=binomial)
anova(f1, f4)
#f1 better
#correlated random slope register by lemma
f5 <- glmer(Variant_bin ~ clogFrequency + cPolysemy + cCentury + Register + Rhyme
            + cLag_variation + (0+cCentury|Lemma) + (1+Register_rs|Lemma) + (0+Rhyme_rs|Lemma)
            + (0+cPolysemy|Author_filtered) + (1|Author_filtered),
            data=df, family=binomial)
anova(f1, f5)
#f1 better
#correlated random slope rhyme by lemma
f6 <- glmer(Variant_bin ~ clogFrequency + cPolysemy + cCentury + Register + Rhyme
            + cLag_variation + (0+cCentury|Lemma) + (0+Register_rs|Lemma) + (1+Rhyme_rs|Lemma)
            + (0+cPolysemy|Author_filtered) + (1|Author_filtered),
            data=df, family=binomial)
anova(f1, f6)
#f1 better
#correlated random slope polysemy by author
f7 <- glmer(Variant_bin ~ clogFrequency + cPolysemy + cCentury + Register + Rhyme
            + cLag_variation + (0+cCentury|Lemma) + (0+Register_rs|Lemma) + (0+Rhyme_rs|Lemma)
            + (1+cPolysemy|Author_filtered) + (1|Lemma),
            data=df, family=binomial)
#f1 better
#add bobyqa
f1a <- glmer(Variant_bin ~ clogFrequency + cPolysemy + cCentury + Register + Rhyme
             + cLag_variation + (0+cCentury|Lemma) + (0+Register_rs|Lemma) + (0+Rhyme_rs|Lemma)
             + (0+cPolysemy|Author_filtered) + (1|Lemma) + (1|Author_filtered),
             data=df, family=binomial, control=glmerControl(optimizer = "bobyqa"))
#=>no effect polysemy

#add region
f8 <- glmer(Variant_bin ~ clogFrequency + cPolysemy + cCentury + Register + Rhyme + Region +
             + cLag_variation + (0+cCentury|Lemma) + (0+Register_rs|Lemma) + (0+Rhyme_rs|Lemma)
             + (0+cPolysemy|Author_filtered) + (1|Lemma) + (1|Author_filtered),
             data=df, family=binomial, control=glmerControl(optimizer = "bobyqa"))

#R-squared and C-value
r.squaredGLMM(f1a)
Preds <- predict(f1a, type="response")
auc(df$Variant_bin, Preds) 
