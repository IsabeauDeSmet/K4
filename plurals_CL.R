#Plurals

#=============================================
#---------------------------------------------
# 1. SETTING WORKING DIRECTORY & LOADING FILES
# --------------------------------------------
#=============================================
setwd("C:/Users/u0112465/OneDrive - KU Leuven/Onderzoek/K4/meervoudsvorming")

d <- read.csv("Dataset_plurals_token_CL.csv")

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

# Fix skew in author -> too many levels with only one attestation
filter.infrequent <- function(words, threshold = 5, dummy = "OTHER") {
  # code from WBRS for recoding infrequent factor levels (default is <= 5
  # observations)
  return (relevel(
    as.factor(
      ifelse(words %in% levels(as.factor(words))[table(words) >= threshold],
             as.character(words), dummy)), dummy))
}
d$Author_filtered <- filter.infrequent(words = d$Author, threshold=5, dummy="OTHER")

d$cCentury <- scale(d$Century, center=TRUE, scale=TRUE)
d$clogFreq_meaning <- scale(d$logFreq_meaning, center=TRUE, scale=FALSE)
d$clogFreq_lemma <- scale(d$logFreq_lemma, center=TRUE, scale=FALSE)
d$clogFreq_lemma_century <- scale(d$logFreq_lemma_century, center=TRUE, scale=FALSE)
d$clag_variation <- scale(d$lag_variation, center=TRUE, scale=TRUE)

#=============================================
#---------------------------------------------
# 3. ANALYSIS TOKEN LEVEL
# --------------------------------------------
#=============================================
#VIF-scores (Zuur et al. 2009) -> oke
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


##Metaphor##
vars <- cbind(d$Variant, d$Metaphor, d$logFreq_lemma, d$cCentury,
              d$clag_variation, d$Register, d$Lemma, d$Author_filtered, d$Rhyme)

corvif(vars)

#full model
fit_a0 <- glmer(Variant ~ Metaphor +  cCentury + Register + Rhyme + clag_variation +
                 clogFreq_lemma + (1+Metaphor|Lemma) + 
                 (1+cCentury|Lemma) + (1+Register|Lemma) + (1+Rhyme|Lemma) 
               + (1+Metaphor|Author_filtered) + (1|Author_filtered) 
               + (1|Lemma), data=d, family=binomial)
#no convergence
#transform factors for no correlation model
d$Metafoor_rs <- model.matrix(fit_a0)[,2]
d$Register_rs <- model.matrix(fit_a0)[,4]
d$Rhyme_rs <- model.matrix(fit_a0)[,5]

#no correlation model
fit_a <- glmer(Variant ~ Metaphor +  cCentury + Register + Rhyme + clag_variation +
                 clogFreq_lemma + (0+Metaphor_rs|Lemma) + 
                 (0+cCentury|Lemma) + (0+Register_rs|Lemma) + (0+Rhyme_rs|Lemma) 
               + (0+Metaphor_rs|Author_filtered) + (1|Author_filtered) 
               + (1|Lemma), data=d, family=binomial)
#no convergence
#drop register by lemma
fit_b <- glmer(Variant ~ Metaphor +  cCentury + Register + Rhyme + clag_variation +
                 clogFreq_lemma + (0+Metaphor_rs|Lemma) + 
                 (0+cCentury|Lemma)  + (0+Rhyme_rs|Lemma) 
               + (0+Metaphor_rs|Author_filtered) + (1|Author_filtered) 
               + (1|Lemma), data=d, family=binomial)
anova(fit_a, fit_b)
#fit_b better, no convergence
#drop metaphor by author
fit_c <- glmer(Variant ~ Metaphor +  cCentury + Register + Rhyme + clag_variation +
                 clogFreq_lemma + (0+Metaphor_rs|Lemma) + 
                 (0+cCentury|Lemma)  + (0+Rhyme_rs|Lemma) 
               + (1|Author_filtered) 
               + (1|Lemma), data=d, family=binomial)
anova(fit_b, fit_c)
#fit_c not significantly worse, but simpler, yet no convergence
#drop metaphor by lemma
fit_d <- glmer(Variant ~ Metaphor +  cCentury + Register + Rhyme + clag_variation +
                 clogFreq_lemma  + 
                 (0+cCentury|Lemma)  + (0+Rhyme_rs|Lemma) 
               + (1|Author_filtered) 
               + (1|Lemma), data=d, family=binomial)
anova(fit_c, fit_d)
#fit_c better
#drop rhyme by lemma
fit_e <- glmer(Variant ~ Metaphor +  cCentury + Register + Rhyme + clag_variation +
                 clogFreq_lemma + (0+Metaphor_rs|Lemma) + 
                 (0+cCentury|Lemma) 
               + (1|Author_filtered) 
               + (1|Lemma), data=d, family=binomial)
anova(fit_c, fit_e)
#fit_c better
#drop century by lemma
fit_f <- glmer(Variant ~ Metaphor +  cCentury + Register + Rhyme + clag_variation +
                 clogFreq_lemma + (0+Metaphor_rs|Lemma)  + (0+Rhyme_rs|Lemma) 
               + (1|Author_filtered) 
               + (1|Lemma), data=d, family=binomial)
anova(fit_c, fit_f)
#fit_c beter
#drop random intercept author
fit_g <- glmer(Variant ~ Metaphor +  cCentury + Register + Rhyme + clag_variation +
                 clogFreq_lemma + (0+Metaphor_rs|Lemma) + 
                 (0+cCentury|Lemma)  + (0+Rhyme_rs|Lemma) 
               + (1|Lemma), data=d, family=binomial)
#fit_c better
#add correlation random slopes again
#correlation metaphor by lemma
fit_h <- glmer(Variant ~ Metaphor +  cCentury + Register + Rhyme + clag_variation +
                 clogFreq_lemma + (1+Metaphor_rs|Lemma) + 
                 (0+cCentury|Lemma)  + (0+Rhyme_rs|Lemma) 
               + (1|Author_filtered) , data=d, family=binomial)
anova(fit_c, fit_h)
#fit_c better
#correlation century by lemma
fit_i <- glmer(Variant ~ Metaphor +  cCentury + Register + Rhyme + clag_variation +
                 clogFreq_lemma + (0+Metaphor_rs|Lemma) + 
                 (1+cCentury|Lemma)  + (0+Rhyme_rs|Lemma) 
               + (1|Author_filtered) , data=d, family=binomial)
anova(fit_c, fit_i)
#fit_i better
#correlation rhyme by lemma
fit_j <- glmer(Variant ~ Metaphor +  cCentury + Register + Rhyme + clag_variation +
                 clogFreq_lemma + (0+Metaphor_rs|Lemma) + 
                 (1+cCentury|Lemma)  + (1+Rhyme_rs|Lemma) 
               + (1|Author_filtered) , data=d, family=binomial)
anova(fit_i, fit_j)
#fit_i better
#add bobyqa
fit_k <- glmer(Variant ~ Metaphor +  cCentury + Register + Rhyme + clag_variation +
                 clogFreq_lemma + (0+Metaphor_rs|Lemma) + 
                 (1+cCentury|Lemma)  + (0+Rhyme_rs|Lemma) 
               + (1|Author_filtered) , data=d, family=binomial,
               control=glmerControl(optimizer = "bobyqa"))
#=> no effect metafoor

#add region
fit_l <- glmer(Variant ~ Metaphor +  cCentury + Register + Rhyme + clag_variation + Region +
                 clogFreq_lemma + (0+Metaphor_rs|Lemma) + 
                 (1+cCentury|Lemma)  + (0+Rhyme_rs|Lemma) 
               + (1|Author_filtered) , data=d, family=binomial,
               control=glmerControl(optimizer = "bobyqa"))

#R-squared and C-value
r.squaredGLMM(fit_k)
Preds <- predict(fit_k, type="response")
auc(d$Variant, Preds) 

##Original meaning##
vars <- cbind(d$Variant, d$Original_meaning, d$logFreq_lemma, d$cCentury,
              d$clag_variation, d$Register, d$Lemma, d$Author_filtered, d$Rhyme)

corvif(vars)

#full model
fit_A0 <- glmer(Variant ~ Original_meaning +  cCentury + Register + Rhyme + clag_variation +
                  clogFreq_lemma + (1+Original_meaning|Lemma) + 
                  (1+cCentury|Lemma) + (1+Register|Lemma) + (1+Rhyme|Lemma) 
                + (1+Metaphor|Author_filtered) + (1|Author_filtered) 
                + (1|Lemma), data=d, family=binomial)
#no convergence
#transform factors for no correlation model
d$Original_meaning_rs <- model.matrix(fit_A0)[,2]

fit_A <- glmer(Variant ~ Original_meaning +  cCentury + Register + Rhyme + clag_variation +
                 clogFreq_lemma + (0+Original_meaning_rs|Lemma) + 
                 (0+cCentury|Lemma) + (0+Register_rs|Lemma) + (0+Rhyme_rs|Lemma) 
               + (0+Original_meaning_rs|Author_filtered) + (1|Author_filtered) 
               + (1|Lemma), data=d, family=binomial)
#no convergence
#drop register by lemma
fit_B <- glmer(Variant ~ Original_meaning_rs +  cCentury + Register + Rhyme + clag_variation +
                 clogFreq_lemma + (0+Original_meaning_rs|Lemma) + 
                 (0+cCentury|Lemma)  + (0+Rhyme_rs|Lemma) 
               + (0+Original_meaning_rs|Author_filtered) + (1|Author_filtered) 
               + (1|Lemma), data=d, family=binomial)
anova(fit_A, fit_B)
#fit_B better
#drop original meaning by author
fit_C <- glmer(Variant ~ Original_meaning_rs +  cCentury + Register + Rhyme + clag_variation +
                 clogFreq_lemma + (0+Original_meaning_rs|Lemma) + 
                 (0+cCentury|Lemma)  + (0+Rhyme_rs|Lemma) + (1|Author_filtered) 
               + (1|Lemma), data=d, family=binomial)
anova(fit_B, fit_C)
#fit_C better
#drop original meaning by lemma
fit_D <- glmer(Variant ~ Original_meaning_rs +  cCentury + Register + Rhyme + clag_variation +
                 clogFreq_lemma + 
                 (0+cCentury|Lemma)  + (0+Rhyme_rs|Lemma) + (1|Author_filtered) 
               + (1|Lemma), data=d, family=binomial)
anova(fit_C, fit_D)
#fit_C better
# drop rhyme by lemma
fit_E <- glmer(Variant ~ Original_meaning_rs +  cCentury + Register + Rhyme + clag_variation +
                 clogFreq_lemma + (0+Original_meaning_rs|Lemma) + 
                 (0+cCentury|Lemma) + (1|Author_filtered) 
               + (1|Lemma), data=d, family=binomial)
anova(fit_C, fit_E)
#fit_E better
#drop century by lemma
fit_F <- glmer(Variant ~ Original_meaning_rs +  cCentury + Register + Rhyme + clag_variation +
                 clogFreq_lemma + (0+Original_meaning_rs|Lemma) + (1|Author_filtered) 
               + (1|Lemma), data=d, family=binomial)
anova(fit_E, fit_F)
#fit_E better
#drop intercept author
fit_G <- glmer(Variant ~ Original_meaning_rs +  cCentury + Register + Rhyme + clag_variation +
                 clogFreq_lemma + 
                 (0+cCentury|Lemma) + (0+Original_meaning_rs|Lemma) 
               + (1|Lemma), data=d, family=binomial)
anova(fit_F, fit_G)
#fit_E better
#add correlation ran slopes again
#correlation meaning by lemma
fit_H <- glmer(Variant ~ Original_meaning_rs +  cCentury + Register + Rhyme + clag_variation +
                 clogFreq_lemma + (1+Original_meaning_rs|Lemma) + 
                 (0+cCentury|Lemma) + (1|Author_filtered) , data=d, family=binomial)
anova(fit_H, fit_E)
#fit_E beter
fit_I <- glmer(Variant ~ Original_meaning_rs +  cCentury + Register + Rhyme + clag_variation +
                 clogFreq_lemma + (0+Original_meaning_rs|Lemma) + 
                 (1+cCentury|Lemma) + (1|Author_filtered) , data=d, family=binomial)
anova(fit_I, fit_E)
#fit_I better
#add bobyqa
fit_J <- glmer(Variant ~ Original_meaning_rs +  cCentury + Register + Rhyme + clag_variation +
                 clogFreq_lemma + (0+Original_meaning_rs|Lemma) + 
                 (1+cCentury|Lemma) + (1|Author_filtered) , data=d, family=binomial, control=glmerControl(optimizer = "bobyqa"))
#=> no effect original meaning

#add region
fit_K <- glmer(Variant ~ Original_meaning_rs +  cCentury + Register + Rhyme + clag_variation + Region +
                 clogFreq_lemma + (0+Original_meaning_rs|Lemma) + 
                 (1+cCentury|Lemma) + (1|Author_filtered) , data=d, family=binomial, control=glmerControl(optimizer = "bobyqa"))

#R-squared and C-value
r.squaredGLMM(fit_J)
Preds <- predict(fit_J, type="response")
auc(d$Variant, Preds) 

##meaning freq##
#VIF's
vars <- cbind(d$Variant, d$Freq_meaning_cat, d$logFreq_lemma, d$cCentury,
              d$clag_variation, d$Register, d$Lemma, d$Author_filtered, d$Rhyme)
corvif(vars)

#full model
model0 <- glmer(Variant ~ Freq_meaning_cat +  cCentury + Register + Rhyme + clag_variation +
                  clogFreq_lemma + (1+Freq_meaning_cat|Lemma) + 
                  (1+cCentury|Lemma) + (1+Register|Lemma) + (1+Rhyme|Lemma) 
                + (1+Freq_meaning_cat|Author_filtered) + (1|Author_filtered) 
                + (1|Lemma), data=d, family=binomial)
#no convergence
#transform factors for no correlation model
d$Freq_meaning_cat_rs <- model.matrix(model0)[,2]
#no correlation model
model1 <- glmer(Variant ~ Freq_meaning_cat +  cCentury + Register + Rhyme + clag_variation +
                  clogFreq_lemma + (0+Freq_meaning_cat_rs|Lemma) + 
                  (0+cCentury|Lemma) + (0+Register_rs|Lemma) + (0+Rhyme_rs|Lemma) 
                + (0+Freq_meaning_cat_rs|Author_filtered) + (1|Author_filtered) 
                + (1|Lemma), data=d, family=binomial)
#no convergence
#drop register by lemma
model2 <- glmer(Variant ~ Freq_meaning_cat +  cCentury + Register + Rhyme + clag_variation +
                  clogFreq_lemma + (0+Freq_meaning_cat_rs|Lemma) + 
                  (0+cCentury|Lemma)  + (0+Rhyme_rs|Lemma) 
                + (0+Freq_meaning_cat_rs|Author_filtered) + (1|Author_filtered) 
                + (1|Lemma), data=d, family=binomial)
anova(model1, model2)
#model2 better, no convergence
#drop meaning freq by author
model3 <- glmer(Variant ~ Freq_meaning_cat +  cCentury + Register + Rhyme + clag_variation +
                  clogFreq_lemma + (0+Freq_meaning_cat_rs|Lemma) + 
                  (0+cCentury|Lemma)  + (0+Rhyme_rs|Lemma)  + (1|Author_filtered) 
                + (1|Lemma), data=d, family=binomial)
anova(model2, model3)
#model3 better, no convergence
#drop meaning freq by lemma
model4 <- glmer(Variant ~ Freq_meaning_cat +  cCentury + Register + Rhyme + clag_variation +
                  clogFreq_lemma + 
                  (0+cCentury|Lemma)  + (0+Rhyme_rs|Lemma)  + (1|Author_filtered) 
                + (1|Lemma), data=d, family=binomial)
anova(model3, model4)
#model3 better,
#drop rhyme by lemma
model5 <- glmer(Variant ~ Freq_meaning_cat +  cCentury + Register + Rhyme + clag_variation +
                  clogFreq_lemma + (0+Freq_meaning_cat_rs|Lemma) + 
                  (0+cCentury|Lemma) + (1|Author_filtered) 
                + (1|Lemma), data=d, family=binomial)
anova(model3, model5)
#model3 better
#drop century by lemma
model6 <- glmer(Variant ~ Freq_meaning_cat +  cCentury + Register + Rhyme + clag_variation +
                  clogFreq_lemma + (0+Freq_meaning_cat_rs|Lemma) + (0+Rhyme_rs|Lemma)  + (1|Author_filtered) 
                + (1|Lemma), data=d, family=binomial)
anova(model3, model6)
#model3 better
#add correlated random slopes again
#correlation freq meaning by lemma
model7 <- glmer(Variant ~ Freq_meaning_cat +  cCentury + Register + Rhyme + clag_variation +
                  clogFreq_lemma + (1+Freq_meaning_cat_rs|Lemma) + 
                  (0+cCentury|Lemma)  + (0+Rhyme_rs|Lemma)  + (1|Author_filtered), data=d, family=binomial)
anova(model3, model7)
#model3 better
#correlation century by lemma
model8 <- glmer(Variant ~ Freq_meaning_cat +  cCentury + Register + Rhyme + clag_variation +
                  clogFrequentie_lemma + (0+Freq_meaning_cat_rs|Lemma) + 
                  (1+cCentury|Lemma)  + (0+Rhyme_rs|Lemma)  + (1|Author_filtered), data=d, family=binomial)
anova(model3, model8)
#model8 better
#correlation rhyme by lemma
model9 <- glmer(Variant ~ Freq_meaning_cat +  cCentury + Register + Rhyme + clag_variation +
                  clogFreq_lemma + (0+Freq_meaning_cat_rs|Lemma) + 
                  (1+cCentury|Lemma)  + (1+Rhyme_rs|Lemma)  + (1|Author_filtered), data=d, family=binomial)
#add bobyqa
model10 <- glmer(Variant ~ Freq_meaning_cat +  cCentury + Register + Rhyme + clag_variation +
                   clogFreq_lemma + (0+Freq_meaning_cat_rs|Lemma) + 
                   (1+cCentury|Lemma)  + (0+Rhyme_rs|Lemma)  + (1|Author_filtered), data=d, family=binomial, control=glmerControl(optimizer = "bobyqa"))
#=> no effect freq meaning

#add region
model11 <- glmer(Variant ~ Freq_meaning_cat +  cCentury + Register + Rhyme + clag_variation + Region +
                   clogFreq_lemma + (0+Freq_meaning_cat_rs|Lemma) + 
                   (1+cCentury|Lemma)  + (0+Rhyme_rs|Lemma)  + (1|Author_filtered), data=d, family=binomial, control=glmerControl(optimizer = "bobyqa"))

#R-squared and C-value
r.squaredGLMM(model10)
Preds <- predict(model10, type="response")
auc(d$Variant, Preds) 

##Concreteness##
#VIF's
vars <- cbind(d$Variant, d$Concreteness, d$logFreq_lemma, d$cCentury,
              d$clag_variation, d$Register, d$Lemma, d$Author_filtered, d$Rhyme)
corvif(vars)

#full model
modela0 <- glmer(Variant ~ Concreteness +  cCentury + Register + Rhyme + clag_variation +
                  clogFreq_lemma + (1+Concreteness|Lemma) + 
                  (1+cCentury|Lemma) + (1+Register|Lemma) + (1+Rhyme|Lemma) 
                + (1+Concreteness|Author_filtered) + (1|Author_filtered) 
                + (1|Lemma), data=d, family=binomial)
#no convergence
#transform factors for no correlation model
d$Concreteness_rs <- model.matrix(modela0)[,2]
#no correlation model
modela <- glmer(Variant ~ Concreteness +  cCentury + Register + Rhyme + clag_variation +
                  clogFreq_lemma + (0+Concreteness_rs|Lemma) + 
                  (0+cCentury|Lemma) + (0+Register_rs|Lemma) + (0+Rhyme_rs|Lemma) 
                + (0+Concreteness_rs|Author_filtered) + (1|Author_filtered) 
                + (1|Lemma), data=d, family=binomial)
#no convergence
#drop register by lemma
modelb <- glmer(Variant ~ Fysiek +  cCentury + Register + Rhyme + clag_variation +
                  clogFreq_lemma + (0+Concreteness_rs|Lemma) + 
                  (0+cCentury|Lemma)  + (0+Rhyme_rs|Lemma) 
                + (0+Concreteness_rs|Author_filtered) + (1|Author_filtered) 
                + (1|Lemma), data=d, family=binomial)
anova(modela, modelb)
#modelb better
#drop concreteness by lemma
modelc <- glmer(Variant ~ Concreteness +  cCentury + Register + Rhyme + clag_variation +
                  clogFreq_lemma + 
                  (0+cCentury|Lemma)  + (0+Rhyme_rs|Lemma) 
                + (0+Concreteness_rs|Author_filtered) + (1|Author_filtered) 
                + (1|Lemma), data=d, family=binomial)
anova(modelb, modelc)
#modelc better
#drop concreteness by author
modeld <- glmer(Variant ~ Fysiek +  cCentury + Register + Rhyme + clag_variation +
                  clogFreq_lemma + 
                  (0+cCentury|Lemma)  + (0+Rhyme_rs|Lemma) 
                + (1|Author_filtered) 
                + (1|Lemma), data=d, family=binomial)
anova(modelc, modeld)
#modelc better
#drop rhyme by lemma
modele <- glmer(Variant ~ Concreteness +  cCentury + Register + Rhyme + clag_variation +
                  clogFreq_lemma + 
                  (0+cCentury|Lemma)  
                + (0+Concreteness_rs|Author_filtered) + (1|Author_filtered) 
                + (1|Lemma), data=d, family=binomial)
anova(modelc, modele)
#modelc better
#drop century by lemma
modelf <- glmer(Variant ~ Concreteness +  cCentury + Register + Rhyme + clag_variation +
                  clogFreq_lemma  + (0+Rhyme_rs|Lemma) 
                + (0+Concreteness_rs|Author_filtered) + (1|Author_filtered) 
                + (1|Lemma), data=d, family=binomial)
anova(modelc, modelf)
#modelc better
#add correlation ran slopes again
#correlation century by lemma
modelg <- glmer(Variant ~ Concreteness +  cCentury + Register + Rhyme + clag_variation +
                  clogFreq_lemma + 
                  (1+cCentury|Lemma)  + (0+Rhyme_rs|Lemma) 
                + (0+Concreteness_rs|Author_filtered) + (1|Author_filtered), data=d, family=binomial)
anova(modelc, modelg)
#modelg better
#correlation rhyme by lemma
modelh <- glmer(Variant ~ Concreteness +  cCentury + Register + Rhyme + clag_variation +
                  clogFreq_lemma + 
                  (1+cCentury|Lemma)  + (1+Rhyme_rs|Lemma) 
                + (0+Concreteness_rs|Author_filtered) + (1|Author_filtered), data=d, family=binomial)µ
anova(modelg, modelh)
#modelg better
#correlation concreteness by author
modeli <- glmer(Variant ~ Concreteness +  cCentury + Register + Rhyme + clag_variation +
                  clogFreq_lemma + 
                  (1+cCentury|Lemma)  + (0+Rhyme_rs|Lemma) 
                + (1+Concreteness_rs|Author_filtered), data=d, family=binomial)
anova(modelg, modeli)
#modeli better
#add bobyqa
modelj <- glmer(Variant ~ Concreteness +  cCentury + Register + Rhyme + clag_variation +
                  clogFreq_lemma + 
                  (1+cCentury|Lemma)  + (0+Rhyme_rs|Lemma) 
                + (1+Concreteness_rs|Author_filtered), data=d, family=binomial, control=glmerControl(optimizer = "bobyqa"))
#=> effect concreteness

#add region
modelk <- glmer(Variant ~ Concreteness +  cCentury + Register + Rhyme + clag_variation +
                  clogFreq_lemma + Region +
                  (1+cCentury|Lemma)  + (0+Rhyme_rs|Lemma) 
                + (1+Concreteness_rs|Author_filtered), data=d, family=binomial, control=glmerControl(optimizer = "bobyqa"))

#R-squared and C-value
r.squaredGLMM(modelj)
Preds <- predict(modelj, type="response")
auc(d$Variant, Preds) 

#=============================================
#---------------------------------------------
# 4. PREP ANALYSIS TYPE LEVEL
# --------------------------------------------
#=============================================

df$clag_variation <- scale(df$lag_variation, center=TRUE, scale=TRUE)
df$cCentury <- scale(df$Century, center=TRUE, scale=TRUE)
df$logFrequency <- log10(df$Frequency)
df$clogFrequency <- scale(df$logFrequency, center=TRUE, scale=FALSE)
df$cAverage_AOA <- scale(df$Average_AOA, center=TRUE, scale=TRUE)
df$cConcreteness <- scale(df$Concreteness, center=TRUE, scale=TRUE)
df$cLength <- scale(df$Length, center=TRUE, scale=TRUE)
df$cPolysemy <- scale(df$Polysemy, center=TRUE, scale=TRUE)

#Decrease number of levels for author
filter.infrequent <- function(words, threshold = 5, dummy = "OTHER") {
  # code from WBRS for recoding infrequent factor levels (default is <= 5
  # observations)
  return (relevel(
    as.factor(
      ifelse(words %in% levels(as.factor(words))[table(words) >= threshold],
             as.character(words), dummy)), dummy))
}
df$Author_filtered <- filter.infrequent(words = df$Auteur, threshold=5, dummy="OTHER")
#=============================================
#---------------------------------------------
# 5. ANALYSIS TYPE LEVEL
# --------------------------------------------
#=============================================

##Age of acquisition##
#VIF's
vars <- cbind(d$Variant, d$cAverage_AOA, d$clogFrequency, d$cCentury,
              d$clag_variation, d$Register, d$Lemma, d$Author_filtered, d$Rhyme)
corvif(vars)

#full model
m0 <- glmer(Variant ~ cCentury + clogFrequency + cAverage_AOA + clag_variation + Register + Rhyme + (1+cCentury|Lemma) 
            + (1+Register_rs|Lemma) + (1+Rhyme_rs|Lemma) + (1+cAverage_AOA|Author_filtered)
            ,data=df, family=binomial)
#no convergence
#drop correlation parameters
#transform factor
df$Rhyme <- model.matrix(m0)[,7]
df$Register <- model.matrix(m0)[,6]

#no correlation model
m1 <- glmer(Variant ~ cCentury + clogFrequency + cAverage_AOA + clag_variation + Register + Rhyme + (0+cCentury|Lemma) 
            + (0+Register_rs|Lemma) + (0+Rhyme_rs|Lemma) + (0+cAverage_AOA|Author_filtered) + (1|Author_filtered)
            + (1|Lemma),
            data=df, family=binomial)
#no convergence
#drop AoA by author
m2 <- glmer(Variant ~ cCentury + clogFrequency + cAverage_AOA + clag_variation + Register + Rhyme + (0+cCentury|Lemma) 
            + (0+Register_rs|Lemma) + (0+Rhyme_rs|Lemma)  + (1|Author_filtered)
            + (1|Lemma),
            data=df, family=binomial)
anova(m1, m2)
#m2 better
#drop rhyme by lemma
m3 <- glmer(Variant ~ cCentury + clogFrequency + cAverage_AOA + clag_variation + Register + Rhyme + (0+cCentury|Lemma) 
            + (0+Register_rs|Lemma)  + (1|Author_filtered)
            + (1|Lemma),
            data=df, family=binomial)
anova(m2, m3)
#m3 not worse, but simpler, yet still no convergence
#drop register by lemma
m4 <- glmer(Variant ~ cCentury + clogFrequency + cAverage_AOA + clag_variation + Register + Rhyme + (0+cCentury|Lemma) 
            + (1|Author_filtered)
            + (1|Lemma),
            data=df, family=binomial)
anova(m3, m4)
#m3 better
#drop century by lemma
m5 <- glmer(Variant ~ cCentury + clogFrequency + cAverage_AOA + clag_variation + Register + Rhyme 
            + (0+Register_rs|Lemma)  + (1|Author_filtered)
            + (1|Lemma),
            data=df, family=binomial)
anova(m3, m5)
#m3 better
#add correlation random slopes again
#correlation century by lemma
m6 <- glmer(Variant ~ cCentury + clogFrequency + cAverage_AOA + clag_variation + Register + Rhyme + (1+cCentury|Lemma) 
            + (0+Register_rs|Lemma)  + (1|Author_filtered),
            data=df, family=binomial)
anova(m3, m6)
#m6 better
m7 <- glmer(Variant ~ cCentury + clogFrequency + cAverage_AOA + clag_variation + Register + Rhyme + (1+cCentury|Lemma) 
            + (1+Register_rs|Lemma)  + (1|Author_filtered),
            data=df, family=binomial)
anova(m6, m7)
#m6 better
#add bobyqa
m6a <- glmer(Variant ~ cCentury + clogFrequency + cAverage_AOA + clag_variation + Register + Rhyme + (1+cCentury|Lemma) 
             + (0+Register_rs|Lemma)  + (1|Author_filtered),
             data=df, family=binomial,
             control=glmerControl(optimizer = "bobyqa"))
#=> no effect of Age of acquisition

#add region
m8 <- glmer(Variant ~ cCentury +  Region + clogFrequency + cAverage_AOA + clag_variation + Register + Rhyme + (1+cCentury|Lemma) 
            + (0+Register_rs|Lemma)  + (1|Author_filtered),
            data=df, family=binomial)

#R-squared and C-value
r.squaredGLMM(m6a)
Preds <- predict(m6a, type="response")
auc(d$Variant, Preds) 


##Concreteness##
#VIF's
vars <- cbind(d$Variant, d$cConcreteness, d$clogFrequency, d$cCentury,
              d$clag_variation, d$Register, d$Lemma, d$Author_filtered, d$Rhyme)
corvif(vars)

#full model
m0 <- glmer(Variant ~ cCentury + clogFrequency + cConcreteness + clag_variation + Register + Rhyme + (1+cCentury|Lemma) 
            + (1+Register_rs|Lemma) + (1+Rhyme_rs|Lemma) + (1+cConcreteness|Author_filtered)
            ,data=df, family=binomial)
#no convergence
#drop correlation parameters
m_1 <- glmer(Variant ~ cCentury + clogFrequency + cConcreteness + clag_variation + Register + Rhyme + (0+cCentury|Lemma) 
             + (0+Register_rs|Lemma) + (0+Rhyme_rs|Lemma) + (0+cConcreteness|Author_filtered) + (1|Author_filtered)
             + (1|Lemma),
             data=df, family=binomial)
#no convergence
#drop concreteness by author
m_2 <- glmer(Variant ~ cCentury + clogFrequency + cConcreteness + clag_variation + Register + Rhyme + (0+cCentury|Lemma) 
             + (0+Register_rs|Lemma) + (0+Rhyme_rs|Lemma) + (1|Author_filtered)
             + (1|Lemma),
             data=df, family=binomial)
anova(m_1, m_2)
#m_1 better
#drop rhyme by lemma
m_3 <- glmer(Variant ~ cCentury + clogFrequency + cConcreteness + clag_variation + Register + Rhyme + (0+cCentury|Lemma) 
             + (0+Register_rs|Lemma) + (0+cConcreteness|Author_filtered) + (1|Author_filtered)
             + (1|Lemma),
             data=df, family=binomial)
anova(m_2, m_3)
#m_3 better
#drop register by lemma
m_4 <- glmer(Variant ~ cCentury + clogFrequency + cConcreteness + clag_variation + Register + Rhyme + (0+cCentury|Lemma) 
             + (0+cConcreteness|Author_filtered) + (1|Author_filtered)
             + (1|Lemma),
             data=df, family=binomial)
anova(m_3, m_4)
#m_3 better
#add correlation random slopes again
#add correlation century by lemma
m_5<- glmer(Variant ~ cCentury + clogFrequency + cConcreteness + clag_variation + Register + Rhyme + (1+cCentury|Lemma) 
            + (0+Register_rs|Lemma) + (0+cConcreteness|Author_filtered) + (1|Author_filtered)
            , data=df, family=binomial)
anova(m_3, m_5)
#m_5 better
m_6 <- glmer(Variant ~ cCentury + clogFrequency + cConcreteness + clag_variation + Register + Rhyme + (1+cCentury|Lemma) 
             + (1+Register_rs|Lemma) + (0+cConcreteness|Author_filtered) + (1|Author_filtered)
             ,data=df, family=binomial)
anova(m_5, m_6)
#m_5 better
m_7<- glmer(Variant ~ cCentury + clogFrequency + cConcreteness + clag_variation + Register + Rhyme + (1+cCentury|Lemma) 
            + (0+Register_rs|Lemma) + (1+cConcreteness|Author_filtered) 
            ,data=df, family=binomial)
anova(m_5, m_7)
#m_5 better
#add bobyqa
m_5a <- glmer(Variant ~ cCentury + clogFrequency + cConcreteness + clag_variation + Register + Rhyme + (1+cCentury|Lemma) 
              + (0+Register_rs|Lemma) + (0+cConcreteness|Author_filtered) + (1|Author_filtered)
              ,data=df, family=binomial,
              control=glmerControl(optimizer = "bobyqa"))
#=> effect concreteness


#add region
m_8<- glmer(Variant ~ cCentury + Region + clogFrequency + cConcreteness + clag_variation + Register + Rhyme + (1+cCentury|Lemma) 
            + (0+Register_rs|Lemma) + (0+cConcreteness|Author_filtered) + (1|Author_filtered)
            ,data=df, family=binomial)

#R-squared and C-value
r.squaredGLMM(m_8)
Preds <- predict(m_8, type="response")
auc(d$Variant, Preds) 


##POLYSEMY##
#VIF's
vars <- cbind(d$Variant, d$cPolysemy, d$clogFrequency, d$cCentury,
              d$clag_variation, d$Register, d$Lemma, d$Author_filtered, d$Rhyme)
corvif(vars)

#full model
a0 <- glmer(Variant ~ cCentury + clogFrequency + cPolysemy + clag_variation + Register + Rhyme + (1+cCentury|Lemma) 
            + (1+Register_rs|Lemma) + (1+Rhyme_rs|Lemma) + (1+cPolysemy|Author_filtered)
            ,data=df, family=binomial)
#no convergence
#no correlation model
a1 <- glmer(Variant ~ cCentury + clogFrequency + cPolysemy + clag_variation + Register + Rhyme + (0+cCentury|Lemma) 
            + (0+Register_rs|Lemma) + (0+Rhyme_rs|Lemma) + (0+cPolysemy|Author_filtered) + (1|Author_filtered)
            + (1|Lemma),
            data=df, family=binomial)
#no convergence
#drop polysemy by author
a2 <- glmer(Variant ~ cCentury + clogFrequency + cPolysemy + clag_variation + Register + Rhyme + (0+cCentury|Lemma) 
            + (0+Register_rs|Lemma) + (0+Rhyme_rs|Lemma)  + (1|Author_filtered)
            + (1|Lemma),
            data=df, family=binomial)
anova(a1, a2)
#a1 better
#drop rhyme by lemma
a3 <- glmer(Variant ~ cCentury + clogFrequency + cPolysemy + clag_variation + Register + Rhyme + (0+cCentury|Lemma) 
            + (0+Register_rs|Lemma) + (0+cPolysemy|Author_filtered) + (1|Author_filtered)
            + (1|Lemma),
            data=df, family=binomial)
anova(a1, a3)
#a3 better
#drop regsiter by lemma
a4 <- glmer(Variant ~ cCentury + clogFrequency + cPolysemy + clag_variation + Register + Rhyme + (0+cCentury|Lemma) 
            + (0+cPolysemy|Author_filtered) + (1|Author_filtered)
            + (1|Lemma),
            data=df, family=binomial)
#a3 better
#add correlation random slopes again
#correlation century by lemma
a5 <- glmer(Variant ~ cCentury + clogFrequency + cPolysemy + clag_variation + Register + Rhyme + (1+cCentury|Lemma) 
            + (0+Register_rs|Lemma) + (0+cPolysemy|Author_filtered) + (1|Author_filtered)
            ,data=df, family=binomial)
anova(a3, a5)
#a5 better
a6 <- glmer(Variant ~ cCentury + clogFrequency + cPolysemy + clag_variation + Register + Rhyme + (1+cCentury|Lemma) 
            + (1+Register_rs|Lemma) + (0+cPolysemy|Author_filtered) + (1|Author_filtered)
            , data=df, family=binomial)
anova(a5, a6)
#a5 beter
a7 <- glmer(Variant ~ cCentury + clogFrequency + cPolysemy + clag_variation + Register + Rhyme + (1+cCentury|Lemma) 
            + (0+Register_rs|Lemma) + (1+cPolysemy|Author_filtered), data=df, family=binomial)
#a5 better
#add bobyqa
a5a <- glmer(Variant ~ cCentury + clogFrequency + cPolysemy + clag_variation + Register + Rhyme + (1+cCentury|Lemma) 
             + (0+Register_rs|Lemma) + (0+cPolysemy|Author_filtered) + (1|Author_filtered)
             ,data=df, family=binomial,
             control=glmerControl(optimizer = "bobyqa"))
#=> no effect polysemy

#add region
a8 <- glmer(Variant ~ cCentury + Regio + clogFrequency + cPolysemy + Region + clag_variation + Register + Rhyme + (1+cCentury|Lemma) 
            + (0+Register_rs|Lemma) + (0+cPolysemy|Author_filtered) + (1|Author_filtered)
            ,data=df, family=binomial)

#R-squared and C-value
r.squaredGLMM(a5a)
Preds <- predict(a5a, type="response")
auc(d$Variant, Preds) 