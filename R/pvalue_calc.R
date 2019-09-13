#' Calculate the p values from longitudinal data
#' @param file This is the data file that has been melted to so that predictors
#'             are created in advance
#' @param p_table_name The name you'd like your p_value table to have
#' @return None
#' @examples pvalue_calc(file, p_table_name)
#' @export
#'
#' @import lme4
#' @import orddom
#' @import reshape
#' @import lmtest

pvalue_calc <- function(file, p_table_name){
  variables <- unique (melt_data$variable)
  N <- length (variables)
  Ps <- data.frame (variable = variables, cP12 = rep (NA, N), cwP12 = rep (NA, N),
                    cd12 = rep (NA, N), cP13 = rep (NA, N), cwP13 = rep (NA, N),
                    cd13 = rep (NA, N),  cP23 = rep (NA, N), cwP23 = rep (NA, N),
                    cd23 = rep (NA, N), RcP12 = rep (NA, N), RcwP12 = rep (NA, N),
                    Rcd12 = rep (NA, N), RcP13 = rep (NA, N), RcwP13 = rep (NA, N),
                    Rcd13 = rep (NA, N),  RcP23 = rep (NA, N), RcwP23 = rep (NA, N),
                    Rcd23 = rep (NA, N), NRcP12 = rep (NA, N), NRcwP12 = rep (NA, N),
                    NRcd12 = rep (NA, N), NRcP13 = rep (NA, N), NRcwP13 = rep (NA, N),
                    NRcd13 = rep (NA, N),  NRcP23 = rep (NA, N), NRcwP23 = rep (NA, N),
                    NRcd23 = rep (NA, N), cPRNR1=rep(NA,N), cwPRNR1=rep(NA,N), cdRNR1=rep(NA,N),
                    pModel = rep(NA,N), pPost12  = rep(NA,N), pPost13  = rep(NA,N),
                    pPost23  = rep(NA,N), RpModel = rep(NA,N), RpPost12  = rep(NA,N),
                    RpPost13  = rep(NA,N), RpPost23  = rep(NA,N), cPRNR2=rep(NA,N),
                    cwPRNR2=rep(NA,N), cdRNR2=rep(NA,N), cPRNR3=rep(NA,N), cwPRNR3=rep(NA,N),
                    cdRNR3=rep(NA,N))
  for (i in 1:N) {
    # loop through all variables
    print(i)
    aVariable = variables [i]
    print(aVariable)
    # corona block
    subdatac <- subset (melt_data, variable == aVariable & corona.dash == "c")
    sub1c <- subset (subdatac, Case == 1)
    sub2c <- subset (subdatac, Case == 2)
    sub3c <- subset (subdatac, Case == 3)
    #Select the data present in every timepoint
    subdatac <- subset (subdatac, Patient %in% intersect (unique (sub1c$Patient), intersect (unique (sub2c$Patient), unique (sub3c$Patient))))
    sub1c <- subset (subdatac, Case == 1)
    sub2c <- subset (subdatac, Case == 2)
    sub3c <- subset (subdatac, Case == 3)

    subdatac$value <- as.numeric (subdatac$value)
    sub1c$value <- as.numeric (sub1c$value)
    sub2c$value <- as.numeric (sub2c$value)
    sub3c$value <- as.numeric (sub3c$value)

    if (! (aVariable %in% varsNotInCase2)) { # rule out the variables not defined in case 2
      # cliffs delta orddom paired [11,1]^=cliffs delta (orddom output much more measures)
      d12 <- as.numeric (orddom::orddom (x = (sub1c [! is.na (sub1c$value) & ! is.na (sub2c$value), ])$value, y = (sub2c [! is.na (sub1c$value) & ! is.na (sub2c$value), ])$value, paired = T)[11,1])
      # wilcoxon
      pw12 <- wilcox.test (sub1c$value, sub2c$value, paired = T)$p.value
      m1 <- lm (data = subset (subdatac, Case == 1 | Case == 2), value ~ Patient + age + sex + ACE_INHIBITOR + ANTIDEPRESSANT + ANTIDIABETIC + ANTI_PLATELET + AT2_BLOCKER + BETA_BLOCKER + CALCIUMANTAGONIST + CORTICOSTEROID + DIURETIC + IMIDAZOLE_AGONIST + METFORMIN + PPI + STATIN + THYROID + VITAMIN_D + XANTHAN_OXIDASE_INHIBITOR)
      # exclude case
      m2 <- lm (data = subset (subdatac, Case == 1 | Case == 2), value ~ Patient + Case + age + sex + ACE_INHIBITOR + ANTIDEPRESSANT + ANTIDIABETIC + ANTI_PLATELET + AT2_BLOCKER + BETA_BLOCKER + CALCIUMANTAGONIST + CORTICOSTEROID + DIURETIC + IMIDAZOLE_AGONIST + METFORMIN + PPI + STATIN + THYROID + VITAMIN_D + XANTHAN_OXIDASE_INHIBITOR)
      # likelyhoodratio
      p12 <- lmtest::lrtest (m1, m2)$"Pr(>Chisq)" [2]

    }

    # cliff's delta
    d13 <- as.numeric (orddom::orddom (x = (sub1c [! is.na (sub1c$value) & ! is.na (sub3c$value), ])$value, y = (sub3c [! is.na (sub1c$value) & ! is.na (sub3c$value), ])$value, paired = T)[11,1])
    # wilcoxon
    pw13 <- wilcox.test (sub1c$value, sub3c$value, paired = T)$p.value

    m1 <- lm (data = subset (subdatac, Case == 1 | Case == 3), value ~ Patient + age + sex + ACE_INHIBITOR + ANTIDEPRESSANT + ANTIDIABETIC + ANTI_PLATELET + AT2_BLOCKER + BETA_BLOCKER + CALCIUMANTAGONIST + CORTICOSTEROID + DIURETIC + IMIDAZOLE_AGONIST + METFORMIN + PPI + STATIN + THYROID + VITAMIN_D + XANTHAN_OXIDASE_INHIBITOR)
    m2 <- lm (data = subset (subdatac, Case == 1 | Case == 3), value ~ Patient + Case + age + sex + ACE_INHIBITOR + ANTIDEPRESSANT + ANTIDIABETIC + ANTI_PLATELET + AT2_BLOCKER + BETA_BLOCKER + CALCIUMANTAGONIST + CORTICOSTEROID + DIURETIC + IMIDAZOLE_AGONIST + METFORMIN + PPI + STATIN + THYROID + VITAMIN_D + XANTHAN_OXIDASE_INHIBITOR)
    p13 <- lmtest::lrtest (m1, m2)$"Pr(>Chisq)" [2]


    if (! (aVariable %in% varsNotInCase2)) { # rule out the variables not defined in case 2
      d23 <- as.numeric (orddom::orddom (x = (sub2c [! is.na (sub2c$value) & ! is.na (sub3c$value), ])$value, y = (sub3c [! is.na (sub2c$value) & ! is.na (sub3c$value), ])$value, paired = T)[11,1])
      pw23 <- wilcox.test (sub2c$value, sub3c$value, paired = T)$p.value

      m1 <- lm (data = subset (subdatac, Case == 2 | Case == 3), value ~ Patient + age + sex + ACE_INHIBITOR + ANTIDEPRESSANT + ANTIDIABETIC + ANTI_PLATELET + AT2_BLOCKER + BETA_BLOCKER + CALCIUMANTAGONIST + CORTICOSTEROID + DIURETIC + IMIDAZOLE_AGONIST + METFORMIN + PPI + STATIN + THYROID + VITAMIN_D + XANTHAN_OXIDASE_INHIBITOR)
      m2 <- lm (data = subset (subdatac, Case == 2 | Case == 3), value ~ Patient + Case + age + sex + ACE_INHIBITOR + ANTIDEPRESSANT + ANTIDIABETIC + ANTI_PLATELET + AT2_BLOCKER + BETA_BLOCKER + CALCIUMANTAGONIST + CORTICOSTEROID + DIURETIC + IMIDAZOLE_AGONIST + METFORMIN + PPI + STATIN + THYROID + VITAMIN_D + XANTHAN_OXIDASE_INHIBITOR)
      p23 <- lmtest::lrtest (m1, m2)$"Pr(>Chisq)" [2]

    }

    Ps [i, 2] <- p12
    Ps [i, 3] <- pw12
    Ps [i, 4] <- d12

    Ps [i, 5] <- p13
    Ps [i, 6] <- pw13
    Ps [i, 7] <- d13

    Ps [i, 8] <- p23
    Ps [i, 9] <- pw23
    Ps [i, 10] <- d23

    print("passed block 1")

    #### Then repeat all for responders only. Identical to above except indices in Ps, and the subsets below requiring status to be responder/nonresponder
    # corona block
    subdatac <- subset (melt_data, variable == aVariable & corona.dash == "c" & Status135 == "R")
    sub1c <- subset (subdatac, Case == 1)
    sub2c <- subset (subdatac, Case == 2)
    sub3c <- subset (subdatac, Case == 3)
    #Select the data present in every timepoint
    subdatac <- subset (subdatac, Patient %in% intersect (unique (sub1c$Patient), intersect (unique (sub2c$Patient), unique (sub3c$Patient))))
    sub1c <- subset (subdatac, Case == 1)
    sub2c <- subset (subdatac, Case == 2)
    sub3c <- subset (subdatac, Case == 3)
    subdatac$value <- as.numeric (subdatac$value)
    sub1c$value <- as.numeric (sub1c$value)
    sub2c$value <- as.numeric (sub2c$value)
    sub3c$value <- as.numeric (sub3c$value)

    if (! (aVariable %in% varsNotInCase2)) { # rule out the variables not defined in case 2
      # cliff's delta
      d12 <- as.numeric (orddom::orddom (x = (sub1c [! is.na (sub1c$value) & ! is.na (sub2c$value), ])$value, y = (sub2c [! is.na (sub1c$value) & ! is.na (sub2c$value), ])$value, paired = T)[11,1])
      # wilcoxon
      pw12 <- wilcox.test (sub1c$value, sub2c$value, paired = T)$p.value

      m1 <- lm (data = subset (subdatac, Case == 1 | Case == 2), value ~ Patient + age + sex + ACE_INHIBITOR + ANTIDEPRESSANT + ANTIDIABETIC + ANTI_PLATELET + AT2_BLOCKER + BETA_BLOCKER + CALCIUMANTAGONIST + CORTICOSTEROID + DIURETIC + IMIDAZOLE_AGONIST + METFORMIN + PPI + STATIN + THYROID + VITAMIN_D + XANTHAN_OXIDASE_INHIBITOR)
      m2 <- lm (data = subset (subdatac, Case == 1 | Case == 2), value ~ Patient + Case + age + sex + ACE_INHIBITOR + ANTIDEPRESSANT + ANTIDIABETIC + ANTI_PLATELET + AT2_BLOCKER + BETA_BLOCKER + CALCIUMANTAGONIST + CORTICOSTEROID + DIURETIC + IMIDAZOLE_AGONIST + METFORMIN + PPI + STATIN + THYROID + VITAMIN_D + XANTHAN_OXIDASE_INHIBITOR)
      p12 <- lmtest::lrtest (m1, m2)$"Pr(>Chisq)" [2]

    }

    # cliff's delta
    d13 <- as.numeric (orddom::orddom (x = (sub1c [! is.na (sub1c$value) & ! is.na (sub3c$value), ])$value, y = (sub3c [! is.na (sub1c$value) & ! is.na (sub3c$value), ])$value, paired = T)[11,1])
    # wilcoxon
    pw13 <- wilcox.test (sub1c$value, sub3c$value, paired = T)$p.value

    m1 <- lm (data = subset (subdatac, Case == 1 | Case == 3), value ~ Patient + age + sex + ACE_INHIBITOR + ANTIDEPRESSANT + ANTIDIABETIC + ANTI_PLATELET + AT2_BLOCKER + BETA_BLOCKER + CALCIUMANTAGONIST + CORTICOSTEROID + DIURETIC + IMIDAZOLE_AGONIST + METFORMIN + PPI + STATIN + THYROID + VITAMIN_D + XANTHAN_OXIDASE_INHIBITOR)
    m2 <- lm (data = subset (subdatac, Case == 1 | Case == 3), value ~ Patient + Case + age + sex + ACE_INHIBITOR + ANTIDEPRESSANT + ANTIDIABETIC + ANTI_PLATELET + AT2_BLOCKER + BETA_BLOCKER + CALCIUMANTAGONIST + CORTICOSTEROID + DIURETIC + IMIDAZOLE_AGONIST + METFORMIN + PPI + STATIN + THYROID + VITAMIN_D + XANTHAN_OXIDASE_INHIBITOR)
    p13 <- lmtest::lrtest (m1, m2)$"Pr(>Chisq)" [2]


    if (! (aVariable %in% varsNotInCase2)) { # rule out the variables not defined in case 2
      # cliff's delta
      d23 <- as.numeric (orddom::orddom (x = (sub2c [! is.na (sub2c$value) & ! is.na (sub3c$value), ])$value, y = (sub3c [! is.na (sub2c$value) & ! is.na (sub3c$value), ])$value, paired = T)[11,1])
      # wilcoxon
      pw23 <- wilcox.test (sub2c$value, sub3c$value, paired = T)$p.value

      m1 <- lm (data = subset (subdatac, Case == 2 | Case == 3), value ~ Patient + age + sex + ACE_INHIBITOR + ANTIDEPRESSANT + ANTIDIABETIC + ANTI_PLATELET + AT2_BLOCKER + BETA_BLOCKER + CALCIUMANTAGONIST + CORTICOSTEROID + DIURETIC + IMIDAZOLE_AGONIST + METFORMIN + PPI + STATIN + THYROID + VITAMIN_D + XANTHAN_OXIDASE_INHIBITOR)
      m2 <- lm (data = subset (subdatac, Case == 2 | Case == 3), value ~ Patient + Case + age + sex + ACE_INHIBITOR + ANTIDEPRESSANT + ANTIDIABETIC + ANTI_PLATELET + AT2_BLOCKER + BETA_BLOCKER + CALCIUMANTAGONIST + CORTICOSTEROID + DIURETIC + IMIDAZOLE_AGONIST + METFORMIN + PPI + STATIN + THYROID + VITAMIN_D + XANTHAN_OXIDASE_INHIBITOR)
      p23 <- lmtest::lrtest (m1, m2)$"Pr(>Chisq)" [2]

    }

    Ps [i, 11] <- p12
    Ps [i, 12] <- pw12
    Ps [i, 13] <- d12

    Ps [i, 14] <- p13
    Ps [i, 15] <- pw13
    Ps [i, 16] <- d13

    Ps [i, 17] <- p23
    Ps [i, 18] <- pw23
    Ps [i, 19] <- d23
    print("passed block 2")

    #### Responders block done

    #### Then repeat all for non-responders,identical to above except indices in Ps, and the subsets below requiring status to be responder/nonresponder
    # corona block
    subdatac <- subset (melt_data, variable == aVariable & corona.dash == "c" & Status135 == "non-R")
    sub1c <- subset (subdatac, Case == 1)
    sub2c <- subset (subdatac, Case == 2)
    sub3c <- subset (subdatac, Case == 3)
    subdatac <- subset (subdatac, Patient %in% intersect (unique (sub1c$Patient), intersect (unique (sub2c$Patient), unique (sub3c$Patient))))
    sub1c <- subset (subdatac, Case == 1)
    sub2c <- subset (subdatac, Case == 2)
    sub3c <- subset (subdatac, Case == 3)

    subdatac$value <- as.numeric (subdatac$value)
    sub1c$value <- as.numeric (sub1c$value)
    sub2c$value <- as.numeric (sub2c$value)
    sub3c$value <- as.numeric (sub3c$value)


    if (! (aVariable %in% varsNotInCase2)) { # rule out the variables not defined in case 2

      # cliff's delta
      d12 <- as.numeric (orddom::orddom (x = (sub1c [! is.na (sub1c$value) & ! is.na (sub2c$value), ])$value, y = (sub2c [! is.na (sub1c$value) & ! is.na (sub2c$value), ])$value, paired = T)[11,1])
      # wilcoxon
      pw12 <- wilcox.test (sub1c$value, sub2c$value, paired = T)$p.value

      m1 <-lm (data = subset (subdatac, Case == 1 | Case == 2), value ~ Patient + age + sex + ACE_INHIBITOR + ANTIDEPRESSANT + ANTIDIABETIC + ANTI_PLATELET + AT2_BLOCKER + BETA_BLOCKER + CALCIUMANTAGONIST + CORTICOSTEROID + DIURETIC + IMIDAZOLE_AGONIST + METFORMIN + PPI + STATIN + THYROID + VITAMIN_D + XANTHAN_OXIDASE_INHIBITOR)
      m2 <- lm (data = subset (subdatac, Case == 1 | Case == 2), value ~ Patient + Case + age + sex + ACE_INHIBITOR + ANTIDEPRESSANT + ANTIDIABETIC + ANTI_PLATELET + AT2_BLOCKER + BETA_BLOCKER + CALCIUMANTAGONIST + CORTICOSTEROID + DIURETIC + IMIDAZOLE_AGONIST + METFORMIN + PPI + STATIN + THYROID + VITAMIN_D + XANTHAN_OXIDASE_INHIBITOR)
      p12 <- lmtest::lrtest (m1, m2)$"Pr(>Chisq)" [2]

    }

    # cliff's delta
    d13 <- as.numeric (orddom::orddom (x = (sub1c [! is.na (sub1c$value) & ! is.na (sub3c$value), ])$value, y = (sub3c [! is.na (sub1c$value) & ! is.na (sub3c$value), ])$value, paired = T)[11,1])
    # wilcoxon
    pw13 <- wilcox.test (sub1c$value, sub3c$value, paired = T)$p.value

    m1 <- lm (data = subset (subdatac, Case == 1 | Case == 3), value ~ Patient + age + sex + ACE_INHIBITOR + ANTIDEPRESSANT + ANTIDIABETIC + ANTI_PLATELET + AT2_BLOCKER + BETA_BLOCKER + CALCIUMANTAGONIST + CORTICOSTEROID + DIURETIC + IMIDAZOLE_AGONIST + METFORMIN + PPI + STATIN + THYROID + VITAMIN_D + XANTHAN_OXIDASE_INHIBITOR)
    m2 <- lm (data = subset (subdatac, Case == 1 | Case == 3), value ~ Patient + Case + age + sex + ACE_INHIBITOR + ANTIDEPRESSANT + ANTIDIABETIC + ANTI_PLATELET + AT2_BLOCKER + BETA_BLOCKER + CALCIUMANTAGONIST + CORTICOSTEROID + DIURETIC + IMIDAZOLE_AGONIST + METFORMIN + PPI + STATIN + THYROID + VITAMIN_D + XANTHAN_OXIDASE_INHIBITOR)
    p13 <- lmtest::lrtest (m1, m2)$"Pr(>Chisq)" [2]


    if (! (aVariable %in% varsNotInCase2)) { # rule out the variables not defined in case 2

      # cliff's delta
      d23 <- as.numeric (orddom::orddom (x = (sub2c [! is.na (sub2c$value) & ! is.na (sub3c$value), ])$value, y = (sub3c [! is.na (sub2c$value) & ! is.na (sub3c$value), ])$value, paired = T)[11,1])
      # wilcoxon
      pw23 <- wilcox.test (sub2c$value, sub3c$value, paired = T)$p.value

      m1 <- lm (data = subset (subdatac, Case == 2 | Case == 3), value ~ Patient + age + sex + ACE_INHIBITOR + ANTIDEPRESSANT + ANTIDIABETIC + ANTI_PLATELET + AT2_BLOCKER + BETA_BLOCKER + CALCIUMANTAGONIST + CORTICOSTEROID + DIURETIC + IMIDAZOLE_AGONIST + METFORMIN + PPI + STATIN + THYROID + VITAMIN_D + XANTHAN_OXIDASE_INHIBITOR)
      m2 <- lm (data = subset (subdatac, Case == 2 | Case == 3), value ~ Patient + Case + age + sex + ACE_INHIBITOR + ANTIDEPRESSANT + ANTIDIABETIC + ANTI_PLATELET + AT2_BLOCKER + BETA_BLOCKER + CALCIUMANTAGONIST + CORTICOSTEROID + DIURETIC + IMIDAZOLE_AGONIST + METFORMIN + PPI + STATIN + THYROID + VITAMIN_D + XANTHAN_OXIDASE_INHIBITOR)
      p23 <- lmtest::lrtest (m1, m2)$"Pr(>Chisq)" [2]

    }

    Ps [i, 20] <- p12
    Ps [i, 21] <- pw12
    Ps [i, 22] <- d12

    Ps [i, 23] <- p13
    Ps [i, 24] <- pw13
    Ps [i, 25] <- d13

    Ps [i, 26] <- p23
    Ps [i, 27] <- pw23
    Ps [i, 28] <- d23
    print("passed block 3")

    #### Non-responders block done

    #### Test between corona responder and non-responder in case 1
    subdatac <- subset (melt_data, variable == aVariable & corona.dash == "c" & Case == 1 & (Status135=="R"|Status135=="non-R") )
    subRc <- subset (subdatac, Status135 == "R")
    subNRc <- subset (subdatac, Status135 == "non-R")

    subdatac$value <- as.numeric (subdatac$value)
    subRc$value <- as.numeric (subRc$value)
    subNRc$value <- as.numeric (subNRc$value)

    #cliffs delta
    cdRNR1 <- as.numeric (orddom::orddom (x = (subRc [! is.na (subRc$value), ])$value, y = (subNRc [! is.na (subNRc$value), ])$value, paired = F)[13,1])
    # wilcox
    cwpRNR1 <- wilcox.test (subRc$value, subNRc$value, paired = F)$p.value

    # nested model comparison
    m1 <- lm (data = subdatac, value ~ age + sex + ACE_INHIBITOR + ANTIDEPRESSANT + ANTIDIABETIC + ANTI_PLATELET + AT2_BLOCKER + BETA_BLOCKER + CALCIUMANTAGONIST + CORTICOSTEROID + DIURETIC + IMIDAZOLE_AGONIST + METFORMIN + PPI + STATIN + THYROID + VITAMIN_D + XANTHAN_OXIDASE_INHIBITOR)
    # model 2 = model 1 + case
    m2 <- lm (data = subdatac, value ~ Status135 + age + sex + ACE_INHIBITOR + ANTIDEPRESSANT + ANTIDIABETIC + ANTI_PLATELET + AT2_BLOCKER + BETA_BLOCKER + CALCIUMANTAGONIST + CORTICOSTEROID + DIURETIC + IMIDAZOLE_AGONIST + METFORMIN + PPI + STATIN + THYROID + VITAMIN_D + XANTHAN_OXIDASE_INHIBITOR)
    #likelyhood compare models
    cpRNR1 <- lmtest::lrtest (m1, m2)$"Pr(>Chisq)" [2]

    Ps [i, 29] <- cpRNR1
    Ps [i, 30] <- cwpRNR1
    Ps [i, 31] <- cdRNR1

    # Kruskal Wallis test for the model test on corona patients in all cases
    subdatac <- subset (melt_data, variable == aVariable & corona.dash == "c")

    subdatac$value <- as.numeric (subdatac$value)

    m1 <- lm (data = subset (subdatac), rank (value) ~ Patient + age + sex + ACE_INHIBITOR + ANTIDEPRESSANT + ANTIDIABETIC + ANTI_PLATELET + AT2_BLOCKER + BETA_BLOCKER + CALCIUMANTAGONIST + CORTICOSTEROID + DIURETIC + IMIDAZOLE_AGONIST + METFORMIN + PPI + STATIN + THYROID + VITAMIN_D + XANTHAN_OXIDASE_INHIBITOR)
    m2 <- lm (data = subset (subdatac), rank (value) ~ Patient + Case + age + sex + ACE_INHIBITOR + ANTIDEPRESSANT + ANTIDIABETIC + ANTI_PLATELET + AT2_BLOCKER + BETA_BLOCKER + CALCIUMANTAGONIST + CORTICOSTEROID + DIURETIC + IMIDAZOLE_AGONIST + METFORMIN + PPI + STATIN + THYROID + VITAMIN_D + XANTHAN_OXIDASE_INHIBITOR)
    pModel <- lmtest::lrtest (m1, m2)$"Pr(>Chisq)" [2]

    Ps [i, 32] <- pModel

    Ps [i, 33] <- p.adjust (method = "fdr", c (Ps [i, 3], Ps [i, 6], Ps [i, 9])) [1]
    Ps [i, 34] <- p.adjust (method = "fdr", c (Ps [i, 3], Ps [i, 6], Ps [i, 9])) [2]
    Ps [i, 35] <- p.adjust (method = "fdr", c (Ps [i, 3], Ps [i, 6], Ps [i, 9])) [3]

    # Kruskal Wallis test for the model test on corona patients in responders
    subdatac <- subset (melt_data, variable == aVariable & corona.dash == "c" & Status135 == "R")

    subdatac$value <- as.numeric (subdatac$value)

    m1 <- lm (data = subset (subdatac), rank (value) ~ Patient + age + sex + ACE_INHIBITOR + ANTIDEPRESSANT + ANTIDIABETIC + ANTI_PLATELET + AT2_BLOCKER + BETA_BLOCKER + CALCIUMANTAGONIST + CORTICOSTEROID + DIURETIC + IMIDAZOLE_AGONIST + METFORMIN + PPI + STATIN + THYROID + VITAMIN_D + XANTHAN_OXIDASE_INHIBITOR)
    m2 <- lm (data = subset (subdatac), rank (value) ~ Patient + Case + age + sex + ACE_INHIBITOR + ANTIDEPRESSANT + ANTIDIABETIC + ANTI_PLATELET + AT2_BLOCKER + BETA_BLOCKER + CALCIUMANTAGONIST + CORTICOSTEROID + DIURETIC + IMIDAZOLE_AGONIST + METFORMIN + PPI + STATIN + THYROID + VITAMIN_D + XANTHAN_OXIDASE_INHIBITOR)
    pModel <- lmtest::lrtest (m1, m2)$"Pr(>Chisq)" [2]

    Ps [i, 36] <- pModel

    # Adjust p value with Benjamini Hochberg method to get q value
    Ps [i, 37] <- p.adjust (method = "fdr", c (Ps [i, 12], Ps [i, 15], Ps [i, 18])) [1]
    Ps [i, 38] <- p.adjust (method = "fdr", c (Ps [i, 12], Ps [i, 15], Ps [i, 18])) [2]
    Ps [i, 39] <- p.adjust (method = "fdr", c (Ps [i, 12], Ps [i, 15], Ps [i, 18])) [3]


    #### Test between corona responder and non-responder in case 2
    if (! (aVariable %in% varsNotInCase2)) { # rule out the variables not defined in case 2
      subdatac <- subset (melt_data, variable == aVariable & corona.dash == "c" & Case == 2 & (Status135=="R"|Status135=="non-R") )
      subRc <- subset (subdatac, Status135 == "R")
      subNRc <- subset (subdatac, Status135 == "non-R")

      subdatac$value <- as.numeric (subdatac$value)
      subRc$value <- as.numeric (subRc$value)
      subNRc$value <- as.numeric (subNRc$value)


      #cliffs delta
      #	cdRNR1 <- as.numeric (orddom (x = (subRc [! is.na (subRc$value) & ! is.na (subNRc$value), ])$value, y = (subNRc [! is.na (subNRc$value) & ! is.na (subRc$value), ])$value, paired = F)[13,1])
      cdRNR2 <- as.numeric (orddom::orddom (x = (subRc [! is.na (subRc$value), ])$value, y = (subNRc [! is.na (subNRc$value), ])$value, paired = F)[13,1])
      ## SF 20181019 - corrected the above to match it is not paired

      #wilcox
      cwpRNR2 <- wilcox.test (subRc$value, subNRc$value, paired = F)$p.value

      #nested model comparison
      m1 <- lm (data = subdatac, value ~ age + sex + ACE_INHIBITOR + ANTIDEPRESSANT + ANTIDIABETIC + ANTI_PLATELET + AT2_BLOCKER + BETA_BLOCKER + CALCIUMANTAGONIST + CORTICOSTEROID + DIURETIC + IMIDAZOLE_AGONIST + METFORMIN + PPI + STATIN + THYROID + VITAMIN_D + XANTHAN_OXIDASE_INHIBITOR)
      # model 2 = model 1 + case

      m2 <- lm (data = subdatac, value ~ Status135 + age + sex + ACE_INHIBITOR + ANTIDEPRESSANT + ANTIDIABETIC + ANTI_PLATELET + AT2_BLOCKER + BETA_BLOCKER + CALCIUMANTAGONIST + CORTICOSTEROID + DIURETIC + IMIDAZOLE_AGONIST + METFORMIN + PPI + STATIN + THYROID + VITAMIN_D + XANTHAN_OXIDASE_INHIBITOR)
      #likelyhood compare models
      cpRNR2 <- lmtest::lrtest (m1, m2)$"Pr(>Chisq)" [2]

      Ps [i, 40] <- cpRNR2
      Ps [i, 41] <- cwpRNR2
      Ps [i, 42] <- cdRNR2


    }

    #### Test between corona responder and non-responder in case 3
    subdatac <- subset (melt_data, variable == aVariable & corona.dash == "c" & Case == 3 & (Status135=="R"|Status135=="non-R") )
    subRc <- subset (subdatac, Status135 == "R")
    subNRc <- subset (subdatac, Status135 == "non-R")

    subdatac$value <- as.numeric (subdatac$value)
    subRc$value <- as.numeric (subRc$value)
    subNRc$value <- as.numeric (subNRc$value)

    cdRNR3 <- as.numeric (orddom::orddom (x = (subRc [! is.na (subRc$value), ])$value, y = (subNRc [! is.na (subNRc$value), ])$value, paired = F)[13,1])
    cwpRNR3 <- wilcox.test (subRc$value, subNRc$value, paired = F)$p.value

    # nested model comparison, model 2 = model 1 + case
    m1 <- lm (data = subdatac, value ~ age + sex + ACE_INHIBITOR + ANTIDEPRESSANT + ANTIDIABETIC + ANTI_PLATELET + AT2_BLOCKER + BETA_BLOCKER + CALCIUMANTAGONIST + CORTICOSTEROID + DIURETIC + IMIDAZOLE_AGONIST + METFORMIN + PPI + STATIN + THYROID + VITAMIN_D + XANTHAN_OXIDASE_INHIBITOR)

    m2 <- lm (data = subdatac, value ~ Status135 + age + sex + ACE_INHIBITOR + ANTIDEPRESSANT + ANTIDIABETIC + ANTI_PLATELET + AT2_BLOCKER + BETA_BLOCKER + CALCIUMANTAGONIST + CORTICOSTEROID + DIURETIC + IMIDAZOLE_AGONIST + METFORMIN + PPI + STATIN + THYROID + VITAMIN_D + XANTHAN_OXIDASE_INHIBITOR)
    #likelyhood compare models
    cpRNR3 <- lmtest::lrtest (m1, m2)$"Pr(>Chisq)" [2]

    Ps [i, 43] <- cpRNR3
    Ps [i, 44] <- cwpRNR3
    Ps [i, 45] <- cdRNR3
    Ps$qModel <- p.adjust (method = "fdr", Ps$pModel)
  }
  print("passed for loop")
  write.table (file = paste0(p_table_name, ".Ps"), Ps, sep = "\t", quote = F)
  print("Superb! Your Ps table is now in your current directory.")
}
