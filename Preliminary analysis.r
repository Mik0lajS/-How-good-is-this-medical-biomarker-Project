library(ADNIMERGE)
library(pracma)
library(ggplot2)
library(mgcv)
library(scam)

adnimerge_test <- merge(adnimerge, baipetnmrcfdg[, c("RID", "ORIGPROT", "VISCODE", "HCI")], 
                        by = c("RID", "ORIGPROT", "VISCODE"), all = TRUE)

# Some in baipetnmrcfdg but not in adnimerge: delete them (because most information, like AGE or GENDER is missing)
adnimerge <- adnimerge_test[!is.na(adnimerge_test$PTGENDER),]

# DX: diagnosis at every visit (VISCODE)
# DX.bl: diagnosis at baseline

table(adnimerge[adnimerge$VISCODE == "bl", "DX"])
table(adnimerge[adnimerge$VISCODE == "bl", "DX.bl"])

# CN (cognitive normal) in DX == CN + SMC (Subjective memory complaints) en DX.bl
# MCI (mild cognitive impairment) in DX == EMCI (early mild cognitive impairment) + LMCI (late mild cognitive impairment) in DX.bl
# Dementia en DX == AD (Alzheimer’s disease) en DX.bl
# At baseline in DX (VISCODE == "bl") there are more NAs

# Data at baseline (and only some of the variables)
data_baseline <- adnimerge[adnimerge$VISCODE == "bl", c("RID", "ORIGPROT", "COLPROT", "DX", "DX.bl", "AGE", "PTGENDER", "APOE4", "ABETA", "TAU","PTAU", "HCI")]
data_baseline$PTGENDER <- factor(data_baseline$PTGENDER)
names(data_baseline)
#summary(data_baseline)


# Biomarkers are codified as "character" (limit of detection????)
data_baseline$ABETA_num <- data_baseline$ABETA
data_baseline$ABETA_num[data_baseline$ABETA_num == ">1700"] <- 1700
data_baseline$ABETA_num[data_baseline$ABETA_num == "<200"] <- 200
data_baseline$ABETA_num <- as.numeric(data_baseline$ABETA_num)

data_baseline$TAU_num <- data_baseline$TAU
data_baseline$TAU_num[data_baseline$TAU_num == ">1300"] <- 1300
data_baseline$TAU_num[data_baseline$TAU_num == "<80"] <- 80
data_baseline$TAU_num <- as.numeric(data_baseline$TAU_num)

data_baseline$PTAU_num <- data_baseline$PTAU
data_baseline$PTAU_num[data_baseline$PTAU_num == ">120"] <- 120
data_baseline$PTAU_num[data_baseline$PTAU_num == "<8"] <- 8
data_baseline$PTAU_num <- as.numeric(data_baseline$PTAU_num)


# Create vectors of biomarker data for each category and biomarker

# ABETA
y1_ABETA <- data_baseline$ABETA_num[which(data_baseline$DX.bl == "CN" | data_baseline$DX.bl == "SMC")]
y1_ABETA <- y1_ABETA[!is.na(y1_ABETA)]

y2_ABETA <- data_baseline$ABETA_num[which(data_baseline$DX.bl == "EMCI" | data_baseline$DX.bl == "LMCI")]
y2_ABETA <- y2_ABETA[!is.na(y2_ABETA)]

y3_ABETA <- data_baseline$ABETA_num[which(data_baseline$DX.bl == "AD")]
y3_ABETA <- y3_ABETA[!is.na(y3_ABETA)]

# TAU
y1_TAU <- data_baseline$TAU_num[which(data_baseline$DX.bl == "CN" | data_baseline$DX.bl == "SMC")]
y1_TAU <- y1_TAU[!is.na(y1_TAU)]

y2_TAU <- data_baseline$TAU_num[which(data_baseline$DX.bl == "EMCI" | data_baseline$DX.bl == "LMCI")]
y2_TAU <- y2_TAU[!is.na(y2_TAU)]

y3_TAU <- data_baseline$TAU_num[which(data_baseline$DX.bl == "AD")]
y3_TAU <- y3_TAU[!is.na(y3_TAU)]

# PTAU
y1_PTAU <- data_baseline$PTAU_num[which(data_baseline$DX.bl == "CN" | data_baseline$DX.bl == "SMC")]
y1_PTAU <- y1_PTAU[!is.na(y1_PTAU)]

y2_PTAU <- data_baseline$PTAU_num[which(data_baseline$DX.bl == "EMCI" | data_baseline$DX.bl == "LMCI")]
y2_PTAU <- y2_PTAU[!is.na(y2_PTAU)]

y3_PTAU <- data_baseline$PTAU_num[which(data_baseline$DX.bl == "AD")]
y3_PTAU <- y3_PTAU[!is.na(y3_PTAU)]

# HCI
y1_HCI <- data_baseline$HCI[which(data_baseline$DX.bl == "CN" | data_baseline$DX.bl == "SMC")]
y1_HCI <- y1_HCI[!is.na(y1_HCI)]

y2_HCI <- data_baseline$HCI[which(data_baseline$DX.bl == "EMCI" | data_baseline$DX.bl == "LMCI")]
y2_HCI <- y2_HCI[!is.na(y2_HCI)]

y3_HCI <- data_baseline$HCI[which(data_baseline$DX.bl == "AD")]
y3_HCI <- y3_HCI[!is.na(y3_HCI)]


# Plot densities and CDFs

# ABETA
par(mfrow = c(1, 2))
plot(density(y1_ABETA), col = "blue2", xlab= "ABETA", ylab = "Density", ylim = c(0, 0.003), main = "")
lines(density(y2_ABETA), col = "orange")
lines(density(y3_ABETA), col ="red")
legend("topright", legend = c("CN + SMC", "EMCI + LMCI", "AD"),
       lty = 1, col = c("blue2", "orange", "red"), bty = "n")
plot(ecdf(y1_ABETA), col = "blue2", xlab= "ABETA", ylab = "CDF", main = "")
lines(ecdf(y2_ABETA), col = "orange")
lines(ecdf(y3_ABETA), col = "red")

# TAU
par(mfrow = c(1, 2))
plot(density(y1_TAU), col = "blue2", xlab= "TAU", ylab = "Density", main = "")
lines(density(y2_TAU), col = "orange")
lines(density(y3_TAU), col ="red")
legend("topright", legend = c("CN + SMC", "EMCI + LMCI", "AD"),
       lty = 1, col = c("blue2", "orange", "red"), bty = "n")
plot(ecdf(y1_TAU), col = "blue2", xlab= "TAU", ylab = "CDF", main = "")
lines(ecdf(y2_TAU), col = "orange")
lines(ecdf(y3_TAU), col = "red")

# PTAU
par(mfrow = c(1, 2))
plot(density(y1_PTAU), col = "blue2", xlab= "PTAU", ylab = "Density", main = "")
lines(density(y2_PTAU), col = "orange")
lines(density(y3_PTAU), col ="red")
legend("topright", legend = c("CN + SMC", "EMCI + LMCI", "AD"),
       lty = 1, col = c("blue2", "orange", "red"), bty = "n")
plot(ecdf(y1_PTAU), col = "blue2", xlab= "PTAU", ylab = "CDF", main = "")
lines(ecdf(y2_PTAU), col = "orange")
lines(ecdf(y3_PTAU), col = "red")

# HCI
par(mfrow = c(1, 2))
plot(density(y1_HCI), col = "blue2", xlab= "HCI", ylab = "Density", xlim = c(0, 60), main = "")
lines(density(y2_HCI), col = "orange")
lines(density(y3_HCI), col ="red")
legend("topright", legend = c("CN + SMC", "EMCI + LMCI", "AD"),
       lty = 1, col = c("blue2", "orange", "red"), bty = "n")
plot(ecdf(y1_HCI), col = "blue2", xlab= "HCI", ylab = "CDF", main = "")
lines(ecdf(y2_HCI), col = "orange")
lines(ecdf(y3_HCI), col = "red")




# change gender and APOE4 to factor
data_baseline$PTGENDER <- factor(data_baseline$PTGENDER)
data_baseline$APOE4 <- factor(data_baseline$APOE4, levels = c(0,1,2),
                              labels = c("No e4", "One e4", "Two e4"))
summary(data_baseline)
dim(data_baseline)
#help("complete.cases")

# drop NA from important columns
cols_to_check <- c("DX.bl", "AGE", "PTGENDER", "APOE4", "ABETA", "TAU", "PTAU", "HCI")
data_baseline_dropna <- data_baseline[complete.cases(data_baseline[, cols_to_check]), ]

dim(data_baseline_dropna)

# ABETA
y1_ABETA <- data_baseline_dropna$ABETA_num[which(data_baseline_dropna$DX.bl == "CN" | data_baseline_dropna$DX.bl == "SMC")]
y2_ABETA <- data_baseline_dropna$ABETA_num[which(data_baseline_dropna$DX.bl == "EMCI" | data_baseline_dropna$DX.bl == "LMCI")]
y3_ABETA <- data_baseline_dropna$ABETA_num[which(data_baseline_dropna$DX.bl == "AD")]

# TAU
y1_TAU <- data_baseline_dropna$TAU_num[which(data_baseline_dropna$DX.bl == "CN" | data_baseline_dropna$DX.bl == "SMC")]
y2_TAU <- data_baseline_dropna$TAU_num[which(data_baseline_dropna$DX.bl == "EMCI" | data_baseline_dropna$DX.bl == "LMCI")]
y3_TAU <- data_baseline_dropna$TAU_num[which(data_baseline_dropna$DX.bl == "AD")]

# PTAU
y1_PTAU <- data_baseline_dropna$PTAU_num[which(data_baseline_dropna$DX.bl == "CN" | data_baseline_dropna$DX.bl == "SMC")]
y2_PTAU <- data_baseline_dropna$PTAU_num[which(data_baseline_dropna$DX.bl == "EMCI" | data_baseline_dropna$DX.bl == "LMCI")]
y3_PTAU <- data_baseline_dropna$PTAU_num[which(data_baseline_dropna$DX.bl == "AD")]

# HCI
y1_HCI <- data_baseline_dropna$HCI[which(data_baseline_dropna$DX.bl == "CN" | data_baseline_dropna$DX.bl == "SMC")]
y2_HCI <- data_baseline_dropna$HCI[which(data_baseline_dropna$DX.bl == "EMCI" | data_baseline_dropna$DX.bl == "LMCI")]
y3_HCI <- data_baseline_dropna$HCI[which(data_baseline_dropna$DX.bl == "AD")]


# create new column with 1, 2 or 3 depending on cognitive status
data_baseline_dropna$class <- ifelse(
  data_baseline_dropna$DX.bl %in% c("CN", "SMC"), 1,
  ifelse(
    data_baseline_dropna$DX.bl %in% c("EMCI", "LMCI"), 2,
    3
  )
)

data_baseline_dropna$ABETA_neg <- -data_baseline_dropna$ABETA_num

# Obtaining TPF triples for ROC surface
obtain_TPF_triples <- function(y1, y2, y3){
  
  # obtain a grid of thresholds
  grid <- sort(unique(c(y1, y2, y3)))
  
  n <- length(grid)
  n_pairs <-  n * (n-1) / 2
  
  t1 <- t2 <- numeric(n_pairs)
  
  # find all viable pairs of thresholds
  k <- 1
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      t1[k] <- grid[i]
      t2[k] <- grid[j]
      k <- k + 1
    }
  }
  
  # define empirical cdfs
  F1<- ecdf(y1)
  F2 <- ecdf(y2)
  F3 <- ecdf(y3)
  
  # precompute ecdf values 
  TPF1_vals <- F1(grid)
  TPF2_vals <- F2(grid)
  TPF3_vals <- F3(grid)
  
  # vectors with indices of ecdf values 
  t1_idx <- match(t1, grid)
  t2_idx <- match(t2, grid)
  
  TPF1 <- TPF1_vals[t1_idx]
  TPF3 <- 1 - TPF3_vals[t2_idx]
  TPF2 <- TPF2_vals[t2_idx] - TPF2_vals[t1_idx]
  
  return(data.frame(
    t1 = t1,
    t2 = t2,
    TPF1 = TPF1,
    TPF2 = TPF2,
    TPF3 = TPF3
  ))
}


obtain_TPFs_functional <- function(y1, y2, y3, len_grid = 50){
# Returns: list with TPF1, TPF3 vectors and TPF2_mat matrix for plotting 
# the ROC surface with persp
  
  p1 <- p3 <- seq(0, 1, len = len_grid)
  
  F2e <- ecdf(y2)
  
  roce <- matrix(0, length(p1), length(p3))
  for (i in 1:length(p1)) {
    for (j in 1:length(p3)) {
      
      c1 = quantile(y1, p1[i], type = 1)
      c2 = quantile(y3, (1 - p3[j]), type = 1)
      
      if (c1 < c2) 
        roce[i, j] <- pmax(0, F2e(c2) - F2e(c1))
      else roce[i, j] <- 0
    }
  }
  return(list(TPF1 = p1, TPF3 = p3, TPF2_mat = roce))
}


# VUS estimate - Mann-Whitney statistic
Mann_Whitney <- function(y1, y2, y3){
  n1 <- length(y1)
  n2 <- length(y2)
  n3 <- length(y3)
  
  vus <- 0
  for(i in 1:n1){
    for(j in 1:n2){
      for(k in 1:n3){
        x1 <- y1[i]
        x2 <- y2[j]
        x3 <- y3[k]
        
        if(x1 < x2 && x2 < x3){
          vus <- vus + 1
        } else if(x1 < x2 && x2 == x3){
          vus <- vus + 0.5
        } else if(x1 == x2 && x2 < x3){
          vus <- vus + 0.5
        } else if(x1 == x2 && x2 == x3){
          vus <- vus + 1/6
        }
      }
    }
  }
  
  vus <- vus / (n1 * n2 * n3)
  return (vus)
}



# Numerical integration with trapz
VUS_trapz <- function(y1, y2, y3) {
  result <- obtain_TPFs_functional(y1, y2, y3)
  
  inner_integrals <- numeric(length(result$TPF1))
  
  # integrate along TPF3
  for (i in 1:length(result$TPF1)) {
    inner_integrals[i] <- trapz(result$TPF3, result$TPF2_mat[i, ])
  }
  
  # integrate along TPF1
  vus <- trapz(result$TPF1, inner_integrals)
  
  return(vus)
}



# ============================================
# ROC Surfaces plotting

par(mfrow = c(2, 2))

# VUS for each biomarker
print(Mann_Whitney_vectorised(y1_TAU, y2_TAU, y3_TAU))
print(Mann_Whitney_vectorised(y1_PTAU, y2_PTAU, y3_PTAU))
print(Mann_Whitney_vectorised(y1_HCI, y2_HCI, y3_HCI))
print(Mann_Whitney_vectorised(y3_ABETA, y2_ABETA, y1_ABETA))
print(Mann_Whitney_vectorised(-y1_ABETA, -y2_ABETA, -y3_ABETA))

print(VUS_trapz(y1_TAU, y2_TAU, y3_TAU))
print(VUS_trapz(y1_PTAU, y2_PTAU, y3_PTAU))
print(VUS_trapz(y1_HCI, y2_HCI, y3_HCI))
print(VUS_trapz(y3_ABETA, y2_ABETA, y1_ABETA))
print(VUS_trapz(-y1_ABETA, -y2_ABETA, -y3_ABETA))


par(mfrow = c(2, 2))

# HCI
vals_for_plot_HCI <- obtain_TPFs_functional(y1_HCI, y2_HCI, y3_HCI)
persp(vals_for_plot_HCI$TPF1, vals_for_plot_HCI$TPF3, vals_for_plot_HCI$TPF2_mat,
      phi = 30, theta = 60, xlab = "\n\n p1", ylab = "\n\n p3",
      zlab = "\n\n p2", ticktype = "simple", cex.lab = 1.4)
mtext("HCI, VUS = 0.5341", side = 3, line = 1, adj = 0, cex = 1.5)  

# TAU
vals_for_plot_TAU <- obtain_TPFs_functional(y1_TAU, y2_TAU, y3_TAU)
persp(vals_for_plot_TAU$TPF1, vals_for_plot_TAU$TPF3, vals_for_plot_TAU$TPF2_mat,
      phi = 30, theta = 60, xlab = "\n\n p1", ylab = "\n\n p3",
      zlab = "\n\n p2", ticktype = "simple", cex.lab = 1.4)
mtext("C Tau, VUS = 0.3418", side = 3, line = 1, adj = 0, cex = 1.5)

# PTAU
vals_for_plot_PTAU <- obtain_TPFs_functional(y1_PTAU, y2_PTAU, y3_PTAU)
persp(vals_for_plot_PTAU$TPF1, vals_for_plot_PTAU$TPF3, vals_for_plot_PTAU$TPF2_mat,
      phi = 30, theta = 60, xlab = "\n\n p1", ylab = "\n\n p3",
      zlab = "\n\n p2", ticktype = "simple", cex.lab = 1.4)
mtext("CSF pTau, VUS=0.3506", side = 3, line = 1, adj = 0, cex = 1.5)

# ABETA
vals_for_plot_ABETA <- obtain_TPFs_functional(-y1_ABETA, -y2_ABETA, -y3_ABETA)
persp(vals_for_plot_ABETA$TPF1, vals_for_plot_ABETA$TPF3, vals_for_plot_ABETA$TPF2_mat,
      phi = 30, theta = 60, xlab = "\n\n p1", ylab = "\n\n p3",
      zlab = "\n\n p2", ticktype = "simple", cex.lab = 1.4)
mtext("CSF Abeta, VUS= 0.3911", side = 3, line = 1, adj = 0, cex = 1.5)

# =============================================
# GAMs 
# =============================================

# modelling biomarker

# function implementing the 2 step GAM estimation procedure
GAM_estimation <- function(biomarker, main_df = data_baseline_dropna){
  
  fmla <- as.formula(paste(biomarker, "~ s(AGE, bs = 'ps', k = 10) + PTGENDER + APOE4"))
  
  mod1 <- gam(fmla, data = subset(main_df, class == 1), method = "GCV.Cp")
  mod2 <- gam(fmla, data = subset(main_df, class == 2), method = "GCV.Cp")
  mod3 <- gam(fmla, data = subset(main_df, class == 3), method = "GCV.Cp")
  
  # Modelling variance
  
  # add column storing log((yji - mji_pred)^2) to the df
  df <- main_df
  df$log_sq_residual <- NA
  
  idx1 <- df$class == 1
  df$log_sq_residual[idx1] <- log((df[[biomarker]][idx1] - predict(mod1, df[idx1, ]))^2)
  
  idx2 <- df$class == 2
  df$log_sq_residual[idx2] <- log((df[[biomarker]][idx2] - predict(mod2, df[idx2, ]))^2)
  
  idx3 <- df$class == 3
  df$log_sq_residual[idx3] <- log((df[[biomarker]][idx3] - predict(mod3, df[idx3, ]))^2)
  
  logvar_mod1 <- gam(log_sq_residual ~ s(AGE, bs = "ps", k = 10) + PTGENDER + APOE4, method = "GCV.Cp",
                     data = subset(df, class == 1))
  logvar_mod2 <- gam(log_sq_residual ~ s(AGE, bs = "ps", k = 10) + PTGENDER + APOE4, method = "GCV.Cp", 
                     data = subset(df, class == 2))
  logvar_mod3 <- gam(log_sq_residual ~ s(AGE, bs = "ps", k = 10) + PTGENDER + APOE4, method = "GCV.Cp", 
                     data = subset(df, class == 3))
  
  logvar_preds_1 <- as.numeric(predict(logvar_mod1, df[idx1, ]))
  logvar_preds_2 <- as.numeric(predict(logvar_mod2, df[idx2, ]))
  logvar_preds_3 <- as.numeric(predict(logvar_mod3, df[idx3, ]))
  
  theta_1 <- sum(exp(df$log_sq_residual[idx1]) * exp(logvar_preds_1)) / sum(((exp(logvar_preds_1))^2))
  theta_2 <- sum(exp(df$log_sq_residual[idx2]) * exp(logvar_preds_2)) / sum(((exp(logvar_preds_2))^2))
  theta_3 <- sum(exp(df$log_sq_residual[idx3]) * exp(logvar_preds_3)) / sum(((exp(logvar_preds_3))^2))
  
  
  SD_mod1 <- function(x){
    return(sqrt(theta_1 * exp(as.numeric(predict(logvar_mod1, x)))))
  }
  SD_mod2 <- function(x){
    return(sqrt(theta_2 * exp(as.numeric(predict(logvar_mod2, x)))))
  }
  SD_mod3 <- function(x){
    return(sqrt(theta_3 * exp(as.numeric(predict(logvar_mod3, x)))))
  }
  
  return(list(
    mod1 = mod1, mod2 = mod2, mod3 = mod3,
    logvar_mod1 = logvar_mod1, logvar_mod2 = logvar_mod2, logvar_mod3 = logvar_mod3,
    theta_1 = theta_1, theta_2 = theta_2, theta_3 = theta_3,
    SD_mod1 = SD_mod1, SD_mod2 = SD_mod2, SD_mod3 = SD_mod3
  ))
}

GAM_models_data_HCI <- GAM_estimation("HCI")

# ecdfs and other params
get_parameters_ecdfs_GAM <- function(mod1, mod2, mod3, SD_mod1, SD_mod2, SD_mod3) {
  
  idx1 <- data_baseline_dropna$class == 1
  idx2 <- data_baseline_dropna$class == 2
  idx3 <- data_baseline_dropna$class == 3
  
  SD_preds1 <- SD_mod1(data_baseline_dropna[idx1, ])
  SD_preds2 <- SD_mod2(data_baseline_dropna[idx2, ])
  SD_preds3 <- SD_mod3(data_baseline_dropna[idx3, ])
  
  # ecdfs 
  std_res1 <- residuals(mod1) / SD_preds1
  std_res2 <- residuals(mod2) / SD_preds2 
  std_res3 <- residuals(mod3) / SD_preds3
  
  ecdf_std_res1 <- ecdf(std_res1)
  ecdf_std_res2 <- ecdf(std_res2) 
  ecdf_std_res3 <- ecdf(std_res3)
  
  return(list(
    std_res1 = std_res1, std_res2 = std_res2, std_res3 = std_res3,
    ecdf_std_res1 = ecdf_std_res1, ecdf_std_res2 = ecdf_std_res2, ecdf_std_res3 = ecdf_std_res3
  ))
}



# modelling biomarker

# function implementing the 2 step GAM estimation procedure
GAM_estimation <- function(biomarker, main_df = data_baseline_dropna){
  
  fmla <- as.formula(paste(biomarker, "~ s(AGE, bs = 'ps', k = 10) + PTGENDER + APOE4"))
  
  mod1 <- gam(fmla, data = subset(main_df, class == 1), method = "GCV.Cp")
  mod2 <- gam(fmla, data = subset(main_df, class == 2), method = "GCV.Cp")
  mod3 <- gam(fmla, data = subset(main_df, class == 3), method = "GCV.Cp")
  
  # Modelling variance
  
  # add column storing log((yji - mji_pred)^2) to the df
  df <- main_df
  df$log_sq_residual <- NA
  
  idx1 <- df$class == 1
  df$log_sq_residual[idx1] <- log((df[[biomarker]][idx1] - predict(mod1, df[idx1, ]))^2)
  
  idx2 <- df$class == 2
  df$log_sq_residual[idx2] <- log((df[[biomarker]][idx2] - predict(mod2, df[idx2, ]))^2)
  
  idx3 <- df$class == 3
  df$log_sq_residual[idx3] <- log((df[[biomarker]][idx3] - predict(mod3, df[idx3, ]))^2)
  
  logvar_mod1 <- gam(log_sq_residual ~ s(AGE, bs = "ps", k = 10) + PTGENDER + APOE4, method = "GCV.Cp",
                     data = subset(df, class == 1))
  logvar_mod2 <- gam(log_sq_residual ~ s(AGE, bs = "ps", k = 10) + PTGENDER + APOE4, method = "GCV.Cp", 
                     data = subset(df, class == 2))
  logvar_mod3 <- gam(log_sq_residual ~ s(AGE, bs = "ps", k = 10) + PTGENDER + APOE4, method = "GCV.Cp", 
                     data = subset(df, class == 3))
  
  logvar_preds_1 <- as.numeric(predict(logvar_mod1, df[idx1, ]))
  logvar_preds_2 <- as.numeric(predict(logvar_mod2, df[idx2, ]))
  logvar_preds_3 <- as.numeric(predict(logvar_mod3, df[idx3, ]))
  
  theta_1 <- sum(exp(df$log_sq_residual[idx1]) * exp(logvar_preds_1)) / sum(((exp(logvar_preds_1))^2))
  theta_2 <- sum(exp(df$log_sq_residual[idx2]) * exp(logvar_preds_2)) / sum(((exp(logvar_preds_2))^2))
  theta_3 <- sum(exp(df$log_sq_residual[idx3]) * exp(logvar_preds_3)) / sum(((exp(logvar_preds_3))^2))
  
  
  SD_mod1 <- function(x){
    return(sqrt(theta_1 * exp(as.numeric(predict(logvar_mod1, x)))))
  }
  SD_mod2 <- function(x){
    return(sqrt(theta_2 * exp(as.numeric(predict(logvar_mod2, x)))))
  }
  SD_mod3 <- function(x){
    return(sqrt(theta_3 * exp(as.numeric(predict(logvar_mod3, x)))))
  }
  
  return(list(
    mod1 = mod1, mod2 = mod2, mod3 = mod3,
    logvar_mod1 = logvar_mod1, logvar_mod2 = logvar_mod2, logvar_mod3 = logvar_mod3,
    theta_1 = theta_1, theta_2 = theta_2, theta_3 = theta_3,
    SD_mod1 = SD_mod1, SD_mod2 = SD_mod2, SD_mod3 = SD_mod3
  ))
}


# returns a function approximating scam smoothed ecdfs and its inverse by interpolation
fast_approx_scam_cdf <- function(model, x_data, n_grid = 400) {
  
  x_grid <- seq(min(x_data), max(x_data), length.out = n_grid)
  y_grid <- as.numeric(predict(model, newdata = data.frame(x1 = x_grid, x2 = x_grid, x3 = x_grid)))
  
  y_grid <- pmin(pmax(y_grid, 0), 1)
  
  F_fast <- approxfun(x_grid, y_grid, rule = 2)
  
  
  return(F_fast = F_fast)
}

# ecdfs and other params
get_parameters_ecdfs_GAM <- function(mod1, mod2, mod3, SD_mod1, SD_mod2, SD_mod3) {
  
  idx1 <- data_baseline_dropna$class == 1
  idx2 <- data_baseline_dropna$class == 2
  idx3 <- data_baseline_dropna$class == 3
  
  SD_preds1 <- SD_mod1(data_baseline_dropna[idx1, ])
  SD_preds2 <- SD_mod2(data_baseline_dropna[idx2, ])
  SD_preds3 <- SD_mod3(data_baseline_dropna[idx3, ])
  
  # ecdfs 
  std_res1 <- residuals(mod1) / SD_preds1
  std_res2 <- residuals(mod2) / SD_preds2 
  std_res3 <- residuals(mod3) / SD_preds3
  
  ecdf_std_res1 <- ecdf(std_res1)
  ecdf_std_res2 <- ecdf(std_res2) 
  ecdf_std_res3 <- ecdf(std_res3)
  
  return(list(
    std_res1 = std_res1, std_res2 = std_res2, std_res3 = std_res3,
    ecdf_std_res1 = ecdf_std_res1, ecdf_std_res2 = ecdf_std_res2, ecdf_std_res3 = ecdf_std_res3
  ))
}

# returns smoothed ecdfs and params
get_parameters_smoothed_ecdfs_GAM <- function(mod1, mod2, mod3, SD_mod1, SD_mod2, SD_mod3, data_frame = data_baseline_dropna) {
  
  idx1 <- data_frame$class == 1
  idx2 <- data_frame$class == 2
  idx3 <- data_frame$class == 3
  
  SD_preds1 <- SD_mod1(data_frame[idx1, ])
  SD_preds2 <- SD_mod2(data_frame[idx2, ])
  SD_preds3 <- SD_mod3(data_frame[idx3, ])
  
  std_res1 <- residuals(mod1) / SD_preds1
  std_res2 <- residuals(mod2) / SD_preds2 
  std_res3 <- residuals(mod3) / SD_preds3
  
  ecdf_std_res1 <- ecdf(std_res1)
  ecdf_std_res2 <- ecdf(std_res2) 
  ecdf_std_res3 <- ecdf(std_res3)
  
  x1 <- sort(std_res1)
  y1 <- ecdf_std_res1(x1)
  x2 <- sort(std_res2)
  y2 <- ecdf_std_res2(x2)
  x3 <- sort(std_res3)
  y3 <- ecdf_std_res3(x3)
  
  
  scam_cdf1 <- scam(y1 ~ s(x1, bs = "mpi")) 
  scam_cdf2 <- scam(y2 ~ s(x2, bs = "mpi")) 
  scam_cdf3 <- scam(y3 ~ s(x3, bs = "mpi"))
  
  F1_fast <- fast_approx_scam_cdf(scam_cdf1, x1)
  F2_fast <- fast_approx_scam_cdf(scam_cdf2, x2)
  F3_fast <- fast_approx_scam_cdf(scam_cdf3, x3)
  
  
  return(list(
    std_res1 = std_res1, std_res2 = std_res2, std_res3 = std_res3,
    ecdf_std_res2 = ecdf_std_res2,
    scam_cdf2 = scam_cdf2,
    F2_fast = F2_fast, F1_fast = F1_fast, F3_fast = F3_fast
  ))
}

# fitted vs standardised residuals
create_residual_plots <- function(biomarker, GAM_models_data, params_data, data) {
  
  for (i in 1:3) {
    model_name <- paste0("mod", i)
    gam_model <- GAM_models_data[[model_name]]
    
    std_res <- params_data[[paste0("std_res", i)]]
    
    plot(fitted(gam_model), std_res,
         xlab = "Fitted values", 
         ylab = "Standardized residuals",
         main = paste(biomarker, model_name, "- Residuals vs Fitted"))
    
    abline(h = 0, col = "red", lty = 2)
    lines(lowess(fitted(gam_model), std_res), col = "blue")
    qqnorm(std_res, main = paste(biomarker, model_name, "- Q-Q Plot"))
    qqline(std_res, col = "red")
  }
  
  #par(mfrow = c(1, 1))  
}

GAM_models_data_HCI <- GAM_estimation("HCI")
GAM_models_data_TAU <- GAM_estimation("TAU_num")
GAM_models_data_PTAU <- GAM_estimation("PTAU_num")
GAM_models_data_ABETA <- GAM_estimation("ABETA_neg")

params_HCI <- get_parameters_ecdfs_GAM(GAM_models_data_HCI$mod1, GAM_models_data_HCI$mod2, GAM_models_data_HCI$mod3,
                                       GAM_models_data_HCI$SD_mod1, GAM_models_data_HCI$SD_mod2, GAM_models_data_HCI$SD_mod3)

params_PTAU <- get_parameters_ecdfs_GAM(GAM_models_data_PTAU$mod1, GAM_models_data_PTAU$mod2, GAM_models_data_PTAU$mod3,
                                        GAM_models_data_PTAU$SD_mod1, GAM_models_data_PTAU$SD_mod2, GAM_models_data_PTAU$SD_mod3)

params_TAU <- get_parameters_ecdfs_GAM(GAM_models_data_TAU$mod1, GAM_models_data_TAU$mod2, GAM_models_data_TAU$mod3,
                                       GAM_models_data_TAU$SD_mod1, GAM_models_data_TAU$SD_mod2, GAM_models_data_TAU$SD_mod3)

params_ABETA <- get_parameters_ecdfs_GAM(GAM_models_data_ABETA$mod1, GAM_models_data_ABETA$mod2, GAM_models_data_ABETA$mod3,
                                         GAM_models_data_ABETA$SD_mod1, GAM_models_data_ABETA$SD_mod2, GAM_models_data_ABETA$SD_mod3)

create_residual_plots("HCI", GAM_models_data_HCI, params_HCI, data_baseline_dropna)

create_residual_plots("PTAU", GAM_models_data_PTAU, params_PTAU, data_baseline_dropna)

create_residual_plots("TAU", GAM_models_data_TAU, params_TAU, data_baseline_dropna)

create_residual_plots("ABETA", GAM_models_data_ABETA, params_ABETA, data_baseline_dropna)





