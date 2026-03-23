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


# change gender and APOE4 to factor
data_baseline$PTGENDER <- factor(data_baseline$PTGENDER)
data_baseline$APOE4 <- factor(data_baseline$APOE4, levels = c(0,1,2),
                              labels = c("No e4", "One e4", "Two e4"))
summary(data_baseline)
dim(data_baseline)

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

# =============================================
# GAMs 
# =============================================


# function implementing the 2 step GAM estimation procedure
GAM_estimation <- function(biomarker, main_df = data_baseline_dropna){
  
  fmla <- as.formula(paste(biomarker, "~ s(AGE, bs = 'ps', k = 10) + PTGENDER + APOE4"))
  
  mod1 <- gam(fmla, data = subset(main_df, class == 1), method = "GCV.Cp")
  mod2 <- gam(fmla, data = subset(main_df, class == 2), method = "GCV.Cp")
  mod3 <- gam(fmla, data = subset(main_df, class == 3), method = "GCV.Cp")
  
  # Modelling variance
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




# functions obtaining standardised residuals, ecdfs or smoothed ecdfs
# =============================================


# returns a function approximating scam smoothed ecdfs and its inverse by interpolation
# to improve efficiency
fast_approx_scam_cdf <- function(model, x_data, n_grid = 400) {
  
  x_grid <- seq(min(x_data), max(x_data), length.out = n_grid)
  y_grid <- as.numeric(predict(model, newdata = data.frame(x1 = x_grid, x2 = x_grid, x3 = x_grid)))
  
  y_grid <- pmin(pmax(y_grid, 0), 1)
  
  F_fast <- approxfun(x_grid, y_grid, rule = 2)
  
  u <- !duplicated(y_grid)
  
  Finv_fast <- approxfun(y_grid[u], x_grid[u], rule = 2)
  
  return(list(
    F_fast = F_fast, 
    Finv_fast = Finv_fast
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
  
  # ecdfs 
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
  
 
  F2_fast <- fast_approx_scam_cdf(scam_cdf2, x2)$F_fast

  F1_inv_fast <- fast_approx_scam_cdf(scam_cdf1, x1)$Finv_fast
  F3_inv_fast <- fast_approx_scam_cdf(scam_cdf3, x3)$Finv_fast

  
  return(list(
    std_res1 = std_res1, std_res2 = std_res2, std_res3 = std_res3,
    ecdf_std_res2 = ecdf_std_res2,
    scam_cdf2 = scam_cdf2,
    F2_fast = F2_fast, F1_inv_fast = F1_inv_fast, F3_inv_fast = F3_inv_fast
  ))
}




# functions obtaining ROC triples 
# =============================================

# Optimized smoothed ECDF ROC triples function
cov_ROC_triples_GAM_smoothed_ecdf_fast_opt <- function(GAM_models_data_biomarker, params, cov_df, grid_len = 50) {

  mod1 <- GAM_models_data_biomarker$mod1
  mod2 <- GAM_models_data_biomarker$mod2
  mod3 <- GAM_models_data_biomarker$mod3

  SD_mod1 <- GAM_models_data_biomarker$SD_mod1
  SD_mod2 <- GAM_models_data_biomarker$SD_mod2
  SD_mod3 <- GAM_models_data_biomarker$SD_mod3
  
  F1_inv_fast <- params$F1_inv_fast
  F3_inv_fast <- params$F3_inv_fast
  F2 <- params$F2_fast

  pred1 <- as.numeric(predict(mod1, cov_df))
  pred2 <- as.numeric(predict(mod2, cov_df))
  pred3 <- as.numeric(predict(mod3, cov_df))
  
  SD_pred1 <- SD_mod1(cov_df)
  SD_pred2 <- SD_mod2(cov_df)
  SD_pred3 <- SD_mod3(cov_df)
  
  TPF1 <- TPF3 <- seq(0.001, 0.999, length.out = grid_len)
  
  c1_vec <- pred1 + SD_pred1 * F1_inv_fast(TPF1)
  c2_vec <- pred3 + SD_pred3 * F3_inv_fast(1 - TPF3)

  c1_std <- (c1_vec - pred2) / SD_pred2
  c2_std <- (c2_vec - pred2) / SD_pred2
  
  F2_c1 <- F2(c1_std)
  F2_c2 <- F2(c2_std)

  TPF2_mat <- matrix(0, nrow = grid_len, ncol = grid_len)
  Youden_vals <- matrix(0, nrow = grid_len, ncol = grid_len)
  
  for (i in 1:grid_len) {
    for (j in 1:grid_len) {
      if (c1_vec[i] < c2_vec[j]) {
        TPF2_mat[i, j] <- pmax(0, F2_c2[j] - F2_c1[i])
        Youden_vals[i, j] <- TPF1[i] + TPF2_mat[i, j] + TPF3[j] - 1
      }
      else {TPF2_mat[i, j] = 0 
      Youden_vals[i, j] <- -100}
    }
  }
  
  return(list(
    TPF1 = TPF1, TPF3 = TPF3, 
    TPF2_mat = TPF2_mat, 
    Youden_vals = Youden_vals, 
    c1_vec = c1_vec, c2_vec = c2_vec
  ))
}


# function obtaining VUS 
# =============================================

cov_VUS_trapz_GAM_smoothed_ecdf_fast_opt <- function(GAM_models_data_biomarker, params, cov_df, grid_len = 50) {
  
  roc_surface <- cov_ROC_triples_GAM_smoothed_ecdf_fast_opt(GAM_models_data_biomarker, params, cov_df, grid_len = grid_len)
  
  inner_integrals <- numeric(length(roc_surface$TPF1))
  for (i in 1:length(inner_integrals)) {
    inner_integrals[i] <- trapz(roc_surface$TPF3, roc_surface$TPF2_mat[i, ])
  }
  vus <- trapz(roc_surface$TPF1, inner_integrals)
  
  return(vus)
}

# VUS vs AGE for HCI 
# =============================================
combinations_GAM <- list(
  list(name = "Female, No APOE4", gender = "Female", apoe = "No e4"),
  list(name = "Female, One APOE4", gender = "Female", apoe = "One e4"),
  list(name = "Female, Two APOE4", gender = "Female", apoe = "Two e4"),
  list(name = "Male, No APOE4", gender = "Male", apoe = "No e4"),
  list(name = "Male, One APOE4", gender = "Male", apoe = "One e4"),
  list(name = "Male, Two APOE4", gender = "Male", apoe = "Two e4")
)


age_pred <- seq(60, 85, by = 1)


par(mfrow = c(2, 3))
for (comb in combinations_GAM) {
  vus_values <- numeric(length(age_pred))
  
  mod1 <- GAM_models_data_HCI$mod1
  mod2 <- GAM_models_data_HCI$mod2
  mod3 <- GAM_models_data_HCI$mod3
  
  SD_mod1 <- GAM_models_data_HCI$SD_mod1
  SD_mod2 <- GAM_models_data_HCI$SD_mod2
  SD_mod3 <- GAM_models_data_HCI$SD_mod3
  
  params <- get_parameters_ecdfs_GAM(mod1, mod2, mod3, SD_mod1, SD_mod2, SD_mod3)
  
  for (i in 1:length(age_pred)) {
    current_cov_df <- data.frame(
      AGE = age_pred[i],
      PTGENDER = factor(comb$gender, levels = levels(data_baseline$PTGENDER)),
      APOE4 = factor(comb$apoe, levels = levels(data_baseline$APOE4))
    )
    
    vus_values[i] <- cov_VUS_trapz_GAM_opt(GAM_models_data_HCI, params, current_cov_df)
  }
  
  plot(age_pred, vus_values, type = "l", lwd = 2,
       main = comb$name,
       xlab = "Age", ylab = "VUS",
       ylim = c(0, 1))
  grid()
}


# functions obtaining VUS 95% bootstrap CI 
# =============================================

bootstrap_VUS_GAM_smoothed <- function(biomarker, cov_df, age_range = age_pred, n_bootstrap = 1000) {
  bootstrap_results <- matrix(NA, nrow = n_bootstrap, ncol = length(age_range))
  
  idx1 <- which(data_baseline_dropna$class == 1)
  idx2 <- which(data_baseline_dropna$class == 2)
  idx3 <- which(data_baseline_dropna$class == 3)
  
  for (b in 1:n_bootstrap) {
    s1 <- sample(idx1, length(idx1), replace = TRUE)
    s2 <- sample(idx2, length(idx2), replace = TRUE)
    s3 <- sample(idx3, length(idx3), replace = TRUE)
    
    df_boot <- data_baseline_dropna
    df_boot[idx1,] <- data_baseline_dropna[s1,]
    df_boot[idx2,] <- data_baseline_dropna[s2,]
    df_boot[idx3,] <- data_baseline_dropna[s3,]
    
    GAM_models_data_boot <- try(GAM_estimation(biomarker, main_df = df_boot), silent = TRUE)
    if(inherits(GAM_models_data_boot, "try-error")) next
    
    test_pred <- try(predict(GAM_models_data_boot$mod1, cov_df), silent = TRUE)
    if(inherits(test_pred, "try-error")) next
    
    params <- try(get_parameters_smoothed_ecdfs_GAM(GAM_models_data_boot$mod1, GAM_models_data_boot$mod2, GAM_models_data_boot$mod3, 
                                                   GAM_models_data_boot$SD_mod1, GAM_models_data_boot$SD_mod2, GAM_models_data_boot$SD_mod3, data_frame = df_boot), silent = TRUE)
    if (inherits(params, "try-error")) next
    
    for (i in 1:length(age_pred)) {
      current_cov_df <- cov_df
      current_cov_df$AGE <- age_pred[i]
      bootstrap_results[b, i] <- cov_VUS_trapz_GAM_smoothed_ecdf_fast_opt(GAM_models_data_boot, params, current_cov_df)
    }
    cat("completed bootstrap", b)
  }
  
  lower_band <- numeric(length(age_range))
  upper_band <- numeric(length(age_range))
  
  for (i in 1:length(age_range)) {
    col_values <- bootstrap_results[, i]
    valid_values <- col_values[!is.na(col_values)]
    lower_band[i] <- quantile(valid_values, 0.025, na.rm = TRUE)
    upper_band[i] <- quantile(valid_values, 0.975, na.rm = TRUE)
  }
  
  return(list(lower = lower_band, upper = upper_band, bootstrap_results = bootstrap_results))
}




# VUS vs AGE bootstrap bands for all covariates
# =============================================

GAM_models_data_HCI <- GAM_estimation("HCI")
GAM_models_data_TAU <- GAM_estimation("TAU_num")
GAM_models_data_PTAU <- GAM_estimation("PTAU_num")
GAM_models_data_ABETA <- GAM_estimation("ABETA_neg")


par(mfrow = c(2, 3))
for (comb in combinations_GAM) {
  curr_cov_df <- data.frame(
    AGE = 0,
    PTGENDER = factor(comb$gender, levels = levels(data_baseline$PTGENDER)),
    APOE4 = factor(comb$apoe, levels = levels(data_baseline$APOE4))
  )
  
  bands <- bootstrap_VUS_GAM_smoothed("HCI", curr_cov_df, age_pred, n_bootstrap =1000)
  
  params <- get_parameters_smoothed_ecdfs_GAM(GAM_models_data_HCI$mod1, GAM_models_data_HCI$mod2, 
                                              GAM_models_data_HCI$mod3, GAM_models_data_HCI$SD_mod1, 
                                              GAM_models_data_HCI$SD_mod2, GAM_models_data_HCI$SD_mod3)
  
  vus_original <- numeric(length(age_pred))
  for (i in 1:length(age_pred)) {
    current_cov_df <- data.frame(
      AGE = age_pred[i],
      PTGENDER = factor(comb$gender, levels = levels(data_baseline$PTGENDER)),
      APOE4 = factor(comb$apoe, levels = levels(data_baseline$APOE4))
    )
    
    vus_original[i] <- cov_VUS_trapz_GAM_smoothed_ecdf_fast_opt(GAM_models_data_HCI, params, current_cov_df)
  }
  
  plot(age_pred, vus_original, type = "l", lwd = 2, col = "blue",
       main = comb$name, xlab = "Age", ylab = "VUS", ylim = c(0, 1))
  
  polygon(
    c(age_pred, rev(age_pred)),
    c(bands$lower, rev(bands$upper)),
    col = rgb(0, 0, 1, 0.2),
    border = NA
  )
  
  lines(age_pred, vus_original, lwd = 2, col = "blue")
  grid()
}





