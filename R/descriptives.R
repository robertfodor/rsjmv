# This file contains the main descriptives function and its internal helpers.

# --- Internal CI Helper Functions ---
# (These are not exported and are only available within the package)

ci_bootstrap <- function(x, statistic_fun, ci_level, boot_reps) {
  x_complete <- na.omit(x)
  if (length(x_complete) < 2) return(c(NA, NA))
  stat_func_wrapper <- function(data, indices) statistic_fun(data[indices])
  boot_res <- boot::boot(data = x_complete, statistic = stat_func_wrapper, R = boot_reps)
  boot_ci <- tryCatch(
    boot::boot.ci(boot_res, conf = ci_level, type = "perc"),
    error = function(e) NULL
  )
  if (is.null(boot_ci)) return(c(NA, NA))
  return(c(boot_ci$percent[4], boot_ci$percent[5]))
}

ci_mean_normal <- function(x, ci_level) {
  n <- sum(!is.na(x))
  if (n < 2) return(c(NA, NA))
  mean_val <- mean(x, na.rm = TRUE)
  sem <- sd(x, na.rm = TRUE) / sqrt(n)
  z_crit <- qnorm(1 - (1 - ci_level) / 2)
  return(c(mean_val - z_crit * sem, mean_val + z_crit * sem))
}

ci_mean_t <- function(x, ci_level) {
  if (sum(!is.na(x)) < 2) return(c(NA, NA))
  t_test_res <- t.test(x, conf.level = ci_level)
  return(t_test_res$conf.int)
}

ci_sd_chisq <- function(x, ci_level) {
  x_complete <- na.omit(x)
  n <- length(x_complete)
  if (n < 2) return(c(NA, NA))
  s <- sd(x_complete)
  alpha <- 1 - ci_level
  chiSqLower <- stats::qchisq(p = alpha / 2, df = n - 1, lower.tail = TRUE)
  chiSqUpper <- stats::qchisq(p = 1 - alpha / 2, df = n - 1, lower.tail = TRUE)
  lowerBound <- sqrt((n - 1) * s^2 / chiSqUpper)
  upperBound <- sqrt((n - 1) * s^2 / chiSqLower)
  return(c(lowerBound, upperBound))
}

ci_var_chisq <- function(x, ci_level) {
  x_complete <- na.omit(x)
  n <- length(x_complete)
  if (n < 2) return(c(NA, NA))
  v <- var(x_complete)
  alpha <- 1 - ci_level
  chiSqLower <- stats::qchisq(p = alpha / 2, df = n - 1, lower.tail = TRUE)
  chiSqUpper <- stats::qchisq(p = 1 - alpha / 2, df = n - 1, lower.tail = TRUE)
  lowerBound <- (n - 1) * v / chiSqUpper
  upperBound <- (n - 1) * v / chiSqLower
  return(c(lowerBound, upperBound))
}


#' @title Calculate Detailed Descriptive Statistics
#' @description Generates a comprehensive descriptive statistics table from a data.frame.
#' @param data An input data.frame.
#' @param ci_level Confidence level for intervals (e.g., 0.95).
#' @param mean_ci_method Method for mean CI: "t", "normal", "bootstrap", "none".
#' @param sd_ci_method Method for SD CI: "chisq", "bootstrap", "none".
#' @param var_ci_method Method for variance CI: "chisq", "bootstrap", "none".
#' @param include_se Logical. Include standard errors.
#' @param include_quantiles Logical. Include quantiles (Min, Max, Q1, Q3, IQR, Median).
#' @param include_skewkurt Logical. Include skewness and kurtosis (requires 'e1071').
#' @param include_normality Logical. Include Shapiro-Wilk test.
#' @param boot_reps Number of bootstrap replications if used (requires 'boot').
#' @param output_kable Logical. If TRUE (default), prints a formatted kable table.
#' @param kable_caption Character. The caption for the kable table.
#' @param kable_digits Integer. Number of digits for the kable table.
#' @return If output_kable is TRUE, prints a kable object and returns the
#'   data.frame invisibly. If FALSE, returns the data.frame visibly.
#' @export
descriptives <- function(data,
                                   ci_level = 0.95,
                                   mean_ci_method = c("t", "normal", "bootstrap", "none"),
                                   sd_ci_method = c("chisq", "bootstrap", "none"),
                                   var_ci_method = c("chisq", "bootstrap", "none"),
                                   include_se = TRUE,
                                   include_quantiles = TRUE,
                                   include_skewkurt = TRUE,
                                   include_normality = TRUE,
                                   boot_reps = 1000,
                                   output_kable = TRUE,
                                   kable_caption = "Descriptive Statistics",
                                   kable_digits = 3) {
  mean_ci_method <- match.arg(mean_ci_method)
  sd_ci_method <- match.arg(sd_ci_method)
  var_ci_method <- match.arg(var_ci_method)

  process_variable <- function(x) {
    n <- sum(!is.na(x))
    n_missing <- sum(is.na(x))
    results <- list(N=n, Missing=n_missing, Mean=NA, SE_Mean=NA, CI_Mean_Lower=NA, CI_Mean_Upper=NA, SD=NA, SE_SD=NA, SD_CI_Lower=NA, SD_CI_Upper=NA, Variance=NA, SE_Variance=NA, Var_CI_Lower=NA, Var_CI_Upper=NA, Median=NA, Min=NA, Max=NA, IQR=NA, Q1=NA, Q3=NA, Skewness=NA, Skewness_SE=NA, Kurtosis=NA, Kurtosis_SE=NA, Shapiro_W=NA, Shapiro_p=NA)
    if (is.numeric(x) && n >= 2) {
      results$Mean <- mean(x, na.rm=TRUE)
      results$SD <- sd(x, na.rm=TRUE)
      results$Variance <- var(x, na.rm=TRUE)
      results$Median <- median(x, na.rm=TRUE)
      if (include_quantiles) {
        quants <- quantile(x, probs=c(0.25, 0.75), na.rm=TRUE, type=7)
        results$Min <- min(x, na.rm=TRUE)
        results$Max <- max(x, na.rm=TRUE)
        results$IQR <- IQR(x, na.rm=TRUE, type=7)
        results$Q1 <- quants[1]
        results$Q3 <- quants[2]
      }
      if (include_se) {
        results$SE_Mean <- results$SD / sqrt(n)
        results$SE_SD <- results$SD / sqrt(2 * (n - 1))
        results$SE_Variance <- results$Variance * sqrt(2 / (n - 1))
      }
      if (mean_ci_method != "none") {
        ci_m <- switch(mean_ci_method, "t"=ci_mean_t(x, ci_level), "normal"=ci_mean_normal(x, ci_level), "bootstrap"=ci_bootstrap(x, mean, ci_level, boot_reps), c(NA, NA))
        results$CI_Mean_Lower <- ci_m[1]
        results$CI_Mean_Upper <- ci_m[2]
      }
      if (sd_ci_method != "none") {
        ci_s <- switch(sd_ci_method, "chisq"=ci_sd_chisq(x, ci_level), "bootstrap"=ci_bootstrap(x, sd, ci_level, boot_reps), c(NA, NA))
        results$SD_CI_Lower <- ci_s[1]
        results$SD_CI_Upper <- ci_s[2]
      }
      if (var_ci_method != "none") {
        ci_v <- switch(var_ci_method, "chisq"=ci_var_chisq(x, ci_level), "bootstrap"=ci_bootstrap(x, var, ci_level, boot_reps), c(NA, NA))
        results$Var_CI_Lower <- ci_v[1]
        results$Var_CI_Upper <- ci_v[2]
      }
      if (include_skewkurt) {
        results$Skewness <- e1071::skewness(x, na.rm=TRUE, type=2)
        results$Skewness_SE <- sqrt(6*n*(n-1)/((n-2)*(n+1)*(n+3)))
        results$Kurtosis <- e1071::kurtosis(x, na.rm=TRUE, type=2)
        results$Kurtosis_SE <- sqrt(24*n*(n-1)^2/((n-3)*(n-2)*(n+3)*(n+5)))
      }
      if (include_normality && n > 2 && n <= 5000) {
        sw_test <- tryCatch(shapiro.test(x), error=function(e) NULL)
        if (!is.null(sw_test)) {
          results$Shapiro_W <- sw_test$statistic
          results$Shapiro_p <- sw_test$p.value
        }
      }
    }
    final_results <- results
    if (!include_se) final_results[c("SE_Mean", "SE_SD", "SE_Variance", "Skewness_SE", "Kurtosis_SE")] <- NULL
    if (mean_ci_method == "none") final_results[c("CI_Mean_Lower", "CI_Mean_Upper")] <- NULL
    if (sd_ci_method == "none") final_results[c("SD_CI_Lower", "SD_CI_Upper")] <- NULL
    if (var_ci_method == "none") final_results[c("Var_CI_Lower", "Var_CI_Upper")] <- NULL
    if (!include_quantiles) final_results[c("Min", "Max", "IQR", "Q1", "Q3", "Median")] <- NULL
    if (!include_skewkurt) final_results[c("Skewness", "Skewness_SE", "Kurtosis", "Kurtosis_SE")] <- NULL
    if (!include_normality) final_results[c("Shapiro_W", "Shapiro_p")] <- NULL
    return(unlist(final_results))
  }
  stats_list <- lapply(data, process_variable)
  stats_df <- as.data.frame(do.call(rbind, stats_list))
  stats_df <- cbind(Variable=names(data), stats_df)
  rownames(stats_df) <- NULL
  if (output_kable) {
    if (!requireNamespace("knitr", quietly=TRUE) || !requireNamespace("kableExtra", quietly=TRUE)) {
      warning("Packages 'knitr' and 'kableExtra' are required for kable output. Returning data.frame instead.", call.=FALSE)
      return(stats_df)
    }
    header_map <- c("Variable"="Variable", "N"="N", "Missing"="Missing", "Mean"="Mean", "SE_Mean"="SE (Mean)", "CI_Mean_Lower"="CI Lower (Mean)", "CI_Mean_Upper"="CI Upper (Mean)", "SD"="SD", "SE_SD"="SE (SD)", "SD_CI_Lower"="CI Lower (SD)", "SD_CI_Upper"="CI Upper (SD)", "Variance"="Variance", "SE_Variance"="SE (Var)", "Var_CI_Lower"="CI Lower (Var)", "Var_CI_Upper"="CI Upper (Var)", "Median"="Median", "Min"="Min", "Max"="Max", "IQR"="IQR", "Q1"="Q1 (25%)", "Q3"="Q3 (75%)", "Skewness"="Skewness", "Skewness_SE"="SE (Skew)", "Kurtosis"="Kurtosis", "Kurtosis_SE"="SE (Kurt)", "Shapiro_W"="Shapiro W", "Shapiro_p"="Shapiro p")
    current_headers <- names(stats_df)
    readable_headers <- unname(header_map[current_headers])
    kable_obj <- knitr::kable(stats_df, col.names=readable_headers, caption=kable_caption, digits=kable_digits, align=c("l", rep("r", length(current_headers) - 1))) |> kableExtra::kable_styling(bootstrap_options=c("striped", "hover", "condensed"), full_width=FALSE)
    print(kable_obj)
    return(invisible(stats_df))
  } else {
    return(stats_df)
  }
}
