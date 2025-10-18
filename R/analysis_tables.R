#' @title Create an Advanced Correlation Table
#' @description Conducts multivariate normality tests and creates a publication-ready
#'   correlation table, automatically selecting Pearson or Spearman based on normality.
#' @param data The input data.frame.
#' @param vars A character vector of variable names to include.
#' @param force_method Force a correlation method: "auto", "pearson", or "spearman".
#' @param stars Logical. If TRUE, use significance stars. If FALSE, use 'r, p=' format.
#' @param caption The title for the correlation table.
#' @param normality_p_threshold The significance threshold for normality tests.
#' @return Invisibly returns a list with detailed analysis results. Prints tables to the console.
#' @export
create_correlation_table <- function(data, vars, force_method = "auto", stars = TRUE, caption = NULL, normality_p_threshold = 0.05) {
  if (!is.data.frame(data)) stop("'data' must be a data.frame.")
  if (!all(vars %in% names(data))) {
    missing_vars <- vars[!vars %in% names(data)]
    stop(paste("The following variables were not found:", paste(missing_vars, collapse = ", ")))
  }
  data_subset <- data[, vars, drop = FALSE] |> na.omit()
  if (nrow(data_subset) < 3) stop("Not enough complete cases (minimum 3 required).")

  # Multivariate normality test for all variables
  mv_norm_result_df <- tryCatch({
    if (nrow(data_subset) <= 5000) {
      test <- mvShapiroTest::mvShapiro.Test(as.matrix(data_subset))
      data.frame(Statistic = "Shapiro-Wilk W", Value = round(test$statistic, 4), `p-value` = test$p.value, check.names = FALSE)
    } else {
      data.frame(Statistic = "Shapiro-Wilk W", Error = "Sample size > 5000, test not performed.", check.names = FALSE)
    }
  }, error = function(e) data.frame(Statistic = "Shapiro-Wilk W", Error = stringr::str_squish(e$message), check.names = FALSE))
  print(knitr::kable(mv_norm_result_df, caption = "Multivariate Normality (All Variables)") |> kableExtra::kable_styling())

  # Pairwise normality tests
  pairwise_vars <- combn(vars, 2, simplify = FALSE)
  pairwise_normality_list <- lapply(pairwise_vars, function(v) {
    test_result <- suppressWarnings(mvShapiroTest::mvShapiro.Test(as.matrix(data_subset[, v])))
    data.frame(Variable_Pair = paste(v, collapse = " & "), p_value = test_result$p.value)
  })
  pairwise_normality_df <- do.call(rbind, pairwise_normality_list)
  print(knitr::kable(pairwise_normality_df, caption = "Pairwise Multivariate Normality", col.names = c("Variable Pair", "p-value")) |> kableExtra::kable_styling())

  # Select correlation method
  recommended_method <- ifelse(any(pairwise_normality_df$p_value < normality_p_threshold, na.rm = TRUE), "spearman", "pearson")
  cor_method <- if (force_method == "auto") recommended_method else force_method
  if (force_method != "auto" && cor_method != recommended_method) {
    warning(paste("Normality tests suggest '", recommended_method, "', but '", cor_method, "' was forced.", sep=""))
  }

  # Calculate correlations
  cor_results <- psych::corr.test(data_subset, method = cor_method, adjust = "none")
  formatted_table <- .format_correlation_matrix(cor_results$r, cor_results$p, vars, stars)

  # Render table
  if (is.null(caption)) caption <- paste("Correlation Matrix (Method: ", cor_method, ")", sep="")
  table_footer <- if (stars) "Note: * p < .05, ** p < .01, *** p < .001" else ""
  print(knitr::kable(formatted_table, caption = caption) |> kableExtra::kable_styling() |> kableExtra::footnote(general = table_footer, general_title = ""))

  invisible(list(multivariate_normality = mv_norm_result_df, pairwise_normality = pairwise_normality_df, used_method = cor_method, correlation_table = formatted_table, raw_results = cor_results))
}

# Internal helper to format the correlation matrix
.format_correlation_matrix <- function(r_matrix, p_matrix, vars, stars) {
  n_vars <- length(vars)
  formatted_table <- matrix("", nrow = n_vars, ncol = n_vars, dimnames = list(vars, vars))
  for (i in 1:n_vars) {
    for (j in 1:n_vars) {
      if (i < j) {
        r_val <- r_matrix[i, j]
        p_val <- p_matrix[i, j]
        if (stars) {
          sig_stars <- dplyr::case_when(p_val < 0.001 ~ "***", p_val < 0.01 ~ "**", p_val < 0.05 ~ "*", TRUE ~ "")
          formatted_table[j, i] <- paste0(format(round(r_val, 3), nsmall = 3), sig_stars)
        } else {
          formatted_table[j, i] <- paste0(format(round(r_val, 3), nsmall = 3), ", ", format.pval(p_val, digits = 3, eps = 0.001, add_p = TRUE))
        }
      } else if (i == j) {
        formatted_table[i, j] <- "-"
      }
    }
  }
  return(formatted_table)
}


#' @title Create a Formatted Hierarchical Regression Table
#' @description Takes multiple linear models and presents them as a single
#'   hierarchical regression table with model comparison statistics.
#' @param ... One or more `lm` model objects.
#' @param caption The main title for the regression table.
#' @return A kable object representing the formatted table.
#' @export
lmtable <- function(..., caption = "Hierarchical Linear Regression Results") {
  models <- list(...)
  model_names <- paste("Model", seq_along(models))
  names(models) <- model_names
  N <- tryCatch(nobs(models[[1]]), error = function(e) "N/A")

  model_stats <- purrr::map_dfr(models, broom::glance, .id = "model_name")
  model_summary <- model_stats |>
    dplyr::mutate(
      R2_change = r.squared - dplyr::lag(r.squared, default = 0),
      df_change = df - dplyr::lag(df, default = 0),
      F_change = dplyr::if_else(df_change > 0, (R2_change / df_change) / ((1 - r.squared) / df.residual), NA_real_),
      p_model_change = dplyr::if_else(df_change > 0, pf(F_change, df_change, df.residual, lower.tail = FALSE), NA_real_)
    ) |>
    dplyr::select(model_name, R2=r.squared, Adj_R2=adj.r.squared, R2_change, F_omnibus=statistic, p_omnibus=p.value, F_change, p_model_change)

  model_coeffs <- purrr::map_dfr(models, ~ broom::tidy(.), .id = "model_name") |>
    dplyr::filter(term != "(Intercept)") |>
    dplyr::left_join(
      purrr::map_dfr(models, ~ effectsize::standardize_parameters(.), .id = "model_name") |> dplyr::select(model_name, term=Parameter, beta=Std_Coefficient),
      by = c("model_name", "term")
    ) |>
    dplyr::select(model_name, term, b=estimate, beta, t=statistic, p=p.value)

  final_table_data <- model_coeffs |>
    dplyr::left_join(model_summary, by = "model_name") |>
    dplyr::group_by(model_name) |>
    dplyr::mutate(dplyr::across(c(R2, Adj_R2, R2_change, F_omnibus, p_omnibus, F_change, p_model_change), ~ ifelse(dplyr::row_number() == 1, .x, NA))) |>
    dplyr::ungroup() |>
    dplyr::transmute(
      Model = model_name,
      Variable = term,
      b = format(round(b, 2), nsmall = 2),
      `β` = format(round(beta, 2), nsmall = 2),
      `t, p` = paste0(format(round(t, 2), nsmall = 2), ", ", format.pval(p, digits = 3, eps = 0.001)),
      `R²` = ifelse(is.na(R2), "", format(round(R2, 3), nsmall = 3)),
      `Adj. R²` = ifelse(is.na(Adj_R2), "", format(round(Adj_R2, 3), nsmall = 3)),
      `Δ R²` = ifelse(is.na(R2_change) | R2_change == 0, "", format(round(R2_change, 3), nsmall = 3)),
      `Model F` = ifelse(is.na(F_omnibus), "", paste0(format(round(F_omnibus, 2), nsmall = 2), ", ", format.pval(p_omnibus, digits = 3, eps = 0.001))),
      `Δ F` = ifelse(is.na(F_change), "", paste0(format(round(F_change, 2), nsmall = 2), ", ", format.pval(p_model_change, digits = 3, eps = 0.001)))
    )

  full_caption <- paste0(caption, " (N = ", N, ")")
  knitr::kable(final_table_data, caption = full_caption, align = "llccccccccc") |>
    kableExtra::kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = FALSE)
}
