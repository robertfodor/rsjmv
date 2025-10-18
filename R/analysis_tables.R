#' Create an Advanced Correlation Table
#'
#' Conducts multivariate normality tests and creates a publication-ready
#' correlation table, automatically selecting Pearson or Spearman based on normality.
#'
#' @param data The input data.frame.
#' @param vars A character vector of variable names to include.
#' @param force_method Force a correlation method: "auto", "pearson", or "spearman".
#' @param stars Logical. If TRUE, use significance stars. If FALSE, use 'r, p=' format.
#' @param caption The title for the correlation table.
#' @param normality_p_threshold The significance threshold for normality tests.
#' @return Invisibly returns a list with detailed analysis results. Prints tables to the console.
#' @export
create_correlation_table <- function(data, vars, force_method = "auto", stars = TRUE, caption = NULL, normality_p_threshold = 0.05) {
  # This function depends on: dplyr, mvShapiroTest, psych, knitr, kableExtra, stringr, rlang
  if (!is.data.frame(data)) stop("A 'data' argumentumnak data.frame-nek kell lennie.")
  if (!all(vars %in% names(data))) {
    missing_vars <- vars[!vars %in% names(data)]
    stop(paste("A következő változók nem találhatóak meg az adatkeretben:", paste(missing_vars, collapse = ", ")))
  }
  
  data_subset <- data |>
    dplyr::select(dplyr::all_of(vars)) |>
    na.omit()
  
  if (nrow(data_subset) < 3) {
    stop("A hiányzó adatok eltávolítása után nem maradt elég teljes sor (minimum 3 szükséges).")
  }
  
  # Multivariate normality test
  mv_norm_result_df <- tryCatch({
    if (nrow(data_subset) <= 5000) {
      test <- mvShapiroTest::mvShapiro.Test(as.matrix(data_subset))
      data.frame(Statisztika = "Shapiro-Wilk W", Érték = round(test$statistic, 4), `p-érték` = test$p.value, check.names = FALSE)
    } else {
      data.frame(Statisztika = "Shapiro-Wilk W", Hiba = "Az esetszám a teszt elvégzéséhez túl magas (> 5000).", check.names = FALSE)
    }
  }, error = function(e) {
    data.frame(Statisztika = "Shapiro-Wilk W", Hiba = stringr::str_squish(e$message), check.names = FALSE)
  })
  
  print(knitr::kable(mv_norm_result_df, caption = "Többdimenziós normalitás (minden változó együtt)") |> kableExtra::kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = FALSE))
  
  # Pairwise normality tests
  pairwise_vars <- utils::combn(vars, 2, simplify = FALSE)
  pairwise_normality_list <- lapply(pairwise_vars, function(v) {
    test_data <- data_subset[, v]
    test_result <- suppressWarnings(mvShapiroTest::mvShapiro.Test(as.matrix(test_data)))
    data.frame(Var_Pair = paste(v, collapse = " & "), p_value = test_result$p.value)
  })
  pairwise_normality_df <- do.call(rbind, pairwise_normality_list)
  
  print(knitr::kable(pairwise_normality_df, caption = "Páronkénti többdimenziós normalitás teszt eredményei", col.names = c("Változópár", "p-érték")) |>
          kableExtra::kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = FALSE))
  
  recommended_method <- ifelse(any(pairwise_normality_df$p_value < normality_p_threshold, na.rm = TRUE), "spearman", "pearson")
  cor_method <- if (force_method == "auto") recommended_method else force_method
  
  if (force_method != "auto" && cor_method != recommended_method) {
    warning(paste0("Figyelem: A páronkénti normalitásvizsgálatok alapján a javasolt módszer a '", recommended_method, "', de a függvény a '", cor_method, "' módszer használatára lett kényszerítve."))
  }
  
  cor_results <- psych::corr.test(data_subset, method = cor_method, adjust = "none")
  cor_matrix <- cor_results$r
  cor_p_values <- cor_results$p
  
  n_vars <- length(vars)
  formatted_table <- matrix("", nrow = n_vars, ncol = n_vars, dimnames = list(vars, vars))
  for (i in 1:n_vars) {
    for (j in 1:n_vars) {
      if (i < j) {
        r_val <- cor_matrix[i, j]
        p_val <- cor_p_values[i, j]
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
  
  if (is.null(caption)) {
    method_note <- if (force_method != "auto") " (kényszerített módszer)" else ""
    caption <- paste0("Korrelációs mátrix (módszer: ", cor_method, ")", method_note)
  }
  table_footer <- if (stars) "Jelölés: * p < .05, ** p < .01, *** p < .001" else ""
  
  print(knitr::kable(formatted_table, caption = caption) |>
          kableExtra::kable_styling(bootstrap_options = c("striped", "hover", "condensed"), full_width = FALSE) |>
          kableExtra::footnote(general = table_footer, general_title = ""))
  
  invisible(list(multivariate_normality = mv_norm_result_df, pairwise_normality = pairwise_normality_df, recommended_method = recommended_method, used_method = cor_method, formatted_correlation_table = formatted_table, raw_correlation_results = cor_results))
}

#' Create a Formatted Hierarchical Regression Table
#'
#' Takes multiple linear models and presents them as a single
#' hierarchical regression table with model comparison statistics.
#'
#' @param ... One or more `lm` model objects.
#' @param caption The main title for the regression table.
#' @return A kable object representing the formatted table.
#' @export
lmtable <- function(..., caption = "Hierarchikus lineáris regresszió modell mutatószámai") {
  # This function depends on: purrr, broom, dplyr, effectsize, knitr, kableExtra
  models <- list(...)
  model_names <- paste("Modell", seq_along(models))
  names(models) <- model_names
  N <- tryCatch(stats::nobs(models[[1]]), error = function(e) "N/A")
  
  model_stats <- purrr::map_dfr(models, broom::glance, .id = "model_name")
  
  model_summary <- model_stats |>
    dplyr::mutate(
      r_squared_prev = dplyr::lag(.data$r.squared, default = 0),
      R2_change = .data$r.squared - r_squared_prev,
      df_change = .data$df - dplyr::lag(.data$df, default = 0),
      F_change = dplyr::case_when(
        df_change > 0 ~ (R2_change / df_change) / ((1 - .data$r.squared) / .data$df.residual),
        TRUE ~ NA_real_
      ),
      p_model_change = dplyr::case_when(
        df_change > 0 ~ stats::pf(F_change, df_change, .data$df.residual, lower.tail = FALSE),
        TRUE ~ NA_real_
      )
    ) |>
    dplyr::select(
      model_name, R2 = "r.squared", Adj_R2 = "adj.r.squared",
      "R2_change", F_omnibus = "statistic", p_omnibus = "p.value",
      "F_change", "p_model_change"
    )
  
  model_coeffs_raw <- purrr::map_dfr(models, ~ broom::tidy(.), .id = "model_name")
  std_betas_raw <- purrr::map_dfr(models, ~ effectsize::standardize_parameters(.), .id = "model_name")
  
  std_betas <- std_betas_raw |>
    dplyr::filter(.data$Parameter != "(Intercept)") |>
    dplyr::rename(term = "Parameter", beta = "Std_Coefficient") |>
    dplyr::select("model_name", "term", "beta")
  
  model_coeffs <- model_coeffs_raw |>
    dplyr::filter(.data$term != "(Intercept)") |>
    dplyr::left_join(std_betas, by = c("model_name", "term")) |>
    dplyr::select(model_name, term, b = "estimate", "beta", t = "statistic", p = "p.value")
  
  final_table_data <- model_coeffs |>
    dplyr::left_join(model_summary, by = "model_name") |>
    dplyr::group_by(.data$model_name) |>
    dplyr::mutate(dplyr::across(
      .cols = c("R2", "Adj_R2", "R2_change", "F_omnibus", "p_omnibus", "F_change", "p_model_change"),
      .fns = ~ ifelse(dplyr::row_number() == 1, .x, NA)
    )) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      `t, p =` = paste0(format(round(.data$t, 2), nsmall = 2), ", ", format.pval(.data$p, digits = 3, eps = 0.001, add_p = TRUE)),
      b = format(round(.data$b, 3), nsmall = 2),
      beta = format(round(.data$beta, 3), nsmall = 2),
      `R²` = ifelse(is.na(.data$R2), "", format(round(.data$R2, 3), nsmall = 3)),
      `Adj. R²` = ifelse(is.na(.data$Adj_R2), "", format(round(.data$Adj_R2, 3), nsmall = 3)),
      `Δ R²` = ifelse(is.na(.data$R2_change) | .data$R2_change <= 1e-9, "", format(round(.data$R2_change, 3), nsmall = 3)),
      `Modell F` = ifelse(is.na(.data$F_omnibus), "", paste0(format(round(.data$F_omnibus, 3), nsmall = 3), ", ", format.pval(.data$p_omnibus, digits = 3, eps = 0.001, add_p = TRUE))),
      `Δ F` = ifelse(is.na(.data$F_change), "", paste0(format(round(.data$F_change, 3), nsmall = 3), ", ", format.pval(.data$p_model_change, digits = 3, eps = 0.001, add_p = TRUE)))
    ) |>
    dplyr::select(
      `Modell` = "model_name", `Változó` = "term", "b", `β` = "beta",
      `t, p =`, `R²`, `Adj. R²`, `Δ R²`, `Modell F`, `Δ F`
    )
  
  full_caption <- paste0(caption, " (N = ", N, ")")
  
  knitr::kable(final_table_data, caption = full_caption, align = "llccccccccc")
}


#' Identify Out-of-Bounds Values (Outliers)
#'
#' Uses the 1.5 * IQR rule to identify values that are potential outliers.
#'
#' @param x A numeric vector.
#' @param na.rm Logical. Should missing values be removed? Defaults to TRUE.
#' @return A logical vector where TRUE indicates the value is out of bounds.
#' @export
out_of_bounds <- function(x, na.rm = TRUE) {
  q <- stats::quantile(x, probs = c(0.25, 0.75), na.rm = na.rm)
  iqr <- stats::IQR(x, na.rm = na.rm)
  lower_bound <- q[1] - 1.5 * iqr
  upper_bound <- q[2] + 1.5 * iqr
  x < lower_bound | x > upper_bound
}