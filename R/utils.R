#' @title Identify Out-of-Bounds Values (Outliers)
#' @description Uses the 1.5 * IQR rule to identify values that are potential outliers.
#' @param x A numeric vector.
#' @param na.rm Logical. Should missing values be removed? Defaults to TRUE.
#' @return A logical vector where TRUE indicates the value is out of bounds.
#' @export
out_of_bounds <- function(x, na.rm = TRUE) {
  q <- quantile(x, probs = c(0.25, 0.75), na.rm = na.rm)
  iqr <- IQR(x, na.rm = na.rm)
  lower_bound <- q[1] - 1.5 * iqr
  upper_bound <- q[2] + 1.5 * iqr
  x < lower_bound | x > upper_bound
}

# Create some sample data
set.seed(123)
sample_df <- data.frame(
  age = rnorm(100, 35, 5),
  extraversion = rnorm(100, 50, 10),
  neuroticism = rnorm(100, 45, 12),
  life_satisfaction = rnorm(100, 60, 15)
)

# Use your custom functions
descriptives(sample_df)

create_correlation_table(
  data = sample_df,
  vars = c("age", "extraversion", "neuroticism", "life_satisfaction")
)

# Build and display a hierarchical regression model
model1 <- lm(life_satisfaction ~ age, data = sample_df)
model2 <- lm(life_satisfaction ~ age + extraversion + neuroticism, data = sample_df)
lmtable(model1, model2)
