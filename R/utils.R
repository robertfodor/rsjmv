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
```

### Step 3: Build and Install the Package

With the files in place, open your `psystats.Rproj` in RStudio. The following commands, run in the console, will build the documentation, check for issues, and install the package on your system.

1.  **Generate Documentation:** This command reads the `#'` comments (called 'Roxygen' comments) and creates the `NAMESPACE` file and the official documentation files.

```R
devtools::document()
```

2.  **Check the Package:** This is a crucial quality control step. It runs a series of automated checks to ensure your package is well-formed. It will flag any potential problems.

```R
devtools::check()
```

3.  **Install the Package:** This command bundles your package and installs it into your R library, just like any package from CRAN.

```R
devtools::install()
```

### Step 4: Use Your New Library

Once installed, you can load and use your package in any R script or session.

1.  Restart your R session (Ctrl+Shift+F10 in RStudio) for a clean start.
2.  Load your library and use your functions.

```R
# Load your custom library
library(psystats)

# Create some sample data
set.seed(123)
sample_df <- data.frame(
  age = rnorm(100, 35, 5),
  extraversion = rnorm(100, 50, 10),
  neuroticism = rnorm(100, 45, 12),
  life_satisfaction = rnorm(100, 60, 15)
)

# Use your custom functions
calculate_descriptives(sample_df)

create_correlation_table(
  data = sample_df,
  vars = c("age", "extraversion", "neuroticism", "life_satisfaction")
)

# Build and display a hierarchical regression model
model1 <- lm(life_satisfaction ~ age, data = sample_df)
model2 <- lm(life_satisfaction ~ age + extraversion + neuroticism, data = sample_df)
lmtable(model1, model2)
