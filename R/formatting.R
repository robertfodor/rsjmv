#' Custom APA Theme for ggplot2
#'
#' Applies the APA theme from the `jtools` package and removes panel and plot backgrounds for a cleaner look.
#'
#' @return A ggplot2 theme object with APA styling and blank backgrounds.
#' @export
#' @examples
#' library(ggplot2)
#' library(jtools)
#' ggplot(mtcars, aes(mpg, wt)) +
#'   geom_point() +
#'   theme_rfapa()
theme_rfapa <- function() {
    jtools::theme_apa() +
        ggplot2::theme(
            panel.background = ggplot2::element_blank(),
            plot.background  = ggplot2::element_blank()
        )
}