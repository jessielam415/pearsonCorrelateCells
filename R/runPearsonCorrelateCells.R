#' Launch Shiny App for pearsonCorrelateCells
#'
#' A function that launches the Shiny app for pearsonCorrelateCells.
#' The purpose of this app will be added
#'
#' @return No return value but open up a Shiny page.
#'
#' @examples
#' \dontrun{
#'
#' pearsonCorrelateCells::runPearsonCorrelateCells()
#' }
#'
#' @references
#' Grolemund, G. (2015). Learn Shiny - Video Tutorials. \href{https://shiny.rstudio.com/tutorial/}{Link}
#'
#' @export
#' @importFrom shiny runApp

runPearsonCorrelateCells <- function() {
  appDir <- system.file("shiny-scripts",
                        package = "pearsonCorrelateCells")
  actionShiny <- shiny::runApp(appDir, display.mode = "normal")
  return(actionShiny)
}
# [END]

