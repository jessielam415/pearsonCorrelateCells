#' Launch Shiny App for pearsonCorrelateCells
#'
#' A function that launches the Shiny app for pearsonCorrelateCells.
#' The purpose of this app is to demonstrates the capabilities of the
#' pearsonCorrelateCells package. Users can upload gene expression csv files
#' and visualize the correlation significance between the gene expression in the
#' uploaded csv file and each spot in the 10X Genomics Visium Sagittal Mouse
#' Brain dataset (Lab S, 2019)
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

