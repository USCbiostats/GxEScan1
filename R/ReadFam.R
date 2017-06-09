#' Function to read in family data
#' 
#' Function to read in a plink formatted family file
#' @param filename
#' Name of file to read
#' @param continuous
#' Indicator if phenotype is continuous
#' @param zeroOne
#' Indicator if dichotomous phenotype is coded 0,1 instead of 1,2
#' @param pmissing
#' Value indiating a missing value
#' @return
#' Dataframe with family, individual ID, order records were read in, indicator if
#' subject is male, and phenotype.
#' @export
ReadFamily <- function(filename, continuous = FALSE, zeroOne = FALSE, pmissing = NA) {
  Family <- read.table(filename,
                       col.names = c("FID","IID","PID", "MID", "Sex","Phenotype1"),
                       colClasses = c("FID" = "character",
                                      "IID" = "character",
                                      "PID" = "character",
                                      "MID" = "character",
                                      "Sex" = "integer",
                                      "Phenotype1" = "double"
                                     ),
                       na.strings = pmissing
                      )
  Family$FileOrder = seq(1,nrow(Family))
  Family$Male <- with(Family, ifelse(Sex == 1, 1., ifelse(Sex == 2, 0., NA)))
  if (continuous == FALSE) {
    if (zeroOne == TRUE)
      Family$Phenotype <- with(Family, ifelse(Phenotype1 == 0, 1., ifelse(Phenotype1 == 1, 0., NA)))
    else
      Family$Phenotype <- with(Family, ifelse(Phenotype1 == 1, 0., ifelse(Phenotype1 == 2, 1., NA)))
  } else {
    Family$Phenotype = Family$Phenotype1
  }
  
  drops <- c("PID", "MID", "Sex","Phenotype1")
  return(Family[,!(names(Family) %in% drops)])
}
