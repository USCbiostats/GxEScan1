#' Function to read in covariates
#' 
#' Function to read in covariate information from a file in plink format
#' @param filename
#' Name of covariate file
#' @param header
#' Inidicator if the file has a header line. The default value is TRUE.
#' @param used
#' List of covariates to keep. When no header is specified this is a list
#' of column numbers starring with from the third column. When a header
#' exists, the value is a list of column names. If not value of used is
#' entered, all covariates in the file will be used.
#' @param interaction
#' Column that contains covariate that will be used in iteraction term. This
#' value will be stored in the file column of the output dataframe. GxEScan
#' tests the covariate in the last column for gene environment interactions.
#' When the file has a header this value is the column name. Otherwise it is
#' the column number. If not interaction is supplied, the first variable in
#' used list is treated as the interacting covariate.
#' @param cmissing
#' Value that is be treated as a missing value. Only one missing value can be
#' supplied. If no value is entered NA is used.
#' @return
#' A dataframe with the first two columns being family and individual ID, FID and IID.
#' The remaining column are the value read in and kept as specified by the used
#' parameter. The last column will be the interacting covariate.
#' @export
ReadCovariates <- function(filename, header=TRUE, used, interaction, cmissing=NA) {

  Covariates <- read.table(filename,header=header,
                           na.strings = cmissing)

  if (header == FALSE) {
    colnames(Covariates) <- c("FID", "IID", paste("C",1:(ncol(Covariates)-2),sep = ""))
    if (!missing(interaction))
      interaction <- paste("C",interaction, sep ="")
    if (!missing(used))
      used <- paste("C",used, sep = "")
  }
  
  if (missing(used))
    used <- colnames(Covariates)[1:(ncol(Covariates)-2)]
  
  if (missing(interaction))
    interaction <- colnames(Covariates)[3];

  used <- c(used[ifelse(used == interaction, FALSE, TRUE)],interaction)

  used <- c("FID","IID",used)
  
  Covariates <- Covariates[,used]
  return(Covariates)
}
