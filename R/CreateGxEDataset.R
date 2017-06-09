#' Function to create a dataset for GxEScan
#' 
#' Function to take a family data frame from ReadFamily and a covariate data frame 
#' from ReadCovariate and create a list containing the matrices needed by the
#' GxEScan routine
#' 
#' @param family
#' Data frame is the format returned by ReadFamily
#' @param covariates
#' Data frame in the format returned by ReadCovariates
#' @param sex
#' Indicator to include sex variable from the family file as a covariate
#' @param sexInteraction
#' Indicator to use the sex variable from the family file as the interaction covariate.
#' If sexInteraction is set to TRUE, the value of parameter sex is ignored.
#' @param phenotype
#' Name of column in covariate data frame to use as the phenotype. If this value is
#' set to NA then the phenotype in the family file is used.
#' @export
CreateGxEScanDataset <- function(family, covariates, sex=FALSE, sexInteraction = FALSE, phenotype = NA) {
  if (sexInteraction == TRUE)
    sex = TRUE
  if (!missing(covariates)) {
    famcov <- merge(family, covariates, by=c("FID","IID"))
    famcov <- famcov[order(famcov$FileOrder),]
    famcov <- famcov[is.na(famcov$FileOrder) == FALSE,]
    drops <- c("FID", "IID", "FileOrder")
    famcov <- famcov[,!(names(famcov) %in% drops)]
  } else {
    if (sex == FALSE)
      return (NULL)
    if (is.na(phenotype) == FALSE)
      return (NULL)
    famcov <- fam
    y <- data.matrix(famcov$Phenotype)
    x <- data.matrix(famcov$Male)
    return (list(Y=y,X=x))
  }
# If sex isn't used, drop it from the data frame
  if (sex == FALSE)
    famcov <- famcov[,-which(names(famcov) %in% c("Male"))]
# Get the phenotype and remove it from the data frame
  if (is.na(phenotype) == TRUE) {
    y <- data.matrix(famcov[,"Phenotype"])
    famcov <- famcov[,-which(names(famcov) %in% c("Phenotype"))]
  } else {
    if (ncol(famcov) == 2)
      return (NULL)
    famcov <- famcov[,-which(names(famcov) %in% c("Phenotype"))]
    y <- data.matrix(famcov[,phenotype])
    famcov <- famcov[,-which(names(famcov) %in% c(phenotype))]
  }
# if sex is the interaction covariate move it to the end  
  if (sexInteraction == TRUE) {
    used <- c(colnames(famcov)[2:ncol(famcov)], colnames(famcov)[1])
    famcov <- famcov[,used]
  }
  
  x <- data.matrix(famcov);
  return (list(Y=y,X=x))
}