#include <RcppArmadillo.h>
#include "GxEScan.h"

// [[Rcpp::depends(RcppArmadillo)]]

//' Function to perform a GEWIS scan
//' 
//' Function to perform a GEWIS scan
//' 
//' @param y
//' Vector of outcome values, 0 or 1, all others treated as missing
//' @param x
//' Matric of covariates. Last column is tested for interaction with gene
//' @param GeneticDataFilename
//' Name of file with genetic data in a binary file (measured or dosage data)
//' @param MapFilename
//' Name of map file associated with measured genetic data file
//' @param outFilename
//' Base filename for output files
//' @param dg
//' Perform D|G test
//' @param dgxe
//' Perform Test \deqn{\beta_{x} = 0}
//' @param twodf
//' Perform 2df test betaG = 0, betaGxE = 0
//' @param threedf
//' Perform 3df test betaG = 0, betaGxE = 0, betaD = 0
//' @param ge
//' Perform G|E test
//' @param caseOnly
//' Perform G|E test on cases
//' @param controlOnly
//' Perform G|E test on cases
//' @param dgge
//' Perform 2df test betaG = 0, betaD = 0
//' @return
//' A list containint 5 matrices
//' Allele frequencies in all subjects, cases, and controls
//' The number of subjects with complete data
//' Tne number of cases with complete data
//' The parameter estimates
//' The Z statistic for each estimate
//' The 2df chi-squared statistic for beta_g = 0 and beta_GxE = 0
//' @export
// [[Rcpp::export]]
Rcpp::List ScanSNPs(const arma::vec& y, const arma::mat& x, std::string GeneticDataFilename, std::string MapFilename, std::string outFilename,
                    bool dg, bool dgxe, bool twodf, bool threedf, bool ge, bool caseOnly, bool controlOnly, bool dgge) {
//  unsigned dg, unsigned dgxe = 1, unsigned twodf = 1, unsigned threedf = 1, unsigned ge = 1, unsigned caseOnly = 1, unsigned controlOnly = 1, unsigned dgge = 1) {
  CBedFile bedFile;
  CBinaryDosage dosageFile;
  CGeneticData *geneticData;
  CGxEScan gxeScan(y, x);
  Rcpp::List results;
  
  if (bedFile.ReadFile(GeneticDataFilename, y.n_elem, MapFilename) != 0) {
    if (dosageFile.ReadFile(GeneticDataFilename, y.n_elem, MapFilename) != 0) {
      Rcpp::Rcerr << "Unable to read genetic data files" << std::endl;
      return results;
    }
    geneticData = &dosageFile;
  } else {
    geneticData = &bedFile;
  }
  Rcpp::Rcout << "Read binary dosage file" << std::endl;
  gxeScan.SelectTests(dg, dgxe, twodf, threedf, ge, caseOnly, controlOnly, dgge);
  try {
    gxeScan.Scan(&bedFile, outFilename);
  } catch(...) {
    Rcpp::Rcerr << "Error" << std::endl;
    return NULL;
  }
  Rcpp::Rcout << "Scan complete" << std::endl;
  return Rcpp::List::create(
    Rcpp::Named("AlleleFreq") = gxeScan.Frequency(),
    Rcpp::Named("N") = gxeScan.N(),
    Rcpp::Named("NCases") = gxeScan.NCases(),
    Rcpp::Named("Beta") = gxeScan.Beta(),
//    Rcpp::Named("Order") = gxeScan.SortOrder(),
    Rcpp::Named("Z") = gxeScan.Z(),
    Rcpp::Named("ChiSq") = gxeScan.ChiSq());
}

//' Function to read a SNP from a measured genetic data file
//' 
//' Function to read a SNP from a measured genetic data file
//' 
//' @param BedFilename
//' Name of file with measured genetic data in plink format
//' @param numSubjects
//' Number of subjects with data in file
//' @param MapFilename
//' Name of map file associated with measured genetic data file
//' @param SNP
//' Name of SNP
//' @return
//' Chromosome Number
//' SNP Name
//' Location in base pairs
//' Dosage values
//' @export
// [[Rcpp::export]]
Rcpp::List ReadGeneticFile(std::string BedFilename, unsigned int numSubjects, std::string MapFilename, std::string SNP) {
  CBedFile bedFile;
//  arma::uvec snpLoc;
  unsigned int ui, uj;
  
  if (bedFile.ReadFile(BedFilename, numSubjects, MapFilename) != 0) {
    return Rcpp::List::create(Rcpp::Named("Error") = bedFile.ErrorString());
  }
//  snpLoc = arma::find(bedFile.MapFile().SNP() == "rs10000");
  for (ui = 0; ui < bedFile.NumSNPs(); ++ui) {
    if (bedFile.MapFile().SNP()[ui] == SNP)
      break;
  }

  bedFile.GetFirst();
  for (uj = 0; uj < ui; ++ uj)
    bedFile.GetNext();
  return Rcpp::List::create(Rcpp::Named("Chromosome") = bedFile.Chromosome(),
                            Rcpp::Named("SNP") = bedFile.SNPName(),
                            Rcpp::Named("BP") = bedFile.Location(),
                            Rcpp::Named("Dosage") = bedFile.Dosage()
  );
}
