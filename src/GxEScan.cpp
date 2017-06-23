//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>
#include "GxEScan.h"
using namespace Rcpp;

//  Learn more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

CGxEScan::CGxEScan(const arma::vec& y, const arma::mat& x) {
  arma::uvec notMissing;
  arma::uvec ccorder;
  arma::vec ty;
  arma::mat tx;
  int ui;

  m_bDG = true;
  m_bGxE = true;
  m_b2df = true;
  m_b3df = true;
  m_bGE = true;
  m_bCaseOnly = true;
  m_bControlOnly = true;
  m_bDGGE = true;
  
  m_bModel1 = true;
  m_bModel2 = true;
  m_bModel3 = true;
  m_bModel4 = true;
  m_bModel5 = true;
  
  // Number of Covariates - not including intercept  
  m_numCov = x.n_cols;

  // The following identifies subjects with missing data
  // They will removed form the data

  // Copy Y
  ty = y;
  // Replace missing values of Y with 9
  ty(arma::find(ty > 1)).fill(9);
  ty(arma::find(ty < 0)).fill(9);
  ty(arma::find_nonfinite(ty)).fill(9);
  // Find rows in X with values and set corresponding Y to missing
  for (ui = 0; ui < m_numCov; ++ui)
    ty(arma::find_nonfinite(x.col(ui))).fill(9);

  // Find the case indices
  m_indicesCases = arma::find(ty == 1);
  m_numCases = m_indicesCases.n_elem;
  // Find the control indices
  m_indicesControls = arma::find(ty == 0);
  m_numControls = m_indicesControls.n_elem;

  // Save case data
  m_xCases.zeros(m_numCases, m_numCov + 3);
  m_xCases.col(0).ones();
  m_xCases.submat(0, 1, m_numCases - 1, m_numCov) = x.rows(m_indicesCases);
  
  // Save control data
  m_xControls.zeros(m_numControls, m_numCov + 3);
  m_xControls.col(0).ones();
  m_xControls.submat(0, 1, m_numControls - 1, m_numCov) = x.rows(m_indicesControls);
  
  m_geneCount.set_size(2, 3);
}
// Clear values and set to NA
void CGxEScan::InitializeValues() {
  m_betaDG = NA_REAL;
  m_zDG = NA_REAL;
  m_betaGxE = NA_REAL;
  m_zGxE = NA_REAL;
  m_chi2df = NA_REAL;
  m_chi3df = NA_REAL;
  m_betaGE = NA_REAL;
  m_zGE = NA_REAL;
  m_betaCaseOnly = NA_REAL;
  m_zCaseOnly = NA_REAL;
  m_betaControlOnly = NA_REAL;
  m_zControlOnly = NA_REAL;
  m_chiDGGE = NA_REAL;
  m_betaOut.fill(NA_REAL);
  m_zOut.fill(NA_REAL);
  m_chiSqOut.fill(NA_REAL);
}
// Allocate memory needed
int CGxEScan::AllocateMemory(unsigned int numSNPs) {
  unsigned int ui;
  unsigned int nb, nc;
  
  nb = 0;
  if (m_bDG)
    ++nb;
  if (m_bGxE)
    ++nb;
  if (m_bGE)
    ++nb;
  if (m_bCaseOnly)
    ++nb;
  if (m_bControlOnly)
    ++nb;
  nc = 0;
  if (m_b2df)
    ++nc;
  if (m_b3df)
    ++nc;
  if (m_bDGGE)
    ++nc;
  m_numCovSq = ((m_numCov + 3) * (m_numCov + 4)) / 2;
  try {
    m_xSqCases.set_size(m_numCases, m_numCovSq);
    m_xSqControls.set_size(m_numControls, m_numCovSq);
    // Logistic Regression
    m_betaInit.set_size(m_numCov + 3);
    m_beta.set_size(m_numCov + 3);
    m_betaDiff.set_size(m_numCov + 3);
    m_scoreConstant.set_size(m_numCov + 3);
    m_score.set_size(m_numCov + 3);
    m_abx1.set_size(m_numCases);
    m_abx2.set_size(m_numControls); 
    m_info.zeros(m_numCov + 3, m_numCov + 3);
    m_invInfo.zeros(m_numCov + 3, m_numCov + 3);
    // Polytomous Logistic Regression
    m_PLRbeta.set_size(m_numCov + 1, 2);
    m_PLRbetaDiff.set_size(2*(m_numCov + 1));
    m_PLRscoreConstant.set_size(2*(m_numCov + 1));
    m_PLRscore.set_size(2*(m_numCov + 1));
    m_PLRabx1.set_size(m_numCases, 2);
    m_PLRabx2.set_size(m_numControls, 2);
    m_PLRsumabx1.set_size(m_numCases);
    m_PLRsumabx2.set_size(m_numControls);
    m_PLRinfo.set_size(2*(m_numCov + 1), 2*(m_numCov + 1));
    m_PLRinvInfo.set_size(2*(m_numCov + 1), 2*(m_numCov + 1));
    // The results
    m_freq.zeros(numSNPs, 3);
    m_sort.set_size(numSNPs);
    m_betaOut.set_size(numSNPs, nb);
    m_zOut.set_size(numSNPs, nb);
    m_chiSqOut.set_size(numSNPs, nc);
    m_nOut.set_size(numSNPs);
    m_nCasesOut.set_size(numSNPs);
    m_defIndices.set_size(m_numCovSq);
  } catch (...) {
    Rcpp::Rcerr << "Error allocating memory" << std::endl;
    return 1;
  }
  for (ui = 0; ui < m_numCovSq; ++ui)
    m_defIndices(ui) = ui;
  return 0;
}
// Test if logistic model can converge
int CGxEScan::MinMaxTest() {
  arma::mat minCases;
  arma::mat maxCases;
  arma::mat minControls;
  arma::mat maxControls;
  unsigned int ui;
  
  // Find minimum and maximum of each covariate
  minCases = arma::min(m_xCases.cols(1, m_numCov));
  maxCases = arma::max(m_xCases.cols(1, m_numCov));
  minControls = arma::min(m_xControls.cols(1, m_numCov));
  maxControls = arma::max(m_xControls.cols(1, m_numCov));
  // Test if minimum in cases is less than control or vice versa
  for (ui = 0; ui < m_numCov; ++ui) {
    if ((minCases(ui) > maxControls(ui)) || (maxCases(ui) < minControls(ui))) {
      Rcpp::Rcerr << "At least one covariate can identify outcome" << std::endl;
      return 1;
    }
  }
  return 0;
}
// Standardize the covariates - help Newton-Raphson
int CGxEScan::Standardize() {
  int numRec = m_numCases + m_numControls;
  unsigned int ui;
  arma::uvec zeroStd;
  
  m_meanX = (arma::sum(m_xCases) + arma::sum(m_xControls)) / numRec;
  m_stdDevX = sqrt((arma::sum(square(m_xCases)) + arma::sum(m_xControls % m_xControls)) / numRec - square(m_meanX));
  zeroStd = find(m_stdDevX.cols(1, m_numCov) == 0);
  // Check that each covariate isn't a constant
  if (zeroStd.n_elem > 0) {
    Rcpp::Rcerr << "At least one covariate has standard deviation of 0" << std::endl;
    return 1;
  }
  // Subtract the mean and divide by the standard deviation
  for (ui = 0; ui < m_numCov;) {
    ++ui;
    m_xCases.col(ui) -= m_meanX(ui);
    m_xCases.col(ui) /= m_stdDevX(ui);
    m_xControls.col(ui) -= m_meanX(ui);
    m_xControls.col(ui) /= m_stdDevX(ui);
  }
  // Changed to make it easier to get the correct estimate of the intercept
  m_meanX(0) = -1;
  m_stdDevX(0) = 1;
  return 0;
}
// Save the individual XTX matrix - avoids recalulating for each SNP
void CGxEScan::XSquared() {
  unsigned int ui, uj, um, un;

  m_xgIndices.set_size(m_numCov + 2);
  for (ui = 0, um = 0; ui < m_numCov + 1; ++ui, um += m_numCov + 4) {
    um -= ui;
    for (uj = ui, un = um; uj < m_numCov + 1; ++uj, ++un) {
      m_xSqCases.col(un) = m_xCases.col(ui) % m_xCases.col(uj);
      m_xSqControls.col(un) = m_xControls.col(ui) % m_xControls.col(uj);
    }
    m_xgIndices(ui) = un;
  }
  m_xgIndices(ui) = m_xgIndices(ui - 1) + 2;
}
// Update the XTX data when genetic data is added
void CGxEScan::UpdateXSquared() {
  unsigned int ui;
  
  for (ui = 0; ui < m_numCov + 2; ++ui) {
    m_xSqCases.col(m_xgIndices(ui)) = m_xCases.col(ui) % m_xCases.col(m_numCov + 1);
    m_xSqCases.col(m_xgIndices(ui) + 1) = m_xCases.col(ui) % m_xCases.col(m_numCov + 2);
    //    m_xSqCases.col(m_xgIndices(ui)) = m_xCases.col(m_numCov + 1);
    m_xSqControls.col(m_xgIndices(ui)) = m_xControls.col(ui) % m_xControls.col(m_numCov + 1);
    m_xSqControls.col(m_xgIndices(ui) + 1) = m_xControls.col(ui) % m_xControls.col(m_numCov + 2);
  }
  --ui;
  m_xSqCases.col(m_xgIndices(ui) + 2) = m_xCases.col(m_numCov + 2) % m_xCases.col(m_numCov + 2);
  m_xSqControls.col(m_xgIndices(ui) + 2) = m_xControls.col(m_numCov + 2) % m_xControls.col(m_numCov + 2);
  
}

// Fit the model with the covariates only
int CGxEScan::CovariatesOnly() {
  arma::mat aL;
  
  m_notMissingCases = arma::find_finite(m_xCases.col(0));
  m_notMissingControls = arma::find_finite(m_xControls.col(0));
  m_beta.zeros();
  if (LogReg(m_numCov) != 0)
    return 1;
  m_betaInit = m_beta;
  m_covOnlyBeta = m_beta.subvec(0, m_numCov) / m_stdDevX.submat(0, 0, 0, m_numCov).t();
  m_covOnlyStat = m_beta.subvec(0, m_numCov) / arma::sqrt(m_invInfo.submat(0, 0, m_numCov, m_numCov).diag());
  aL = -(m_meanX.submat(0, 0, 0, m_numCov) / m_stdDevX.submat(0, 0, 0, m_numCov));
  m_covOnlyBeta(0) = arma::dot(aL.t(), m_betaInit.subvec(0, m_numCov));
  m_covOnlyStat(0) = arma::as_scalar(m_covOnlyBeta(0) * sqrt(arma::inv_sympd(aL * m_invInfo.submat(0, 0, m_numCov, m_numCov) * aL.t())));
  return 0;
}
// Add the dosage data to the case and control datasets
void CGxEScan::AddGeneticData(const arma::vec &dosage, unsigned int un) {
  double n1, n2;
  
  m_xCases.col(m_numCov + 1) = dosage.elem(m_indicesCases);
  m_xCases.col(m_numCov + 2) = dosage.elem(m_indicesCases) % m_xCases.col(m_numCov);
  m_notMissingCases = arma::find_finite(m_xCases.col(m_numCov + 1));
  m_xControls.col(m_numCov + 1) = dosage.elem(m_indicesControls);
  m_xControls.col(m_numCov + 2) = dosage.elem(m_indicesControls) % m_xControls.col(m_numCov);
  m_notMissingControls = arma::find_finite(m_xControls.col(m_numCov + 1));
  m_nCasesOut(un) = m_notMissingCases.size();
  m_nOut(un) = m_notMissingCases.size() + m_notMissingControls.size();
  n1 = arma::as_scalar(sum(m_xCases.submat(m_notMissingCases, m_defIndices.subvec(m_numCov + 1, m_numCov + 1))));
  n2 = arma::as_scalar(sum(m_xControls.submat(m_notMissingControls, m_defIndices.subvec(m_numCov + 1, m_numCov + 1))));
  m_freq(un, 0) = (n1 + n2) /(2*(m_notMissingCases.size() + m_notMissingControls.size()));
  m_freq(un, 1) = n1 /(2*m_notMissingCases.size());
  m_freq(un, 2) = n2 /(2*m_notMissingControls.size());
  UpdateXSquared();
}
// Assign probabilities
void CGxEScan::GetProbs(const CGeneticData *geneticData) {
  arma::uvec pIndices;
  
  pIndices.set_size(2);
  pIndices(0) = m_numCov + 1;
  pIndices(1) = m_numCov + 2;
  if (geneticData->Measured() == true) {
    m_xCases.col(m_numCov + 1).zeros();
    m_xCases.submat(find(geneticData->Dosage().elem(m_indicesCases) == 1), pIndices.subvec(0, 0)).ones();
    m_xCases.col(m_numCov + 2).zeros();
    m_xCases.submat(find(geneticData->Dosage().elem(m_indicesCases) == 2), pIndices.subvec(1, 1)).ones();
    m_xControls.col(m_numCov + 1).zeros();
    m_xControls.submat(find(geneticData->Dosage().elem(m_indicesControls) == 1), pIndices.subvec(0, 0)).ones();
    m_xControls.col(m_numCov + 2).zeros();
    m_xControls.submat(find(geneticData->Dosage().elem(m_indicesControls) == 2), pIndices.subvec(1, 1)).ones();
  }
}
// Polytomous Logistic Regression routines
void CGxEScan::CalculatePLRScoreConstants(unsigned int nc, bool bCases, bool bControls) {
  m_PLRscoreConstant.zeros();
  if (bCases == true) {
    m_PLRscoreConstant.subvec(0, nc) += m_xCases.submat(m_notMissingCases, m_defIndices.subvec(0,nc)).t() * m_xCases.submat(m_notMissingCases, m_defIndices.subvec(nc + 1, nc + 1));
    m_PLRscoreConstant.subvec(nc + 1, nc + nc + 1) += m_xCases.submat(m_notMissingCases, m_defIndices.subvec(0,nc)).t() * m_xCases.submat(m_notMissingCases, m_defIndices.subvec(nc + 2, nc + 2));
  }
  if (bControls == true) {
    m_PLRscoreConstant.subvec(0, nc) += m_xControls.submat(m_notMissingControls, m_defIndices.subvec(0,nc)).t() * m_xControls.submat(m_notMissingControls, m_defIndices.subvec(nc + 1, nc + 1));
    m_PLRscoreConstant.subvec(nc + 1, nc + nc + 1) += m_xControls.submat(m_notMissingControls, m_defIndices.subvec(0,nc)).t() * m_xControls.submat(m_notMissingControls, m_defIndices.subvec(nc + 2, nc + 2));
  }
}
// Polytomous Score
int CGxEScan::CalculatePLRScore(unsigned int nc, bool bCases, bool bControls) {
  if (bCases == true) {
    m_PLRabx1.rows(m_notMissingCases) = arma::exp(m_xCases.submat(m_notMissingCases, m_defIndices.subvec(0, nc))  * m_PLRbeta.rows(0, nc));
    m_PLRsumabx1.elem(m_notMissingCases) = m_PLRabx1.submat(m_notMissingCases, m_defIndices.subvec(0,0)) + m_PLRabx1.submat(m_notMissingCases, m_defIndices.subvec(1,1)) + 1;
    m_PLRabx1.submat(m_notMissingCases, m_defIndices.subvec(0, 0)) /= m_PLRsumabx1.elem(m_notMissingCases);
    m_PLRabx1.submat(m_notMissingCases, m_defIndices.subvec(1, 1)) /= m_PLRsumabx1.elem(m_notMissingCases);
  }
  if (bControls == true) {
    m_PLRabx2.rows(m_notMissingControls) = arma::exp(m_xControls.submat(m_notMissingControls, m_defIndices.subvec(0, nc))  * m_PLRbeta.rows(0, nc));
    m_PLRsumabx2.elem(m_notMissingControls) = m_PLRabx2.submat(m_notMissingControls, m_defIndices.subvec(0,0)) + m_PLRabx2.submat(m_notMissingControls, m_defIndices.subvec(1,1)) + 1;
    m_PLRabx2.submat(m_notMissingControls, m_defIndices.subvec(0, 0)) /= m_PLRsumabx2.elem(m_notMissingControls);
    m_PLRabx2.submat(m_notMissingControls, m_defIndices.subvec(1, 1)) /= m_PLRsumabx2.elem(m_notMissingControls);
  }
  m_PLRscore = m_PLRscoreConstant;
  if (bCases == true) {
    m_PLRscore.subvec(0, nc) -= m_xCases.submat(m_notMissingCases, m_defIndices.subvec(0,nc)).t() * m_PLRabx1.submat(m_notMissingCases, m_defIndices.subvec(0,0));
    m_PLRscore.subvec(nc + 1, nc + nc + 1) -= m_xCases.submat(m_notMissingCases, m_defIndices.subvec(0,nc)).t() * m_PLRabx1.submat(m_notMissingCases, m_defIndices.subvec(1,1));
  }
  if (bControls == true) {
    m_PLRscore.subvec(0, nc) -= m_xControls.submat(m_notMissingControls, m_defIndices.subvec(0,nc)).t() * m_PLRabx2.submat(m_notMissingControls, m_defIndices.subvec(0,0));
    m_PLRscore.subvec(nc + 1, nc + nc + 1) -= m_xControls.submat(m_notMissingControls, m_defIndices.subvec(0,nc)).t() * m_PLRabx2.submat(m_notMissingControls, m_defIndices.subvec(1,1));
  }
  return 0;
}
// Polytomous Information
int CGxEScan::CalculatePLRInformation(unsigned int nc, bool bCases, bool bControls) {
  unsigned int uj, uk, um, un; //ui
  
  if (bCases == true) {
    m_PLRsumabx1.elem(m_notMissingCases) = -m_PLRabx1.submat(m_notMissingCases, m_defIndices.subvec(0,0)) % m_PLRabx1.submat(m_notMissingCases, m_defIndices.subvec(1,1));
    m_PLRabx1.rows(m_notMissingCases) -= m_PLRabx1.rows(m_notMissingCases) % m_PLRabx1.rows(m_notMissingCases);
  }
  if (bControls == true) {
    m_PLRsumabx2.elem(m_notMissingControls) = -m_PLRabx2.submat(m_notMissingControls, m_defIndices.subvec(0,0)) % m_PLRabx2.submat(m_notMissingControls, m_defIndices.subvec(1,1));
    m_PLRabx2.rows(m_notMissingControls) -= m_PLRabx2.rows(m_notMissingControls) % m_PLRabx2.rows(m_notMissingControls);
  }

  m_PLRinfo.zeros();
  for (uj = 0, um = 0; uj < nc + 1; ++uj, um += m_numCov + 4) {
    um -= uj;
    for (uk = uj, un = um; uk < nc + 1; ++uk, ++un) {
      if (bCases == true) {
        m_PLRinfo(uj, uk) += arma::dot(m_PLRabx1.submat(m_notMissingCases, m_defIndices.subvec(0,0)), m_xSqCases.submat(m_notMissingCases, m_defIndices.subvec(un, un)));
        m_PLRinfo(uj + nc + 1, uk + nc + 1) += arma::dot(m_PLRabx1.submat(m_notMissingCases, m_defIndices.subvec(1,1)), m_xSqCases.submat(m_notMissingCases, m_defIndices.subvec(un, un)));
        m_PLRinfo(uj, uk + nc + 1) += arma::dot(m_PLRsumabx1.elem(m_notMissingCases), m_xSqCases.submat(m_notMissingCases, m_defIndices.subvec(un, un)));
      }
      if (bControls == true) {
        m_PLRinfo(uj, uk) += arma::dot(m_PLRabx2.submat(m_notMissingControls, m_defIndices.subvec(0,0)), m_xSqControls.submat(m_notMissingControls, m_defIndices.subvec(un, un)));
        m_PLRinfo(uj + nc + 1, uk + nc + 1) += arma::dot(m_PLRabx2.submat(m_notMissingControls, m_defIndices.subvec(1,1)), m_xSqControls.submat(m_notMissingControls, m_defIndices.subvec(un, un)));
        m_PLRinfo(uj, uk + nc + 1) += arma::dot(m_PLRsumabx2.elem(m_notMissingControls), m_xSqControls.submat(m_notMissingControls, m_defIndices.subvec(un, un)));
      }
      m_PLRinfo(uk, uj) = m_PLRinfo(uj, uk);
      m_PLRinfo(uk, uj + nc + 1) = m_PLRinfo(uj, uk + nc + 1);
      m_PLRinfo(uj + nc + 1, uk) = m_PLRinfo(uj, uk + nc + 1);
      m_PLRinfo(uk + nc + 1, uj) = m_PLRinfo(uj, uk + nc + 1);
      m_PLRinfo(uk + nc + 1, uj + nc + 1) = m_PLRinfo(uj + nc + 1, uk + nc + 1);
    }
  }
  return 0;
}
void CGxEScan::CollapsePLRInformation(unsigned int c1, unsigned int c2) {
  m_PLRinfo.cols(c1, c2) += 2 * m_PLRinfo.cols(c1 + m_numCov + 1, c2 + m_numCov + 1);
  m_PLRinfo.rows(c1, c2) += 2 * m_PLRinfo.rows(c1 + m_numCov + 1, c2 + m_numCov + 1);
  m_PLRinfo.cols(c1 + m_numCov + 1, c2 + m_numCov + 1).zeros();
  m_PLRinfo.rows(c1 + m_numCov + 1, c2 + m_numCov + 1).zeros();
  m_PLRinfo.submat(c1 + m_numCov + 1, c1 + m_numCov + 1, c2 + m_numCov + 1, c2 + m_numCov + 1).diag().ones();
}
void CGxEScan::CountGenes() {
  m_geneCount(0, 1) = arma::as_scalar(sum(m_xCases.submat(m_notMissingCases, m_defIndices.subvec(m_numCov + 1, m_numCov + 1))));
  m_geneCount(0, 2) = arma::as_scalar(sum(m_xCases.submat(m_notMissingCases, m_defIndices.subvec(m_numCov + 2, m_numCov + 2))));
  m_geneCount(0, 0) = m_notMissingCases.size() - m_geneCount(0, 1) - m_geneCount(0, 2);
  m_geneCount(1, 1) = arma::as_scalar(sum(m_xControls.submat(m_notMissingControls, m_defIndices.subvec(m_numCov + 1, m_numCov + 1))));
  m_geneCount(1, 2) = arma::as_scalar(sum(m_xControls.submat(m_notMissingControls, m_defIndices.subvec(m_numCov + 2, m_numCov + 2))));
  m_geneCount(1, 0) = m_notMissingControls.size() - m_geneCount(1, 1) - m_geneCount(1, 2);
}
int CGxEScan::PolyLogReg(bool bCases, bool bControls, bool bRestricted) {
  double n1, n2, d;
  unsigned int c1, c2 , c3;
  unsigned int ui;
  
  c2 = m_numCov;
  if (bRestricted == true) {
    c1 = 1;
    c3 = m_numCov + 1;
  } else {
    c1 = m_numCov;
    c3 = m_numCov + m_numCov;    
  }
  m_PLRbeta.zeros();
  n1 = 0;
  n2 = 0;
  d = 0;
  if (bCases == true) {
    n1 += m_geneCount(0, 1);
    n2 += m_geneCount(0, 2);
    d += m_geneCount(0, 0);
  }
  if (bControls == true) {
    n1 += m_geneCount(1, 1);
    n2 += m_geneCount(1, 2);
    d += m_geneCount(1, 0);
  }
  if (d < 5 || n1 < 5 || n2 < 5) {
    Rcpp::Rcerr << "too few 1\t" << d << '\t' << n1 << '\t' << n2 << std::endl;
    return 1;
  }
  m_PLRbeta(0,0) = log(n1 / d);
  m_PLRbeta(0,1) = log(n2 / d);
  CalculatePLRScoreConstants(m_numCov, bCases, bControls);
  for (ui = 0; ui < 25; ++ui) {
    CalculatePLRScore(m_numCov, bCases, bControls);
    m_PLRscore.subvec(c1, c2) += 2 * m_PLRscore.subvec(c1 + m_numCov + 1, c2 + m_numCov + 1);
    CalculatePLRInformation(m_numCov, bCases, bControls);
    CollapsePLRInformation(c1, c2);
    if (max(abs(m_PLRscore.subvec(0, c3))) < 1e-6)
      break;
    m_PLRbetaDiff.zeros();
    try {
      m_PLRbetaDiff.subvec(0, c3) = arma::solve(m_PLRinfo.submat(0, 0, c3, c3), m_PLRscore.subvec(0, c3), arma::solve_opts::no_approx);
    } catch (...) {
      return 1;
    }
    m_PLRbeta.submat(0, 0, m_numCov, 0) += m_PLRbetaDiff.subvec(0, m_numCov);
    m_PLRbeta.submat(0, 1, m_numCov, 1) += m_PLRbetaDiff.subvec(m_numCov + 1, m_numCov + m_numCov + 1);
    m_PLRbeta.submat(c1, 1, c2, 1) = 2 * m_PLRbeta.submat(c1, 0, c2, 0);
    //    m_PLRbeta(0,0) += m_PLRscore(0) / m_PLRinfo(0,0);
  }
  if (ui == 25)
    return 1;
  try {
    m_PLRinvInfo.submat(0, 0, c3, c3) = arma::inv_sympd(m_PLRinfo.submat(0, 0, c3, c3));
  } catch (...) {
    return 1;
  }
  //  CollapsePLRInformation(m_numCov, m_numCov);
  return 0;
}
// Logistic regression
int CGxEScan::LogReg(unsigned int nc, bool hw, bool bCases, bool bControls) {
  unsigned int ui, uj, uk, um, un;
//  arma::uvec nmc;
  
//  nmc = arma::find_finite(m_xSqCases.row(0));
  
  if (hw == true) {
    if (bCases == true) {
      m_scoreConstant.subvec(0, nc) = sum(m_xSqCases.submat(m_notMissingCases, m_xgIndices.subvec(0,nc))).t();
      if (bControls == true)
        m_scoreConstant.subvec(0, nc) += sum(m_xSqControls.submat(m_notMissingControls, m_xgIndices.subvec(0,nc))).t();
    } else {
      m_scoreConstant.subvec(0, nc) = sum(m_xSqControls.submat(m_notMissingControls, m_xgIndices.subvec(0,nc))).t();
    }
    m_scoreConstant *= 0.5;
  } else {
    m_scoreConstant.subvec(0, nc) = sum(m_xCases.submat(m_notMissingCases, m_defIndices.subvec(0, nc))).t();
  }
  for (ui = 0; ui < 25; ++ui) {
    // Calculating the score
    if (bCases == true) {
      m_abx1.elem(m_notMissingCases) = arma::exp(m_xCases.submat(m_notMissingCases, m_defIndices.subvec(0, nc)) * m_beta.subvec(0, nc));
      m_abx1.elem(m_notMissingCases) /= (m_abx1.elem(m_notMissingCases) + 1);
    }
    if (bControls == true) {
      m_abx2.elem(m_notMissingControls) = arma::exp(m_xControls.submat(m_notMissingControls, m_defIndices.subvec(0, nc)) * m_beta.subvec(0, nc));
      m_abx2.elem(m_notMissingControls) /= (m_abx2.elem(m_notMissingControls) + 1);
    }
    m_score.subvec(0, nc) = m_scoreConstant.subvec(0, nc);
    if(bCases == true)
      m_score.subvec(0, nc) -= m_xCases.submat(m_notMissingCases, m_defIndices.subvec(0, nc)).t() * m_abx1.elem(m_notMissingCases);
    if (bControls == true)
      m_score.subvec(0, nc) -= m_xControls.submat(m_notMissingControls, m_defIndices.subvec(0, nc)).t() * m_abx2.elem(m_notMissingControls);
    if (m_score.subvec(0,nc).is_finite() == false)
      return 1;
    // Calculating the information
    
    if(bCases == true)
      m_abx1.elem(m_notMissingCases) -= arma::square(m_abx1.elem(m_notMissingCases));
    if (bControls == true)
      m_abx2.elem(m_notMissingControls) -= arma::square(m_abx2.elem(m_notMissingControls));
    for (uj = 0, um = 0; uj < nc + 1; ++uj, um += m_numCov + 4) {
      um -= uj;
      for (uk = uj, un = um; uk < nc + 1; ++uk, ++un) {
        //      m_info(ui, uj) = sum(m_x2SqCases.col(0));
        if (bCases == true) {
          m_info(uj, uk) = arma::dot(m_abx1.elem(m_notMissingCases), m_xSqCases.submat(m_notMissingCases, m_defIndices.subvec(un, un)));
          if (bControls == true)
            m_info(uj, uk) += arma::dot(m_abx2.elem(m_notMissingControls), m_xSqControls.submat(m_notMissingControls, m_defIndices.subvec(un, un)));
        } else {
          m_info(uj, uk) = arma::dot(m_abx2.elem(m_notMissingControls), m_xSqControls.submat(m_notMissingControls, m_defIndices.subvec(un, un)));
        }
        m_info(uk, uj) = m_info(uj, uk);
      }
    }
    // Check for convergence. Done after information is calculated because most recent information
    // will be needed to calculate statistics
    if (arma::max(arma::abs(m_score.subvec(0, nc))) < 1e-8)
      break;
    // Newton-Raphson step
    try {
      m_betaDiff.subvec(0, nc) = arma::solve(m_info.submat(0, 0, nc, nc), m_score.subvec(0, nc), arma::solve_opts::no_approx);
    } catch (...) {
      return 1;
    }
    m_beta.subvec(0, nc) += m_betaDiff.subvec(0, nc);
  }
  // Did it not converge?
  if (ui == 25)
    return 1;
  try {
    m_invInfo.submat(0, 0, nc, nc) = arma::inv_sympd(m_info.submat(0, 0, nc, nc));
  } catch (...) {
    return 1;
  }
  return 0;
}
int CGxEScan::AlleleTests(unsigned int un) {
//  double b1, b2, b3;
//  double z1, z2, z3;
  
  if (m_bModel3 == true) {
    m_beta.zeros();
    if (LogReg(m_numCov, true) == 1) {
      m_betaGE = NA_REAL;
      m_zGE = NA_REAL;
      m_betaCaseOnly = NA_REAL;
      m_zCaseOnly = NA_REAL;
      m_betaControlOnly = NA_REAL;
      m_zControlOnly = NA_REAL;
      return 1;
    }
    m_betaGE = m_beta(m_numCov) / m_stdDevX(m_numCov);
    m_zGE = m_beta(m_numCov) / sqrt(0.5 * m_invInfo(m_numCov, m_numCov));
  }

  if (m_bModel4 == true) {
    m_beta.zeros();
    if (LogReg(m_numCov, true, true, false) == 0) {
      m_betaCaseOnly = m_beta(m_numCov) / m_stdDevX(m_numCov);
      m_zCaseOnly = m_beta(m_numCov) / sqrt(0.5 * m_invInfo(m_numCov, m_numCov));
    } else {
      m_betaCaseOnly = NA_REAL;
      m_zCaseOnly = NA_REAL;
    }
  }  

  if (m_bModel5 == true) {
    m_beta.zeros();
    if (LogReg(m_numCov, true, false, true) == 0) {
      m_betaControlOnly = m_beta(m_numCov) / m_stdDevX(m_numCov);
      m_zControlOnly = m_beta(m_numCov) / sqrt(0.5 * m_invInfo(m_numCov, m_numCov));
    } else {
      m_betaControlOnly = NA_REAL;
      m_zControlOnly = NA_REAL;
    }
  }  
  return 0;
}
int CGxEScan::PolyLogRegTests(unsigned int un) {
//  double b1, b2, b3;
//  double z1, z2, z3;
  
  if (m_bModel3 == true) {
    if (PolyLogReg() != 0)
      return 1;
    m_betaGE= m_PLRbeta(m_numCov) / m_stdDevX(m_numCov);
    m_zGE = m_PLRbeta(m_numCov) / sqrt(m_PLRinvInfo(m_numCov, m_numCov));
  }
  
  if (m_bModel4 == true) {
    if (PolyLogReg(true, false) != 0)
      return 1;
    m_betaCaseOnly = m_PLRbeta(m_numCov) / m_stdDevX(m_numCov);
    m_zCaseOnly = m_PLRbeta(m_numCov) / sqrt(m_PLRinvInfo(m_numCov, m_numCov));
  }
  
  if (m_bModel5 == true) {
    if (PolyLogReg(false, true) != 0)
      return 1;
    m_betaControlOnly = m_PLRbeta(m_numCov) / m_stdDevX(m_numCov);
    m_zControlOnly = m_PLRbeta(m_numCov) / sqrt(m_PLRinvInfo(m_numCov, m_numCov));
  }

  return 0;
}
int CGxEScan::RestrictedPolyLogRegTests(unsigned int un) {
//  double b1, b2, b3;
//  double z1, z2, z3;

  if (m_bModel3 == true) {
    if (PolyLogReg(true, true, true) != 0)
      return 1;
    m_betaGE = m_PLRbeta(m_numCov) / m_stdDevX(m_numCov);
    m_zGE = m_PLRbeta(m_numCov) / sqrt(m_PLRinvInfo(m_numCov, m_numCov));
  }
  
  if (m_bModel4 == true) {
    if (PolyLogReg(true, false, true) != 0)
      return 1;
    m_betaCaseOnly = m_PLRbeta(m_numCov) / m_stdDevX(m_numCov);
    m_zCaseOnly = m_PLRbeta(m_numCov) / sqrt(m_PLRinvInfo(m_numCov, m_numCov));
  }
  
  if (m_bModel5 == true) {
    if (PolyLogReg(false, true, true) != 0)
      return 1;
    m_betaControlOnly = m_PLRbeta(m_numCov) / m_stdDevX(m_numCov);
    m_zControlOnly = m_PLRbeta(m_numCov) / sqrt(m_PLRinvInfo(m_numCov, m_numCov));
  }

  return 0;
}
// D|G Model
int CGxEScan::DGTest(unsigned int un) {
  m_beta = m_betaInit;
  if (LogReg(m_numCov + 1) != 0) {
    m_beta = m_betaInit;
    return 1;
  }
  m_betaDG = m_beta(m_numCov + 1);
  m_zDG = m_beta(m_numCov + 1) / sqrt(m_invInfo(m_numCov + 1, m_numCov + 1));
  return 0;
}
// D|GxE Model
int CGxEScan::DGxETests(unsigned int un) {
  if (LogReg(m_numCov + 2) == 1)
    return 1;
  m_betaGxE = m_beta(m_numCov + 2) / m_stdDevX(m_numCov);
  m_zGxE = m_beta(m_numCov + 2) / sqrt(m_invInfo(m_numCov + 2, m_numCov + 2));
  m_chi2df = arma::as_scalar(m_beta.subvec(m_numCov + 1, m_numCov + 2).t() * arma::inv_sympd(m_invInfo.submat(m_numCov + 1, m_numCov + 1, m_numCov + 2, m_numCov + 2)) * m_beta.subvec(m_numCov + 1, m_numCov + 2));
  return 0;
}
// G|E Models
int CGxEScan::GETests(unsigned int un) {
  CountGenes();
  if (m_geneCount.min() > 4) {
    if (PolyLogRegTests(un) == 0)
      return 0;
    if (RestrictedPolyLogRegTests(un) == 0)
      return 0;
  }
  if (AlleleTests(un) != 0)
    return 1;
  return 0;
}
// Write results when beta is available
void CGxEScan::WriteBetaResults(std::string outFilename, unsigned int zNum, unsigned int nCol) {
  unsigned int ui;
  unsigned int n;
  unsigned int numSNPs;
  std::ofstream resultFile;
  arma::vec absZ;
  
  numSNPs = m_betaOut.n_rows;
  resultFile.open(outFilename);
  resultFile << "SNPID\tNMISS\tBETA\tSTAT\tP" << std::endl;
  absZ = abs(m_zOut.col(zNum));
  absZ(arma::find_nonfinite(absZ)).fill(-1);
  m_sort = arma::sort_index(absZ, "descend");
  for (ui = 0; ui < numSNPs; ++ui) {
    if (absZ(m_sort(ui)) == -1)
      break;
    switch(nCol) {
    case 0:
      n = m_nOut(m_sort(ui), 0);
      break;
    case 1:
      n = m_nCasesOut(m_sort(ui), 0);
      break;
    case 2:
      n = m_nOut(m_sort(ui), 0) - m_nCasesOut(m_sort(ui));
      break;
    default:
      break;
    }
    resultFile << m_sort(ui) + 1 << '\t' << n << '\t'
               << m_betaOut(m_sort(ui), zNum) << '\t'
               << m_zOut(m_sort(ui), zNum) << '\t'
               << 2 * R::pnorm5(std::abs(m_zOut(m_sort(ui), zNum)), 0, 1, 0, 0) << std::endl;
  }
  resultFile.close();
}
// Write results with chi-squared values
void CGxEScan::WriteChiResults(std::string outFilename, unsigned int df, unsigned int chiNum) {
  unsigned int ui;
  unsigned int numSNPs;
  std::ofstream resultFile;
  arma::vec chiSqSort;
  
  numSNPs = m_betaOut.n_rows;
  resultFile.open(outFilename);
  resultFile << "SNPID\tNMISS\tSTAT\tP" << std::endl;
  chiSqSort = m_chiSqOut.col(chiNum);
  chiSqSort(arma::find_nonfinite(chiSqSort)).fill(-1);
  m_sort = arma::sort_index(chiSqSort, "descend");
  for (ui = 0; ui < numSNPs; ++ui) {
    if (chiSqSort(m_sort(ui)) == -1)
      break;
    resultFile << m_sort(ui) + 1 << '\t'
               << m_nOut(m_sort(ui), 0) << '\t'
               << m_chiSqOut(m_sort(ui), chiNum) << '\t'
               << R::pchisq(m_chiSqOut(m_sort(ui), chiNum), df, 0, 0) << std::endl;
      
//    }
  }
  resultFile.close();
}
// Write files with test results
void CGxEScan::WriteTestResults(std::string outFilename) {
  unsigned int ui;
  
  ui = 0;
  if (m_bDG == true) {
    WriteBetaResults(outFilename + "_CC_DG.gxeout", ui);
    ++ui;
  }
  
  if (m_bGxE == true) {
    WriteBetaResults(outFilename + "_CC_GxE.gxeout", ui);
    ++ui;
  }
  
  if (m_bGE == true) {
    WriteBetaResults(outFilename + "_CC_GE.gxeout", ui);
    ++ui;
  }
  
  if (m_bCaseOnly == true) {
    WriteBetaResults(outFilename + "_Case_GE.gxeout", ui, 1);
    ++ui;
  }
  
  if (m_bControlOnly == true) {
    WriteBetaResults(outFilename + "_Cntl_GE.gxeout", ui, 2);
  }
  
  ui = 0;
  if (m_b2df == true) {
    WriteChiResults(outFilename + "_CC_2df.gxeout", 2, ui);
    ++ui;
  }
  if (m_b3df == true) {
    WriteChiResults(outFilename + "_CC_3df.gxeout", 3, ui);
    ++ui;
  }
  if (m_bDGGE == true) {
    WriteChiResults(outFilename + "_CC_DGGE.gxeout", 2, ui);
  }
}
// Select tests to perform
void CGxEScan::SelectTests(bool dg, bool gxe, bool twodf, bool threedf, bool ge, bool caseOnly, bool controlOnly, bool dgge) {
  m_bDG = dg;
  m_bGxE = gxe;
  m_b2df = twodf;
  m_b3df = threedf;
  m_bGE = ge;
  m_bCaseOnly = caseOnly;
  m_bControlOnly = controlOnly;
  m_bDGGE = dgge;
  m_bModel1 = m_bDG || m_bDGGE;
  m_bModel2 = m_bGxE || m_b2df || m_b3df;
  m_bModel3 = m_b3df || m_bGE;
  m_bModel4 = m_bCaseOnly || m_bDGGE;
  m_bModel5 = m_bControlOnly;
}
// Scan through the SNPs
void CGxEScan::Scan(CGeneticData *geneticData, std::string outFilename) {
  unsigned int ui, uj;
  arma::uvec missing;
  std::ofstream snpOutfile;

  if (AllocateMemory(geneticData->NumSNPs()) != 0)
    return;
  if (MinMaxTest() != 0)
    return;
  if (Standardize() != 0)
    return;
  XSquared();
  if (CovariatesOnly() != 0)
    return;
  
  if (outFilename.length() == 0)
    Rcpp::Rcerr << "No output file" << std::endl;
  else {
    snpOutfile.open(outFilename + ".snpinfo");
    if (!snpOutfile.good())
      return;
    snpOutfile << "SNPID\tCHR\tSNP\tBP\tA1" << std::endl;
  }
  
  InitializeValues();
  geneticData->GetFirst();
  for (ui = 0; ui < geneticData->NumSNPs(); ++ui) {
    AddGeneticData(geneticData->Dosage(), ui);
    if (m_freq(ui, 0) < 0.01 || m_freq(ui) > 0.99) {
      geneticData->GetNext();
      continue;
    }
    snpOutfile << ui + 1 << '\t'
               << geneticData->Chromosome() << '\t'
               << geneticData->SNPName() << '\t'
               << geneticData->Location() << '\t'
               << geneticData->FirstAllele() << std::endl;
    if (m_bModel1 == true)
      DGTest(ui);
    if (m_bModel2 == true && (m_bModel1 == false || !NumericVector::is_na(m_betaDG)))
      DGxETests(ui);
    if (m_bModel3 == true || m_bModel4 == true || m_bModel5 == true) {
      GetProbs(geneticData);
      GETests(ui);
    }
    uj = 0;
    if (m_bDG == true) {
      m_betaOut(ui, uj) = m_betaDG;
      m_zOut(ui, uj) = m_zDG;
      ++uj;
    }
    if (m_bGxE == true) {
      m_betaOut(ui, uj) = m_betaGxE;
      m_zOut(ui, uj) = m_zGxE;
      ++uj;
    }
    if (m_bGE == true) {
      m_betaOut(ui, uj) = m_betaGE;
      m_zOut(ui, uj) = m_zGE;
      ++uj;
    }
    if (m_bCaseOnly == true) {
      m_betaOut(ui, uj) = m_betaCaseOnly;
      m_zOut(ui, uj) = m_zCaseOnly;
      ++uj;
    }
    if (m_bControlOnly == true) {
      m_betaOut(ui, uj) = m_betaControlOnly;
      m_zOut(ui, uj) = m_zControlOnly;
    }
    uj = 0;
    if (m_b2df == true) {
      m_chiSqOut(ui, uj) = m_chi2df;
      ++uj;
    }
    if (m_b3df == true) {
      m_chiSqOut(ui, uj) = m_chi2df + m_zGE * m_zGE;
      ++uj;
    }
    if (m_bDGGE == true) {
      m_chiSqOut(ui, uj) = m_zDG * m_zDG + m_zGE * m_zGE;
    }
    geneticData->GetNext();
  }
  snpOutfile.close();
  WriteTestResults(outFilename);
}
