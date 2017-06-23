#ifndef GXELOGREG_H
#define GXELOGREG_H 1
#ifndef BEDFILE_H
#include "BedFile.h"
#endif
#ifndef BINARYDOSAGE_H
#include "BinaryDosage.h"
#endif
#include <RcppArmadillo.h>

class CGxEScan {
protected:
  // The data
  unsigned int m_numCases;
  unsigned int m_numControls;
  unsigned int m_numCov;
  unsigned int m_numCovSq;
  unsigned int m_numSNPs;
  
  bool m_bDG;
  bool m_bGxE;
  bool m_b2df;
  bool m_b3df;
  bool m_bGE;
  bool m_bCaseOnly;
  bool m_bControlOnly;
  bool m_bDGGE;
  
  bool m_bModel1;
  bool m_bModel2;
  bool m_bModel3;
  bool m_bModel4;
  bool m_bModel5;
  
  double m_betaDG;
  double m_betaGxE;
  double m_betaGE;
  double m_betaCaseOnly;
  double m_betaControlOnly;

  double m_zDG;
  double m_zGxE;
  double m_chi2df;
  double m_chi3df;
  double m_zGE;
  double m_zCaseOnly;
  double m_zControlOnly;
  double m_chiDGGE;
  
  arma::uvec m_indicesCases;
  arma::uvec m_indicesControls;
  arma::uvec m_defIndices;
  arma::uvec m_xgIndices;

  arma::mat m_xCases;  // Covariates for cases with no missing data
  arma::mat m_xControls;  // Covariates for cases with no missing data
  arma::mat m_xSqCases;
  arma::mat m_xSqControls;
  arma::uvec m_notMissingCases;
  arma::uvec m_notMissingControls;
  // Logistic Regression
  arma::vec m_betaInit;
  arma::vec m_beta;
  arma::vec m_betaDiff;
  arma::vec m_scoreConstant;
  arma::vec m_score;
  arma::vec m_abx1;
  arma::vec m_abx2;
  arma::mat m_info;
  arma::mat m_invInfo;
  // Polytomous Logistic Regression  
  arma::mat m_PLRbeta;
  arma::vec m_PLRbetaDiff;
  arma::vec m_PLRscoreConstant;
  arma::vec m_PLRscore;
  arma::mat m_PLRabx1;
  arma::mat m_PLRabx2;
  arma::vec m_PLRsumabx1;
  arma::vec m_PLRsumabx2;
  arma::mat m_PLRinfo;
  arma::mat m_PLRinvInfo;
  // Allele information
  arma::mat m_geneCount;
  // Standardizing values
  arma::mat m_meanX;
  arma::mat m_stdDevX;
  // Allele frequencies
  arma::mat m_freq;
  // The results
  arma::vec m_covOnlyBeta;
  arma::vec m_covOnlyStat;
  arma::uvec m_sort;
  arma::mat m_betaOut;
  arma::mat m_zOut;
  arma::mat m_chiSqOut;
  arma::vec m_nOut;
  arma::vec m_nCasesOut;

  void InitializeValues();  
  int AllocateMemory(unsigned int numSNPs);
  int MinMaxTest();
  int Standardize();
  void XSquared();
  void UpdateXSquared();
  void AddGeneticData(const arma::vec &dosage, unsigned int un);
  int CovariatesOnly();
  void GetProbs(const CGeneticData *geneticData);
  void CalculatePLRScoreConstants(unsigned int nc, bool bCases = true, bool bControls = true);
  int CalculatePLRScore(unsigned int nc, bool bCases = true, bool bControls = true);
  int CalculatePLRInformation(unsigned int nc, bool bCases = true, bool bControls = true);
  void CollapsePLRInformation(unsigned int c1, unsigned int c2);
  int PolyLogReg(bool bCases = true, bool bControls = true, bool bRetricted = false);
  int LogReg(unsigned int nc, bool hw = false, bool bCases = true, bool bControls = true);
  void CountGenes();
  int AlleleTests(unsigned int un);
  int PolyLogRegTests(unsigned int un);
  int RestrictedPolyLogRegTests(unsigned int un);
  int DGTest(unsigned int un);
  int DGxETests(unsigned int un);
  int GETests(unsigned int un);
  void WriteBetaResults(std::string outFileName, unsigned int zNum, unsigned int nCol = 0);
  void WriteChiResults(std::string outFileName, unsigned int df, unsigned int chiNum);
  void WriteTestResults(std::string outFileName);
public:
  CGxEScan(const arma::vec& y, const arma::mat& x);
  virtual ~CGxEScan() {}
  
  void SelectTests(bool dg, bool gxe, bool twodf, bool threedf, bool ge, bool caseOnly, bool controlOnly, bool dgge);
  void Scan(CBinaryGeneticData *geneticData, std::string outFilename);
  void Summary();
  const arma::mat &Cases() const { return m_xCases; }
  const arma::mat &Controls() const { return m_xControls; }
  const arma::mat &CasesSq() const { return m_xSqCases; }
  const arma::mat &ControlsSq() const { return m_xSqControls; }
  const arma::mat &Mean() const { return m_meanX; }
  const arma::mat &StdDev() const { return m_stdDevX; }
  const arma::uvec &SortOrder() const { return m_sort; }
  const arma::mat &Beta() const { return m_betaOut; }
  const arma::mat &Z() const { return m_zOut; }
  const arma::mat &ChiSq() const { return m_chiSqOut; }
  const arma::vec &N() const { return m_nOut; }
  const arma::vec &NCases() const { return m_nCasesOut; }
  const arma::mat &Frequency() const { return m_freq; }
};
#endif