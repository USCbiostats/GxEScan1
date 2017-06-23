#include <RcppArmadillo.h>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cstring>
#include "BinaryDosage.h"

using namespace Rcpp;

// ********************************************************
//                CBinaryDosage
// Class to read in data form binary dosage files
// ********************************************************
// Constructor
CBinaryDosage::CBinaryDosage() : CBinaryGeneticData() {
	m_bMeasured = false;
	m_headerSize = 8;
	m_numSNPsUsed = 0;
	m_versionNumber = 0;
}
// Destructor
CBinaryDosage::~CBinaryDosage() {
}

// Read the file head and determine if file is valid
int CBinaryDosage::ReadFile(const std::string &_filename) {
//	unsigned long long actualSize;
//	unsigned long long expectedSize;
	const char bose[4] = { 'b', 'o', 's', 'e' };
	const char version_1_1[4] = { 0x0, 0x1, 0x0, 0x1 };
	const char version_1_2[4] = { 0x0, 0x1, 0x0, 0x2 };
	const char version_3_1[4] = { 0x0, 0x3, 0x0, 0x1 };
	char header[8];
	unsigned int expNumSub;

	if (OpenFile(_filename, true))
		return 1;

	m_infile.read(header, 4);
	if (memcmp(header, bose, 4)) {
		Rcpp::Rcerr << "File is not a binary dosage file" << std::endl;
		return 1;
	}
	m_infile.read(header, 4);
	if (memcmp(header, version_1_1, 4) == 0) {
		m_bProbabilities = false;
		m_versionNumber = 1.1;
	}
	else if (memcmp(header, version_1_2, 4) == 0) {
		m_bProbabilities = true;
		m_versionNumber = 1.2;
	}
	else if (memcmp(header, version_3_1, 4) == 0) {
	  m_bProbabilities = true;
	  m_versionNumber = 3.1;
	  m_infile.read((char *)&expNumSub, 4);
	  if (expNumSub != m_numSubjects) {
	    Rcpp::Rcerr << "Number of subjects in file doesn't agree with number of subjects indicated" << std::endl;
	    return 1;
	  }
	  m_headerSize = 12;
	}
	else {
		Rcpp::Rcerr << "File is not a binary dosage file" << std::endl;
		return 1;
	}
// Commented out because the code is OS specific
/*
	if (m_versionNumber == 1.1)
		expectedSize = (unsigned long long)(m_numSubjects)* m_numSNPs * 2 + 8;
	else
		expectedSize = (unsigned long long)(m_numSubjects)* m_numSNPs * 4 + 8;

	m_infile.seekg(0, ios_base::end);
	actualSize = FileSize(_filename.c_str());
//	cout << expectedSize << '\t' << actualSize << endl;
	if (actualSize != expectedSize) {
		cout << "Input file not of expected size" << endl;
		return 1;
	}
*/
  m_Dosage.set_size(m_numSubjects);
	
	if (m_versionNumber == 1.1)
		m_Probs.set_size(0, 0);
	else
	  m_Probs.set_size(m_numSubjects, 3);

	m_arraySize = m_numSubjects * sizeof(unsigned short);
	if (m_versionNumber == 1.2)
		m_arraySize += m_numSubjects * sizeof(unsigned short);
	
	m_SNPArray.resize((int)m_arraySize / 2);
	
	return 0;
}

// Read file and Map file
int CBinaryDosage::ReadFile(const std::string &_filename, unsigned int _numSub, const std::string &_mapFileName) {
	if (m_mapFile.ReadFile(_mapFileName)) {
		m_errorString = m_mapFile.ErrorString();
		return 1;
	}
  m_Skipped = m_mapFile.Skipped();
	m_numSubjects = _numSub;
	m_numSNPs = m_mapFile.NumSNPs();
	m_numSNPsUsed = m_mapFile.NumUsed();

	if (m_numSNPsUsed == 0) {
		m_errorString = "No SNPs to read";
		return 1;
	}

	return ReadFile(_filename);
}

int CBinaryDosage::ReOpen(const std::string &_filename, unsigned int _numSub, std::vector<bool> &_skipped,
           std::streampos _stp, unsigned int _currentSNP, double _versionNumber) {
  m_numSubjects = _numSub;
  m_Skipped = _skipped;
  m_numSNPs = m_Skipped.size();
  m_currentSNP = _currentSNP;
  m_versionNumber = _versionNumber;
  if (OpenFile(_filename, true))
    return 1;
  
  m_Dosage.set_size(m_numSubjects);
  
  if (m_versionNumber == 1.1) {
    m_Probs.set_size(0, 0);
    m_bProbabilities = false;
  } else {
    m_Probs.set_size(m_numSubjects, 3);
    m_bProbabilities = true;
  }

  m_arraySize = m_numSubjects * sizeof(unsigned short);
  if (m_versionNumber == 1.2)
    m_arraySize += m_numSubjects * sizeof(unsigned short);
  
  m_SNPArray.resize((int)m_arraySize / 2);

  m_infile.seekg(_stp);
  return 0;
}

// Process SNP if Version 1.1
void CBinaryDosage::ProcessSNP11() {
  m_Dosage = (arma::conv_to<arma::colvec>::from(m_SNPArray));
  m_Dosage.replace(65535, NA_REAL);
  m_Dosage /= 32767;
}

// Process SNP if Version 1.2
void CBinaryDosage::ProcessSNP12() {
  arma::Mat<unsigned short> rd(&m_SNPArray[0], 2, m_numSubjects, false, true);
	
	m_Probs.zeros();
	m_Probs.col(1) = arma::conv_to<arma::colvec>::from(rd.row(0));
	m_Probs.col(2) = arma::conv_to<arma::colvec>::from(rd.row(1));
	m_Probs.replace(65535, NA_REAL);
	m_Probs.col(0) = 65535 - (m_Probs.col(1) + m_Probs.col(2));
	m_Probs /= 65535.;
	m_Dosage = m_Probs.col(1) + m_Probs.col(2) + m_Probs.col(2);
}
// Process SNP if Version 3.1
void CBinaryDosage::ProcessSNP31() {
  unsigned int ui;
  unsigned short sp1;

  m_Probs.zeros();
  for (ui = 0; ui < m_numSubjects; ++ui) {
    if ((m_SNPArray(ui) & 0x8000) != 0) {
      m_SNPArray(ui) &= 0x7FFF;
      m_Dosage(ui) = m_SNPArray(ui) / 10000.;
      m_infile.read((char *)&sp1, 2);
      if ((sp1 & 0x8000) != 0) {
        sp1 &= 0x7FFF;
        m_Probs(ui, 1) = sp1 / 10000.;
        m_infile.read((char *)&sp1, 2);
        m_Probs(ui, 0) = sp1 / 10000.;
        m_infile.read((char *)&sp1, 2);
        m_Probs(ui, 2) = sp1 / 10000.;
      }
      else {
        m_Probs(ui, 1) = sp1 / 10000.;
        m_Probs(ui, 2) = (m_Dosage[ui] - m_Probs(ui, 1)) / 2;
        m_Probs(ui, 0) = 1 - m_Probs(ui, 1) - m_Probs(ui, 2);
      }
    }
    else {
      m_Dosage(ui) = m_SNPArray(ui) / 10000.;
      if (m_Dosage(ui) < 1) {
        m_Probs(ui, 1) = m_Dosage(ui);
        m_Probs(ui, 0) = 1 -  m_Probs(ui, 1);
        m_Probs(ui, 2) = 0;
      }
      else {
        m_Probs(ui, 2) = m_Dosage(ui) - 1;
        m_Probs(ui, 1) = 1 -  m_Probs(ui, 2);
        m_Probs(ui, 0) = 0;
      }
    }
  }
}
// Process the data read in for the last SNP
void CBinaryDosage::ProcessSNP() {
	if (m_versionNumber == 1.1)
		ProcessSNP11();
	else if (m_versionNumber == 1.2)
		ProcessSNP12();
	else
	  ProcessSNP31();
}