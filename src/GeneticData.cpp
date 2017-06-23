#include <RcppArmadillo.h>
#include <vector>
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include "GeneticData.h"

using namespace Rcpp;

//*******************************************************************************
//                      CGeneticData
// Base class for reading files with genetic data
//*******************************************************************************

// ***************************************************************************//
//                     Contstructor/Destructor                                //
// ***************************************************************************//

// Constructor
CGeneticData::CGeneticData() : CReadFile() {
	m_bMeasured = true;
	m_bProbabilities = false;
	m_numSNPs = 0;
	m_numSubjects = 0;
	m_currentSNP = 0;;
}
// Destructor
CGeneticData::~CGeneticData() {
	m_infile.close();
	m_infile.clear();
}

// ***************************************************************************//
//                         Allele Swapping                                    //
// ***************************************************************************//

// Swap the alleles - changes the dosage and swaps P(g=0) with P(g=2)
// Fix - m_pSwapSpace used to be used here. RcppArmadillo has a function
// to swap column values. m_pSwapSpace can be removed form the code entirely (?)
void CGeneticData::SwapAlleles() {
	m_Dosage = 2 - m_Dosage;
	if (m_bMeasured == false && m_bProbabilities == true)
	  m_Probs.swap_cols(0,2);
}

//*******************************************************************************
//                      CGeneticData
// Base class for reading binary files with genetic data
//*******************************************************************************

// ***************************************************************************//
//                     Contstructor/Destructor                                //
// ***************************************************************************//

// Constructor
CBinaryGeneticData::CBinaryGeneticData() : CGeneticData() {
	m_headerSize = 0;
	m_arraySize = 0;
	// Fix - m_pSNPArray may be changed to a std::vector
//	m_pSNPArray = NULL;
}
// Destructor
CBinaryGeneticData::~CBinaryGeneticData() {
}

// ***************************************************************************//
//                         Allele Swapping                                    //
// ***************************************************************************//

// Swap the alleles - changes the dosage and swaps P(g=0) with P(g=2) and allele names
// Changed because of armadillo functions
void CBinaryGeneticData::SwapAlleles() {
	CGeneticData::SwapAlleles();
	m_mapFile.SwapAlleles(m_currentSNP);
}

// ***************************************************************************//
//                         Reading routines                                   //
// ***************************************************************************//

// Reads SNP data from file for current record
int CBinaryGeneticData::ReadSNP() {
  // Are there more SNPs to read?
  if (m_currentSNP >= m_numSNPs)
    return 1;
  // Skipped SNPs indicated as unused in map file
  while (m_Skipped[m_currentSNP]) {
		++m_currentSNP;
		if (m_currentSNP >= m_numSNPs)
			return 1;
// Need to update to different class for each format
// Current version was done to quickly implement version 3.1
//		m_infile.seekg(m_arraySize, std::ios::cur);
		m_infile.read((char *)&m_SNPArray[0], m_arraySize);
		ProcessSNP();
	}
  // Fix - if m_pSNPArray is changed to a std::vector the following line
  // needs to be changed.
  // Would (char *)&m_SNPArray[0] work here? 
	m_infile.read((char *)&m_SNPArray[0], m_arraySize);
  // Fix - m_Missing may be dropped
//	std::fill(m_Missing.begin(), m_Missing.end(), false);
	ProcessSNP();
	return 0;
}
// Get the first SNP - resets file to first record and reads SNP
int CBinaryGeneticData::GetFirst() {
  m_currentSNP = 0;
	m_infile.clear();
	m_infile.seekg(m_headerSize);
	return ReadSNP();
}
// Get the next SNP - moves file to next record and reads SNP
int CBinaryGeneticData::GetNext() {
	++m_currentSNP;
	return ReadSNP();
}
