#include <RcppArmadillo.h>
#include <cstdlib>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include "BedFile.h"

// Fix - Instead of -1 for missing, use value Rcpp considers missing
// If this works, m_Missing can be removed from other sections of the code
// I think I found the value, NA_REAL
const double BinGeno[4] = { 2., NA_REAL, 1., 0 };

//*******************************************************************************
//                      CBedFile
// Class to read in binary genetic files in plink format, .bed
//*******************************************************************************
// Constructor
CBedFile::CBedFile() : CBinaryGeneticData() {
	m_headerSize = 3;
}
// Destructor
CBedFile::~CBedFile() {
}

// Initialize values - this may be legacy code ???
void CBedFile::Initialize() {
	m_numSNPs = 0;
	m_numSubjects = 0;
	m_arraySize = 0;
	m_currentSNP = 0;
	// Initialize shouldn't be called if there is an open file.
	m_infile.close();
}
// Clears out the allocated memory and resets values to defaults - may be legacy code ???
void CBedFile::ClearArrays() {
	m_infile.close();
	m_infile.clear();
	m_Dosage.reset();
	m_Probs.reset();
	// Fix - May be possible to remove m_Missing
//	std::vector<bool>().swap(m_Missing);
}

// Reads in map file and assigns to binary genetic data file
int CBedFile::ReadFile(const std::string &filename, unsigned int numSub, const std::string &mapFilename) {
	if (m_mapFile.ReadFile(mapFilename)) {
		m_errorString = m_mapFile.ErrorString();
		return 1;
	}
	m_Skipped = m_mapFile.Skipped();
	m_numSubjects = numSub;
	m_numSNPs = m_mapFile.NumSNPs();

	return ReadFile(filename);
}
// Reads in binary genetic data file. The complete file is read and stored -- I may want to revist this and only read data as required ???
// This has been updated - Only one SNP is in memory at a time. GetFirst() and GetNext() read SNP data in.
int CBedFile::ReadFile(const std::string &filename) {
	int magic = 0;
	char mode = 0;
	std::ostringstream oss;

	if (OpenFile(filename, true))
		return 1;

	if (m_numSubjects <= 0 || m_numSNPs <= 0) {
		m_errorString += "Number of subjects and number of SNPs must be greater than 0.\n";
		m_infile.close();
		return 1;
	}

	m_infile.read((char *)&magic, 2);
	m_infile.read((char *)&mode, 1);

	if (magic != 0x1B6C) {
		m_errorString += filename + " is not a plink binary genetic file" + '\n';
		m_infile.close();
		return 1;
	}
	if (mode != 1) {
		m_errorString += filename + " is not in SNP-major format" + '\n';
		m_infile.close();
		return 1;
	}

	m_arraySize = (m_numSubjects + 3)/4;
	// Fix m_pSNPArray may be replaced with an armadillo vector of type short
	// This would require using (m_arraySize + 1) / 2 for the vector size
	// Id est, m_SNPArray.resize((m_arraySize + 1) / 2);
	// Want to keep m_arraySize because that is the number of bytes to read in
	m_SNPArray.resize((int)(m_arraySize + 1) / 2);
	m_Dosage.set_size(m_numSubjects);
	m_Dosage.zeros();
	// Fix - m_Missing may be removed
//	m_Missing.resize(m_numSubjects, false);

	return 0;
}

// Process the data for the current record and assign SNP values
void CBedFile::ProcessSNP() {
	const char *x;
	char y;
	unsigned int ui;

	x = (char *)&m_SNPArray[0];
	y = *x;
	for (ui = 0; ui < m_numSubjects;) {
		m_Dosage[ui] = BinGeno[y & 3];
	  // Fix - if missing value can be used, next two lines can be removed
//		if (m_Dosage[ui] < 0)
//			m_Missing[ui] = true;
		++ui;
		if (ui % 4) {
			y = y >> 2;
		}
		else {
			++x;
			y = *x;
		}
	}
}