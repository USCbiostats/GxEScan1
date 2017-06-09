#include <Rcpp.h>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <string>
#include <cmath>
#include "MapFile.h"

using namespace Rcpp;

// ********************************************************
//                   CMapFile
// Class to read in map files associated with dosage files
// Has columns for Chromosome, SNP Name, location bp,
// location CM
// ********************************************************

// ***************************************************************************//
//                     Contstructor/Destructor                                //
// ***************************************************************************//

// Constructor
CMapFile::CMapFile() {
	m_numSNPs = 0;
	m_numUsed = 0;
}
// Destructor
CMapFile::~CMapFile() {
	m_infile.close();
}

// ***************************************************************************//
//                         Memory management                                  //
// ***************************************************************************//

// Resets values to defaults
// Deallocates memory in vectors
void CMapFile::Initialize() {
	std::vector<std::string>().swap(m_Chromosome);
	std::vector<std::string>().swap(m_SNPName);
	std::vector<double>().swap(m_GeneticDistance);
	std::vector<unsigned int>().swap(m_BasePairs);
	std::vector<bool>().swap(m_Skipped);
	m_numSNPs = 0;
	m_numUsed = 0;
	m_infile.close();
}
// Allocate memory needed to store data
int CMapFile::AllocateMemory() {
	if (m_numSNPs == 0) {
		m_errorString = "No SNPs to read in\n";
		Initialize();
		return 1;
	}
	m_Chromosome.resize(m_numSNPs, "");
	m_SNPName.resize(m_numSNPs, "");
	m_GeneticDistance.resize(m_numSNPs, 0);
	m_BasePairs.resize(m_numSNPs, 0);
	m_Skipped.resize(m_numSNPs, false);
	return 0;
}

// ***************************************************************************//
//                 Reading and validating routines                            //
// ***************************************************************************//

// Count the number of lines in the file and verify the entries
int CMapFile::CountRecords() {
	std::istringstream iss;
	std::ostringstream oss;
	std::string junk, strVal;

	m_infile.clear();
	m_infile.seekg(0);

	// Read in the lines
	m_numSNPs = 0;
	getline(m_infile, junk);
	while (m_infile.good() && !m_infile.eof()) {
		++m_numSNPs;
		iss.str(junk);
		iss.clear();
		iss.seekg(0);
		if (Validate(iss, junk)) // Make sure all entries are valid
			return 1;
		iss >> strVal; // This should fail - All entries have been read
		if (!iss.fail()) {
			oss.str("");
			oss << "Too many entries on line " << m_numSNPs << std::endl << junk << std::endl;
			m_errorString = oss.str();
			Initialize();
			return 1;
		}
		getline(m_infile, junk);
	}
	return 0;
}
// Validate the entries on the line read in
// The line read in is stored in inString
// It has already been associated with the istringstream
int CMapFile::Validate(std::istringstream &_iss, std::string &_inString) {
	std::string chr, snp, gdstr, bpstr;
	std::ostringstream oss;
	double gd;
	int chrNum, bp;
	bool missing;

	// There should be at least 4 entries per line
	_iss >> chr >> snp >> gdstr >> bpstr;
	if (_iss.fail()) {
		oss.str("");
		oss << "Error reading from line " << m_numSNPs << std::endl << _inString << std::endl;
		m_errorString = oss.str();
		Initialize();
		return 1;
	}
	// Check if chromosome is valid
	if(ASCII2Integer(chr, chrNum, missing)) { // Is chromosome entry a number?
		// If not
		if (!(chr == "X" || chr == "Y" || chr == "XY" || chr == "MT")) {
			oss.str("");
			oss << "Invalid entry for chromosome on line " << m_numSNPs << "   " << chr << std::endl << "Must be from 0 to 99 inclusive or X, Y, XY, or MT" << std::endl;
			m_errorString = oss.str();
			Initialize();
			return 1;
		}
	} else {	// chromosome entry is a number
		if (chrNum < 0 || chrNum > 99) {
			oss.str("");
			oss << "Invalid entry for chromosome on line " << m_numSNPs << "   " << chr << std::endl << "Must be from 0 to 99 inclusive or X, Y, XY, or MT" << std::endl;
			m_errorString = oss.str();
			Initialize();
			return 1;
		}
	}
	// SNP name is always valid - no need to check it
	// Check genetic distance is numeric or NA
	if (ASCII2Double(gdstr, gd, missing)) {
		oss.str("");
		oss << "Invalid entry for genetic distance on line " << m_numSNPs << "   " << gdstr << std::endl;
		m_errorString = oss.str();
		Initialize();
		return 1;
	}
	// Check base pairs is numeric or NA
	if (ASCII2Integer(bpstr, bp, missing)) {
		oss.str("");
		oss << "Invalid entry for base pairs on line " << m_numSNPs << "   " << bpstr << std::endl;
		m_errorString = oss.str();
		Initialize();
		return 1;
	}
	return 0;
}
// Read in the values from the string entered into the istringstream
// Error checking already done
void CMapFile::ReadValues(std::istringstream &_iss, int _n) {
	std::string gdstr, bpstr, chromosome, snpName;
	int bp;
	double gd;
	bool missing;

	_iss >> chromosome >> snpName >> gdstr >> bpstr;
	m_Chromosome[_n] = chromosome;
	m_SNPName[_n] = snpName;
	ASCII2Double(gdstr, gd, missing);
	m_GeneticDistance[_n] = gd;
	ASCII2Integer(bpstr, bp, missing);
	if (bp < 0) {
		m_Skipped[_n] = true;
		--m_numUsed;
		m_BasePairs[_n] = abs(bp);
	}
	else {
		m_BasePairs[_n] = bp;
	}
}
// No error should occur - error checking already completed
void CMapFile::ReadValues() {
	std::istringstream iss;
	std::string junk;
	unsigned int ui;

	m_infile.clear();
	m_infile.seekg(0);

	// Count the number of lines in the file
	// This is the number of SNPs
	m_numUsed = m_numSNPs;
	for (ui = 0; ui < m_numSNPs; ++ui) {
		getline(m_infile, junk); // Read in the string
		iss.str(junk);
		iss.clear();
		iss.seekg(0);
		ReadValues(iss, ui); // Process the string
	}
}
int CMapFile::ReadFile(const std::string &filename) {
	Initialize();	// Clear out any previous read data
	if (OpenFile(filename)) // Try to open the file
		return 1;

	// Count number of records
	if (CountRecords())
		return 1;

	// Allocate memory
	if (AllocateMemory())
		return 1;

	// Read in the values
	// No errors should occur. Error checking already done
	ReadValues();

	m_infile.close();
	return 0;
}

// ********************************************************
//               CBimFile
// Class to read in map file associated with a binary
// genetic file.
// Has same information as a CMapFile but also has allele values
// Derived from CMapFile
// ********************************************************

// ***************************************************************************//
//                     Contstructor/Destructor                                //
// ***************************************************************************//

// Constructor
CBimFile::CBimFile() : CMapFile() {
}
// Destructor
CBimFile::~CBimFile() {
}

// ***************************************************************************//
//                         Memory management                                  //
// ***************************************************************************//

// Deallocates memory and resets variables to default values - Calls Initialize of CMapFile
void CBimFile::Initialize() {
  std::vector<std::string>().swap(m_FirstAllele);
  std::vector<std::string>().swap(m_SecondAllele);
  
	CMapFile::Initialize();
}
// Allocate memory - Calls AllocateMemory from CMapFile
int CBimFile::AllocateMemory() {
	if (CMapFile::AllocateMemory())
		return 1;
	m_FirstAllele.resize(m_numSNPs, "");
	m_SecondAllele.resize(m_numSNPs, "");
	return 0;
}

// ***************************************************************************//
//                 Reading and validating routines                            //
// ***************************************************************************//

// Validate the entries - Calls Validate of CMapFile
int CBimFile::Validate(std::istringstream &_iss, std::string &_inString) {
	std::ostringstream oss;

	std::string allele;
	if (CMapFile::Validate(_iss, _inString))
		return 1;
	_iss >> allele >> allele;
	if (_iss.fail()) {
		oss.str("");
		oss << "Error reading from line number " << m_numSNPs << std::endl << _inString << std::endl;
		m_errorString = oss.str();
		Initialize();
		return 1;
	}
	return 0;
}
// Read the values - Calls ReadValues of CMapFile
void CBimFile::ReadValues(std::istringstream &_iss, const int _n) {
	CMapFile::ReadValues(_iss, _n);
	_iss >> m_FirstAllele[_n] >> m_SecondAllele[_n];
}

// ***************************************************************************//
//                  Routine to swap allele names                              //
// ***************************************************************************//

// Swap allele values - this is only for one SNP
void CBimFile::SwapAlleles(unsigned int n) {
	std::string swap;

	if (n < m_numSNPs) {
		swap = m_FirstAllele[n];
		m_FirstAllele[n] = m_SecondAllele[n];
		m_SecondAllele[n] = swap;
	}
}