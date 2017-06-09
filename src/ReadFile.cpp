#include <Rcpp.h>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <cstring>
#include "ReadFile.h"

using namespace Rcpp;

// *****************************************************************//
//                    Constructor/Destructor                        //
// *****************************************************************//

// Constructor
CReadFile::CReadFile() {
	m_errorString = "";
}
// Destructor
CReadFile::~CReadFile() {
	m_infile.close(); // Should already be closed
}

// *****************************************************************//
//                           Open File                              //
// *****************************************************************//

// Open the file
// Returns an error if the file can't be opened
int CReadFile::OpenFile(const std::string &_filename, bool _binaryFile) {
	m_infile.close(); // Should not be open
	m_infile.clear();	// Clear error from closing if not open
	if (_binaryFile)
		m_infile.open(_filename.c_str(), std::ios::in | std::ios::binary);
	else
		m_infile.open(_filename.c_str());
	// Check if file successfully opened
	if (!m_infile.good()) {
		m_errorString = "Unable to open file " + _filename + "\n";
		return 1;
	}
	return 0;
}


