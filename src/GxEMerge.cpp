#include <cstring>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <Rcpp.h>
//#include "OS_routines.h"
//#include "NumberStrings.h"

const std::string LogisticFileExtension[8] = { "CC_GxE", "CC_2df", "CC_3df", "CC_DG", "CC_GE", "CC_DGGE", "Case_GE", "Cntl_GE" };
const std::string LinearFileExtension[5] = { "QT_YG", "QT_GxE", "QT_2df", "QT_VH", "QT_YGVH" };
const std::string StatColumnTitle[6] = { "Z", "CHISQ_1DF", "CHISQ_2DF", "CHISQ_3DF", "CHISQ_4DF", "F" };
const std::string UnderscoreString = "_";

class CCombineFiles {
protected:
	char *m_fileContents;
	const char **m_fileNames;
	bool m_bLinear;
	size_t m_numFiles;
	size_t *m_numSNPs;
	size_t m_totalSNPs;

//	int ReadListFile(const std::string &_filename);
	int ProcessSNPFiles(const Rcpp::StringVector &_basefileNames, const std::string &_filename);
//	int ProcessLogisticFiles(const std::string &_outFilename);
//	int ProcessLinearFiles(const std::string &_outFilename);
//	int Process5out(const std::string &_filename, const std::string &_extension, unsigned int _statNo);
//	int Process4out(const std::string &_filename, const std::string &_extension, unsigned int _statNo);
//	int Process3out(const std::string &_filename, const std::string &_extension, unsigned int _statNo);
public:
	CCombineFiles();
	~CCombineFiles();

//	void SetLinear(bool _bLinear = true) { m_bLinear = _bLinear; }
	int ProcessList(bool _bLinear, const Rcpp::StringVector &_basefileNames, const std::string &_outFilename);
};

CCombineFiles::CCombineFiles() {
	m_fileContents = NULL;
	m_fileNames = NULL;
	m_bLinear = false;
	m_numFiles = 0;
	m_numSNPs = NULL;
	m_totalSNPs = 0;
}
CCombineFiles::~CCombineFiles() {
	if (m_fileContents)
		delete[] m_fileContents;
	if (m_fileNames)
		delete[] m_fileNames;
	if (m_numSNPs)
		delete[] m_numSNPs;
}
/*
int CCombineFiles::ReadListFile(const std::string &_filename) {
	long long int filesize;
	unsigned int ui;
	char *lineEnd;

	std::ifstream infile;

	if (m_fileContents)
		delete[] m_fileContents;
	if (m_fileNames)
		delete[] m_fileNames;
	m_numFiles = 0;

	infile.open(_filename.c_str());
	if (!infile.good()) {
		Rcpp::Rcerr << "Unable to open file " << _filename << std::endl;
		return 1;
	}

	filesize = FileSize(_filename.c_str());
	m_fileContents = new char[(unsigned int)filesize + 1];
	memset(m_fileContents, 0, (size_t)filesize + 1);
	infile.read(m_fileContents, filesize);
	infile.close();
	m_numFiles = std::count(m_fileContents, m_fileContents + filesize, 0x0a);

	m_fileNames = new const char*[m_numFiles];
	m_fileNames[0] = m_fileContents;
	lineEnd = m_fileContents;
	for (ui = 1; ui < m_numFiles; ++ui) {
		lineEnd = strchr(lineEnd, 0x0a);
		*lineEnd = 0;

		++lineEnd;
		m_fileNames[ui] = lineEnd;
	}
	lineEnd = strchr(lineEnd, 0x0a);
	*lineEnd = 0;

	Rcpp::Rcout << "File size:\t\t" << filesize << std::endl;
	Rcpp::Rcout << "Number of file names in list:\t" << m_numFiles << std::endl;
//	for (ui = 0; ui < m_numFiles; ++ui)
//		Rcpp::Rcout << m_fileNames[ui] << std::endl;

	return 0;
}
 */
int CCombineFiles::ProcessSNPFiles(const Rcpp::StringVector &_basefileNames, const std::string &_outFilename) {
	std::string snpFilename;
  std::string header;
	std::ifstream snpInfile;
	std::ofstream snpOutfile;
	int ui; // uj;
	unsigned int offset;
//	size_t currentFileSize, maxFileSize;
//	char *fileContents = NULL;
//	const char *contentsStart = NULL;
	std::string outData;
	unsigned int snpid;
	std::string chr, snp, bp, a1;
	std::istringstream iss;

	if (m_numSNPs)
		delete m_numSNPs;
	m_numSNPs = NULL;
  m_numFiles = _basefileNames.length();
	// Check that the files exist
	for (ui = 0; ui < _basefileNames.length(); ++ui) {
		snpFilename = _basefileNames[ui];
		snpFilename += ".snpinfo";
		snpInfile.open(snpFilename.c_str());
		if (!snpInfile.good()) {
			Rcpp::Rcerr << "Unable to open file " << snpFilename << std::endl;
			return 1;
		}
		snpInfile.close();
	}
	
	snpFilename = _outFilename + ".snpinfo";
	snpOutfile.open(snpFilename.c_str());
	if (!snpOutfile.good()) {
		Rcpp::Rcerr << "Unable to open output file " << snpFilename << std::endl;
		return 1;
	}

	m_numSNPs = new size_t[m_numFiles];
	memset(m_numSNPs, 0, m_numFiles * sizeof(size_t));
	
	snpOutfile << "SNPID\tCHR\tSNP\tBP\tA1" << std::endl;
	offset = 0;
	for (ui = 0; ui < _basefileNames.length(); ++ui) {
		snpFilename = _basefileNames[ui];
		snpFilename += ".snpinfo";
		snpInfile.open(snpFilename.c_str());
    getline(snpInfile, header);
    
    snpInfile >> snpid >> chr >> snp >> bp >> a1;
    while(snpInfile.good()) {
      if (snpid > m_numSNPs[ui])
        m_numSNPs[ui] = snpid;
      snpOutfile << offset + snpid << '\t' << chr << '\t' << snp << '\t' << bp << '\t' << a1 << std::endl;
      snpInfile >> snpid >> chr >> snp >> bp >> a1;
    }
		offset += m_numSNPs[ui];
		snpInfile.close();
	}
	m_totalSNPs = offset;
	snpOutfile.close();

	return 0;
}
/*
int CCombineFiles::ProcessLinearFiles(const std::string &_outFilename){
	if (Process5out(_outFilename, LinearFileExtension[0], 0))
		return 1;
	Rcpp::Rcout << LinearFileExtension[0] << " done" << std::endl;
	if (Process5out(_outFilename, LinearFileExtension[1], 0))
		return 1;
	Rcpp::Rcout << LinearFileExtension[1] << " done" << std::endl;
	if (Process4out(_outFilename, LinearFileExtension[2], 2))
		return 1;
	Rcpp::Rcout << LinearFileExtension[2] << " done" << std::endl;
	if (Process3out(_outFilename, LinearFileExtension[3], 5))
		return 1;
	Rcpp::Rcout << LinearFileExtension[3] << " done" << std::endl;
	if (Process3out(_outFilename, LinearFileExtension[4], 4))
		return 1;
	Rcpp::Rcout << LinearFileExtension[4] << " done" << std::endl;
	return 0;
}
int CCombineFiles::ProcessLogisticFiles(const std::string &_outFilename){
	if (Process5out(_outFilename, LogisticFileExtension[0], 0))
		return 1;
	Rcpp::Rcout << LogisticFileExtension[0] << " done" << std::endl;
	if (Process4out(_outFilename, LogisticFileExtension[1], 2))
		return 1;
	Rcpp::Rcout << LogisticFileExtension[1] << " done" << std::endl;
	if (Process4out(_outFilename, LogisticFileExtension[2], 3))
		return 1;
	Rcpp::Rcout << LogisticFileExtension[2] << " done" << std::endl;
	if (Process5out(_outFilename, LogisticFileExtension[3], 0))
		return 1;
	Rcpp::Rcout << LogisticFileExtension[3] << " done" << std::endl;
	if (Process3out(_outFilename, LogisticFileExtension[4], 1))
		return 1;
	Rcpp::Rcout << LogisticFileExtension[4] << " done" << std::endl;
	if (Process3out(_outFilename, LogisticFileExtension[5], 2))
		return 1;
	Rcpp::Rcout << LogisticFileExtension[5] << " done" << std::endl;
	if (Process5out(_outFilename, LogisticFileExtension[6], 0))
		return 1;
	Rcpp::Rcout << LogisticFileExtension[6] << " done" << std::endl;
	if (Process5out(_outFilename, LogisticFileExtension[7], 0))
		return 1;
	Rcpp::Rcout << LogisticFileExtension[7] << " done" << std::endl;
	return 0;
}
int CCombineFiles::Process5out(const std::string &_outFilename, const std::string &_extension, unsigned int _statNo) {
	std::string statFilename;
	std::ifstream *statInfile;
	std::ofstream statOutfile;
	unsigned int ui, uj;
	std::string outData;
	unsigned int *offset;
	unsigned int *fileCount;
	unsigned int *snpid;
	std::string *nmiss, *beta, *stat, *pstr;
	double *p;
	unsigned int lowest;

	statInfile = new std::ifstream[m_numFiles];
	for (ui = 0; ui < m_numFiles; ++ui) {
		statFilename = m_fileNames[ui] + UnderscoreString + _extension;
		statFilename += ".gxeout";
		statInfile[ui].open(statFilename.c_str());
		if (!statInfile[ui].good()) {
			Rcpp::Rcerr << "Unable to open file " << statFilename << std::endl;
			delete[] statInfile;
			return 1;
		}
	}

	statFilename = _outFilename + UnderscoreString + _extension + ".gxeout";
	statOutfile.open(statFilename.c_str());
	if (!statOutfile.good()) {
		Rcpp::Rcerr << "Unable to open output file " << statFilename << std::endl;
		return 1;
	}

	offset = new unsigned int[m_numFiles];
	fileCount = new unsigned int[m_numFiles];
	snpid = new unsigned int[m_numFiles];
	nmiss = new std::string[m_numFiles];
	beta = new std::string[m_numFiles];
	stat = new std::string[m_numFiles];
	pstr = new std::string[m_numFiles];
	p = new double[m_numFiles];

	offset[0] = 0;
	for (ui = 1; ui < (unsigned int)m_numFiles; ++ui)
		offset[ui] = offset[ui - 1] + m_numSNPs[ui - 1];
	memset(fileCount, 0, m_numFiles*sizeof(unsigned int));

	for (ui = 0; ui < m_numFiles; ++ui) {
		statInfile[ui] >> nmiss[ui] >> nmiss[ui] >> beta[ui] >> stat[ui] >> pstr[ui];
		statInfile[ui] >> snpid[ui] >> nmiss[ui] >> beta[ui] >> stat[ui] >> pstr[ui];
		if (pstr[ui] == "NA")
			p[ui] = 9;
		else
			ASCII2Double(pstr[ui], p[ui]);
	}

	statOutfile << "SNPID\tNMISS\tBETA\t" << StatColumnTitle[_statNo] << "\tP" << std::endl;

	for (ui = 0; ui < m_totalSNPs; ++ui) {
//	for (ui = 0; ui < 20; ++ui) {
		for (uj = 0; uj < m_numFiles; ++uj) {
			if (p[uj] != -9)
				break;
		}
		lowest = uj;
		for (; uj < m_numFiles; ++uj) {
			if (p[uj] != -9 && p[uj] < p[lowest])
				lowest = uj;
		}
		statOutfile << snpid[lowest] + offset[lowest] << '\t' << nmiss[lowest] << '\t' << beta[lowest] << '\t' << stat[lowest] << '\t' << pstr[lowest] << std::endl;
		++fileCount[lowest];
		if (fileCount[lowest] < m_numSNPs[lowest]) {
			statInfile[lowest] >> snpid[lowest] >> nmiss[lowest] >> beta[lowest] >> stat[lowest] >> pstr[lowest];
			if (pstr[lowest] == "NA")
				p[lowest] = 9;
			else
				ASCII2Double(pstr[lowest], p[lowest]);
		}
		else {
			p[lowest] = -9;
		}
	}
	for (ui = 0; ui < m_numFiles; ++ui)
		statInfile[ui].close();

	statOutfile.close();
	if (statInfile)
		delete[] statInfile;
	if (offset)
		delete[] offset;
	if (fileCount)
		delete[] fileCount;
	if (snpid)
		delete[] snpid;
	if (nmiss)
		delete[] nmiss;
	if (beta)
		delete[] beta;
	if (stat)
		delete[] stat;
	if (pstr)
		delete[] pstr;
	if (p)
		delete[] p;
	return 0;

}
int CCombineFiles::Process4out(const std::string &_outFilename, const std::string &_extension, unsigned int _statNo) {
	std::string statFilename;
	std::ifstream *statInfile;
	std::ofstream statOutfile;
	unsigned int ui, uj;
	std::string outData;
	unsigned int *offset;
	unsigned int *fileCount;
	unsigned int *snpid;
	std::string *nmiss, *stat, *pstr;
	double *p;
	unsigned int lowest;

	statInfile = new std::ifstream[m_numFiles];
	for (ui = 0; ui < m_numFiles; ++ui) {
		statFilename = m_fileNames[ui] + UnderscoreString + _extension;
		statFilename += ".gxeout";
		statInfile[ui].open(statFilename.c_str());
		if (!statInfile[ui].good()) {
			Rcpp::Rcerr << "Unable to open file " << statFilename << std::endl;
			delete[] statInfile;
			return 1;
		}
	}

	statFilename = _outFilename + UnderscoreString + _extension + ".gxeout";
	statOutfile.open(statFilename.c_str());
	if (!statOutfile.good()) {
		Rcpp::Rcerr << "Unable to open output file " << statFilename << std::endl;
		return 1;
	}

	offset = new unsigned int[m_numFiles];
	fileCount = new unsigned int[m_numFiles];
	snpid = new unsigned int[m_numFiles];
	nmiss = new std::string[m_numFiles];
	stat = new std::string[m_numFiles];
	pstr = new std::string[m_numFiles];
	p = new double[m_numFiles];

	offset[0] = 0;
	for (ui = 1; ui < (unsigned int)m_numFiles; ++ui)
		offset[ui] = offset[ui - 1] + m_numSNPs[ui - 1];
	memset(fileCount, 0, m_numFiles*sizeof(unsigned int));

	for (ui = 0; ui < m_numFiles; ++ui) {
		statInfile[ui] >> nmiss[ui] >> nmiss[ui] >> stat[ui] >> pstr[ui];
		statInfile[ui] >> snpid[ui] >> nmiss[ui] >> stat[ui] >> pstr[ui];
		if (pstr[ui] == "NA")
			p[ui] = 9;
		else
			ASCII2Double(pstr[ui], p[ui]);
	}

	statOutfile << "SNPID\tNMISS\t" << StatColumnTitle[_statNo] << "\tP" << std::endl;

	for (ui = 0; ui < m_totalSNPs; ++ui) {
		//	for (ui = 0; ui < 20; ++ui) {
		for (uj = 0; uj < m_numFiles; ++uj) {
			if (p[uj] != -9)
				break;
		}
		lowest = uj;
		for (; uj < m_numFiles; ++uj) {
			if (p[uj] != -9 && p[uj] < p[lowest])
				lowest = uj;
		}
		statOutfile << snpid[lowest] + offset[lowest] << '\t' << nmiss[lowest] << '\t' << stat[lowest] << '\t' << pstr[lowest] << std::endl;
		++fileCount[lowest];
		if (fileCount[lowest] < m_numSNPs[lowest]) {
			statInfile[lowest] >> snpid[lowest] >> nmiss[lowest] >> stat[lowest] >> pstr[lowest];
			if (pstr[lowest] == "NA")
				p[lowest] = 9;
			else
				ASCII2Double(pstr[lowest], p[lowest]);
		}
		else {
			p[lowest] = -9;
		}
	}
	for (ui = 0; ui < m_numFiles; ++ui)
		statInfile[ui].close();

	statOutfile.close();
	if (statInfile)
		delete[] statInfile;
	if (offset)
		delete[] offset;
	if (fileCount)
		delete[] fileCount;
	if (snpid)
		delete[] snpid;
	if (nmiss)
		delete[] nmiss;
	if (stat)
		delete[] stat;
	if (pstr)
		delete[] pstr;
	if (p)
		delete[] p;
	return 0;

}
int CCombineFiles::Process3out(const std::string &_outFilename, const std::string &_extension, unsigned int _statNo) {
	std::string statFilename;
	std::ifstream *statInfile;
	std::ofstream statOutfile;
	unsigned int ui, uj;
	std::string outData;
	unsigned int *offset;
	unsigned int *fileCount;
	unsigned int *snpid;
	std::string *stat, *pstr;
	double *p;
	unsigned int lowest;

	statInfile = new std::ifstream[m_numFiles];
	for (ui = 0; ui < m_numFiles; ++ui) {
		statFilename = m_fileNames[ui] + UnderscoreString + _extension;
		statFilename += ".gxeout";
		statInfile[ui].open(statFilename.c_str());
		if (!statInfile[ui].good()) {
			Rcpp::Rcerr << "Unable to open file " << statFilename << std::endl;
			delete[] statInfile;
			return 1;
		}
	}

	statFilename = _outFilename + UnderscoreString + _extension + ".gxeout";
	statOutfile.open(statFilename.c_str());
	if (!statOutfile.good()) {
		Rcpp::Rcerr << "Unable to open output file " << statFilename << std::endl;
		return 1;
	}

	offset = new unsigned int[m_numFiles];
	fileCount = new unsigned int[m_numFiles];
	snpid = new unsigned int[m_numFiles];
	stat = new std::string[m_numFiles];
	pstr = new std::string[m_numFiles];
	p = new double[m_numFiles];

	offset[0] = 0;
	for (ui = 1; ui < (unsigned int)m_numFiles; ++ui)
		offset[ui] = offset[ui - 1] + m_numSNPs[ui - 1];
	memset(fileCount, 0, m_numFiles*sizeof(unsigned int));

	for (ui = 0; ui < m_numFiles; ++ui) {
		statInfile[ui] >> stat[ui] >> stat[ui] >> pstr[ui];
		statInfile[ui] >> snpid[ui] >> stat[ui] >> pstr[ui];
		if (pstr[ui] == "NA")
			p[ui] = 9;
		else
			ASCII2Double(pstr[ui], p[ui]);
	}

	statOutfile << "SNPID\t" << StatColumnTitle[_statNo] << "\tP" << std::endl;

	for (ui = 0; ui < m_totalSNPs; ++ui) {
		//	for (ui = 0; ui < 20; ++ui) {
		for (uj = 0; uj < m_numFiles; ++uj) {
			if (p[uj] != -9)
				break;
		}
		lowest = uj;
		for (; uj < m_numFiles; ++uj) {
			if (p[uj] != -9 && p[uj] < p[lowest])
				lowest = uj;
		}
		statOutfile << snpid[lowest] + offset[lowest] << '\t' << stat[lowest] << '\t' << pstr[lowest] << std::endl;
		++fileCount[lowest];
		if (fileCount[lowest] < m_numSNPs[lowest]) {
			statInfile[lowest] >> snpid[lowest] >> stat[lowest] >> pstr[lowest];
			if (pstr[lowest] == "NA")
				p[lowest] = 9;
			else
				ASCII2Double(pstr[lowest], p[lowest]);
		}
		else {
			p[lowest] = -9;
		}
	}
	for (ui = 0; ui < m_numFiles; ++ui)
		statInfile[ui].close();

	statOutfile.close();
	if (statInfile)
		delete[] statInfile;
	if (offset)
		delete[] offset;
	if (fileCount)
		delete[] fileCount;
	if (snpid)
		delete[] snpid;
	if (stat)
		delete[] stat;
	if (pstr)
		delete[] pstr;
	if (p)
		delete[] p;
	return 0;

}
*/
int CCombineFiles::ProcessList(bool _bLogistic, const Rcpp::StringVector &_basefileNames, const std::string &_outFilename) {
	if (ProcessSNPFiles(_basefileNames, _outFilename))
		return 1;
	Rcpp::Rcout << "SNP file completed" << std::endl;
/*
	if (_bLogistic == true)
		return ProcessLinearFiles(_outFilename);
	return ProcessLogisticFiles(_outFilename);
*/
  return 0;
}

//' Function to merge the results from several GxEScans
//' 
//' Function to merge the results from several GxEScans for use with GxEResults
//' to produce a summary for a genome wide scan
//' 
//' @param logistic
//' Indicator if the outcome was logistic, otherwise outcome was linear
//' @param basefileNames
//' List of base file names of results
//' @param outfileName
//' Base filename for the merged results
//' @importFrom Rcpp evalCpp
//' @useDynLib GxEScanR
//' @export
// [[Rcpp::export]]
int GxEMerge(bool logistic, Rcpp::StringVector basefileNames, std::string outfileName) {
	CCombineFiles combineFiles;

  Rcpp::Rcout << outfileName << "\n";
  for (int i = 0; i < basefileNames.size(); ++i)
    Rcpp::Rcout << basefileNames[i] << "\n";
  
  return (combineFiles.ProcessList(logistic, basefileNames, outfileName));
//  return (0);
}
