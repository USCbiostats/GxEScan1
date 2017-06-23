#ifndef BINARYDOSAGE_H
#define BINARYDOSAGE_H 1
#endif
#ifndef GENETICDATA_H
#include "GeneticData.h"
#endif

int VCF_to_BinaryDosage(std::string vcfFilename, std::string outBaseFilename, unsigned int initSub = 10000);

Rcpp::List ExtractDosages(std::string bdosageFilename, std::string mapFilename, unsigned int numSNPs);

Rcpp::List ExtractMoreDosages(Rcpp::List inputs);

class CBinaryDosage : public CBinaryGeneticData {
protected:
	unsigned int m_numSNPsUsed;
	double m_versionNumber;

	void ProcessSNP11();
	void ProcessSNP12();
	void ProcessSNP31();
public:
	CBinaryDosage();
	~CBinaryDosage();

	virtual int ReadFile(const std::string &_filename);
	int ReadFile(const std::string &_filename, unsigned int _numSub, const std::string &_mapFileName);
	int ReOpen(const std::string &_filename, unsigned int _numSub, std::vector<bool> &_skipped,
            std::streampos _stp, unsigned int _currentSNP, double _versionNumber);

	double Version() const { return m_versionNumber; }
	virtual unsigned int NumSNPs() const { return m_numSNPsUsed; }
	virtual void ProcessSNP();
};