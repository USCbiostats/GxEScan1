#ifndef BEDFILE_H
#define BEDFILE_H 1
#ifndef GENETICDATA_H
#include "GeneticData.h"
#endif

// Class to read a .fam file in plink format
class CBedFile : public CBinaryGeneticData {
protected:
	void Initialize();
	void ClearArrays();
public:
	CBedFile();
	virtual ~CBedFile();

	virtual int ReadFile(const std::string &filename);
	int ReadFile(const std::string &filename, unsigned int numSub, const std::string &mapFilename);

	virtual void ProcessSNP();
};

#endif
