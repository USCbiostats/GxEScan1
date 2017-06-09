#ifndef MAPFILE_H
#define MAPFILE_H 1
#ifndef READFILE_H
#include "ReadFile.h"
#endif
#ifndef NUMBERSTRINGS_H
#include "NumberStrings.h"
#endif

class CMapFile : public CReadFile {
protected:
	// Number of SNPs in the file
	unsigned int m_numSNPs;
	// Chromosome value 0-99 or X, Y, XY, or MT
	std::vector<std::string> m_Chromosome;
	// SNP name - usually rs#
	std::vector<std::string> m_SNPName;
	// Genetic distance in morgans -- Not used by GxEScan
	std::vector<double> m_GeneticDistance;
	// Base pairs from 3 prime end, I think - only used in output
	std::vector<unsigned int> m_BasePairs;
	// Skip SNP indicator - occurs when base pair value is less than 0
	std::vector<bool> m_Skipped;
	// Number of SNPs used, i.e., not skipped
	unsigned int m_numUsed;

	// Clear allocate memory and prepare to read new file
	virtual void Initialize();

	// Validate the values in the current string
	virtual int Validate(std::istringstream &_iss, std::string &_inString);
	// Read in the values and store them
	virtual void ReadValues(std::istringstream &_iss, const int _n);

	// Allocate memory for stored values
	virtual int AllocateMemory();
	// Count entries in file
	int CountRecords();
	// Read the entries
	void ReadValues();
public:
	CMapFile();
	~CMapFile();

	// Read the data from the file
	virtual int ReadFile(const std::string &_filename);

	unsigned int NumSNPs() const { return m_numSNPs; }

	const std::vector<std::string> &Chromosome() const { return m_Chromosome; }
	const std::vector<std::string> &SNP() const { return m_SNPName; }
	const std::vector<double> &GeneticDistance() const { return m_GeneticDistance; }
	const std::vector<unsigned int> &BasePairs()	const { return m_BasePairs; }
	const std::vector<bool> &Skipped() const { return m_Skipped; }
	unsigned int NumUsed() const { return m_numUsed;  }
};

// The .bim file is of the same format as a .map file.
// The only difference is a .bim file has two additional columns
class CBimFile : public CMapFile {
protected:
	// First allele listed in the file
	std::vector<std::string> m_FirstAllele;
	// Second allele listed in the file
	std::vector<std::string> m_SecondAllele;

	// Deallocate memory and set variables to default values
	virtual void Initialize();

	// Read the data checking for valid values
	// Does not store values, only verfies entries are valid
	virtual int Validate(std::istringstream &_iss, std::string &_inString);
	// Read in the values and store them
	virtual void ReadValues(std::istringstream &_iss, const int _n);
	// Allocates memory needed to store read in values
	virtual int AllocateMemory();
public:
	CBimFile();
	~CBimFile();

	// Swap first and second allele values
	void SwapAlleles(unsigned int _n);

	const std::vector<std::string> FirstAllele() const { return m_FirstAllele; }
	const std::vector<std::string> SecondAllele() const { return m_SecondAllele; }
};

#endif
