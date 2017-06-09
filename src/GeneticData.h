#ifndef GENETICDATA_H
#define GENETICDATA_H 1
#ifndef MAPFILE_H
#include "MapFile.h"
#endif 

// Base class for reading all genetic data files
class CGeneticData : public CReadFile {
protected:
    // Number of SNPs that are used
	unsigned int m_numSNPs;
	// Number of subjects in genetic dataset - can be different than the family file
	unsigned int m_numSubjects;
	// SNP current in memory
	unsigned int m_currentSNP;
	// Are the values measured or imputed
	bool m_bMeasured;
	// Do we have the probabilities or just the dosage - imputed only
	bool m_bProbabilities;
	// Array to hold dosage and probabilities
	// Changed to aramadillo vector and matrix
	arma::Col <double> m_Dosage;
	arma::Mat <double> m_Probs;
	// Swap space need to swap probability values if alleles are swapped - used to set major/minor allele
	// Fix - No longer needed (?)
	// double *m_pSwapSpace;
public:
	CGeneticData();
	virtual ~CGeneticData();

	virtual void SwapAlleles();

	bool Measured() const { return m_bMeasured; }
	bool Probabilities() const { return m_bProbabilities; }
	const arma::Col <double> &Dosage() const { return m_Dosage; }
	const arma::Mat <double> &Probs() const { return m_Probs; }
	const std::streampos FilePosition() { return m_infile.tellg(); }
	int CurrentSNP() const { return m_currentSNP; }
	//	const std::vector <bool> &Missing() const { return m_Missing; }
	// This is the number of SNPs used
	virtual unsigned int NumSNPs() const { return m_numSNPs; }
	// This is the number of SNPs found in the file - default value is the number of SNPs in map file
	virtual unsigned int NumSNPsInFile() const { return m_numSNPs; }
	// This is the number of subjects in the family file
	virtual unsigned int NumSubjects() const { return m_numSubjects; }
	// This is the number of subjects in the genetic file - default value is number of subjects in the family file
	virtual unsigned int NumSubjectsInFile() const { return m_numSubjects; }
	virtual unsigned int NumSubjectsWithData() const { return m_numSubjects; }
		
	virtual const std::string &SNPName() const = 0;
	virtual const std::string &FirstAllele() const = 0;
	virtual const std::string &SecondAllele() const = 0;
	virtual const std::string &Chromosome() const = 0;
	virtual unsigned int Location() const = 0;
	virtual int GetFirst() = 0;
	virtual int GetNext() = 0;
};

// Base class for reading binary genetic data files
class CBinaryGeneticData : public CGeneticData {
protected:
	// Binary map file - need for binary genetic data
	CBimFile m_mapFile;
  // Pointer to vector containing indicators is SNP is used
  std::vector<bool> m_Skipped;
	// Size of file header
	unsigned int m_headerSize;
	// Size of data array for one SNP
	unsigned int m_arraySize;
	// Array of unprocessed SNP data
	// Fix - May be changed to type arma::Col <short> would also rename m_SNPArray
//	char *m_pSNPArray;
  arma::Col <unsigned short> m_SNPArray;
public:
	CBinaryGeneticData();
	~CBinaryGeneticData();

	virtual void SwapAlleles();

	virtual const std::string &SNPName() const { return m_mapFile.SNP()[m_currentSNP]; }
	virtual const std::string &FirstAllele() const { return m_mapFile.FirstAllele()[m_currentSNP]; }
	virtual const std::string &SecondAllele() const { return m_mapFile.SecondAllele()[m_currentSNP]; }
	virtual const std::string &Chromosome() const { return m_mapFile.Chromosome()[m_currentSNP]; }
	virtual int CurrentSNP() const { return m_currentSNP; }
	virtual unsigned int Location() const { return m_mapFile.BasePairs()[m_currentSNP]; }

	virtual int GetFirst();
	virtual int GetNext();
	virtual void GoTo(std::streampos p) { m_infile.seekg(p); }

	int ReadSNP();
	virtual void ProcessSNP() = 0;
	const CBimFile &MapFile() const { return m_mapFile; }
};
#endif

