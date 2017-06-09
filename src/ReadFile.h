#ifndef READFILE_H
#define READFILE_H 1

// Base class for reading plink files
class CReadFile {
protected:
	std::string m_errorString;
	std::ifstream m_infile;
	virtual int OpenFile(const std::string &_filename, bool _binaryFile = false);

public:
	CReadFile();
	virtual ~CReadFile();

	virtual int ReadFile(const std::string &_filename) = 0;
	const std::string &ErrorString() const { return m_errorString; }
};

#endif