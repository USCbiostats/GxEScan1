#ifndef NUMBERSTRINGS_H
#define NUMBERSTRINGS_H 1

// ******************************** ASCII to Integer functions *************************
int ASCII2Integer(const std::string &_intString, int &_n);
int ASCII2Integer(const std::string &_intString, int &_n, bool &_missing);
int ASCII2Integer(const std::string &_intString, const int _missingValue, int &_n, bool &_missing);
int ASCII2Integer(const std::string &_intString, const std::string &_missingValue, int &_n, bool &_missing);

// ******************************** ASCII to Double functions *************************
int ASCII2Double(const std::string &_doubleString, double &_x);
int ASCII2Double(const std::string &_doubleString, double &_x, bool &_missing);
int ASCII2Double(const std::string &_doubleString, const double _missingValue, double &_x, bool &_missing);
int ASCII2Double(const std::string &_doubleString, const std::string &_missingValue, double &_x, bool &_missing);

// **************************** string trimming function ******************************
std::string Trim(const std::string &_str);

// ***************************** Process Integer String ********************************
int ProcessNumberString(unsigned int _sz, bool *_nlist, unsigned int &_numEntry, const std::string &_ns, std::string &_errorString);

#endif