#include <Rcpp.h>
#include <cstdlib>
#include <sstream>
#include <cstring>
#include <string>
#include "NumberStrings.h"

using namespace Rcpp;

// ******************************** ASCII to Integer functions *************************

// Reads an ASCII string and determines if it is an integer
// Blanks spaces are trimmed of the end of the string.
// A return value of 0 indicates the string was successfully read in as an integer
// otherwise a value of 1 is returned. If an error is returned
// the value for the integer is set to 0.
int ASCII2Integer(const std::string &_intString, int &_n) {
	std::istringstream iss;
	std::string trimString;

	// White space is trimmed off the end of the string
	trimString = Trim(_intString);

	// The complete string must be read
	// If there are characters left over
	// there is an error. Those characters
	// are non numeric, i.e. not 0,1,...,9
	iss.str(trimString);
	iss >> _n;
	if (iss.fail() || !iss.eof()) {
		_n = 0;
		return 1;
	}
	return 0;
}

// Reads an ASCII string and determines if it is an integer
// Blanks spaces are trimmed off the end of the string.
// A return value of 0 indicates the string was successfully read in as an integer
// otherwise a value of 1 is returned. If an error is returned
// the value for the integer is set to 0.
// An error indicator is also returned
int ASCII2Integer(const std::string &_intString, int &_n, bool &_error) {
	std::istringstream iss;
	std::string trimString;

	// White space is trimmed off the end of the string
	trimString = Trim(_intString);

	// The complete string must be read
	// If there are characters left over
	// there is an error.
	iss.str(trimString);
	iss >> _n;
	if (iss.fail() || !iss.eof()) {
		_n = 0;
		_error = true;
		return 1;
	}
	_error = false;
	return 0;
}

// Reads an ASCII string and determines if it is an integer.
// Blanks spaces are trimmed off the end of the string.
// A return value of 0 indicates the string was successfully read in as an integer
// otherwise a value of 1 is returned. If an error is returned
// the value for the integer is set to 0.
// The read in value is compared to a missing value passed to the routine.
// If the read in value is equal to the missing value, the missing indicator is set
// and the integer value is set to 0.
int ASCII2Integer(const std::string &_intString, const int _missingValue, int &_n, bool &_missing) {
	if (ASCII2Integer(_intString, _n, _missing))
		return 1;
	if (_n == _missingValue) {
		_n = 0;
		_missing = true;
	}
	return 0;
}

// Reads an ASCII string and determines if it is an integer.
// Blanks spaces are trimmed off the end of the string.
// A return value of 0 indicates the string was successfully read in as an integer
// otherwise a value of 1 is returned. If an error is returned
// the value for the integer is set to 0.
// The entered string is compared to a missing value string passed to the routine.
// If the entered string is equal to the missing string, the missing indicator is set
// and the integer value is set to 0.
int ASCII2Integer(const std::string &_intString, const std::string &_missingValue, int &_n, bool &_missing) {
	if (Trim(_intString) == _missingValue) {
		_n = 0;
		_missing = true;
		return 0;
	}
	return ASCII2Integer(_intString, _n, _missing);
}

// ******************************** ASCII to Double functions *************************

// Reads an ASCII string and determines if it can be stored as a double
// Blanks spaces are trimmed of the end of the string.
// A return value of 0 indicates the string was successfully read in as a double
// otherwise a value of 1 is returned. If an error is returned
// the value for the double is set to 0.
int ASCII2Double(const std::string &_doubleString, double &_x) {
	std::istringstream iss;
	std::string trimString;

	// White space is trimmed off the end of the string
	trimString = Trim(_doubleString);

	iss.str(trimString);
	iss >> _x;
	// The complete string must be read
	// If there are characters left over
	// there is an error.
	if (iss.fail() || !iss.eof()) {
		_x = 0.;
		return 1;
	}
	return 0;
}

// Reads an ASCII string and determines if it is a double
// Blanks spaces are trimmed off the end of the string.
// A return value of 0 indicates the string was successfully read in as an double
// otherwise a value of 1 is returned. If an error is returned
// the value for the double is set to 0.
// An error indicator is also returned
int ASCII2Double(const std::string &_doubleString, double &_x, bool &_error) {
	std::istringstream iss;
	std::string trimString;

	// White space is trimmed off the end of the string
	trimString = Trim(_doubleString);

	iss.str(trimString);
	iss >> _x;
	// The complete string must be read
	// If there are characters left over
	// there is an error.
	if (iss.fail() || !iss.eof()) {
		_x = 0.;
		_error = true;
		return 1;
	}
	_error = false;
	return 0;
}

// Reads an ASCII string and determines if it is a double.
// Blanks spaces are trimmed off the end of the string.
// A return value of 0 indicates the string was successfully read in as a double
// otherwise a value of 1 is returned. If an error is returned
// the value for the double is set to 0.
// The read in value is compared to a missing value passed to the routine.
// If the read in value is equal to the missing value, the missing indicator is set
// and the double value is set to 0.
int ASCII2Double(const std::string &_doubleString, const double _missingValue, double &_x, bool &_missing) {
	if (ASCII2Double(_doubleString, _x, _missing))
		return 1;
	if (_x == _missingValue) {
		_x = 0.;
		_missing = true;
	}
	return 0;
}

// Reads an ASCII string and determines if it is a double.
// Blanks spaces are trimmed off the end of the string.
// A return value of 0 indicates the string was successfully read in as a double
// otherwise a value of 1 is returned. If an error is returned
// the value for the double is set to 0.
// The entered string is compared to a missing value string passed to the routine.
// If the entered string is equal to the missing string, the missing indicator is set
// and the double value is set to 0.
int ASCII2Double(const std::string &_doubleString, const std::string &_missingValue, double &_x, bool &_missing) {
	if (Trim(_doubleString) == _missingValue) {
		_x = 0;
		_missing = true;
		return 0;
	}
	return ASCII2Double(_doubleString, _x, _missing);
}

// ****************** string trimming function ******************************

// Trims the trailing spaces off a string
// This modifies the passed string
std::string Trim(const std::string &_str) {
	std::string trimmed;
	size_t lastSpace;

	lastSpace = _str.find_last_not_of(" ");
	if (lastSpace == std::string::npos) {
		trimmed = _str;
	}
	else {
		trimmed = _str.substr(0, lastSpace + 1);
	}
	return trimmed;

}

// ************************* Process Integer String **************************

// Processes a comma delimited integer number string and assigns indicators for each
// A range of numbers can indicated using the -, e.g. 2-4 indicates to use 2, 3, and 4
// The numbers must be postive, and the maximum value is specified by sz
// nlist is the array of indicators returned. It must already be declared and have size sz
// The number of values indicated by the string is returned in numEntry
// ns is the string to process
// If an error occurs a value of 1 is returned and an error message is written to the error string
// otherwise a value of 0 is returned.
int ProcessNumberString(const unsigned int _sz, bool *_nlist, unsigned int &_numEntry, const std::string &_ns, std::string &_errorString) {
	std::istringstream iss, iss2;
	std::string junk;
	unsigned int x, y;
	unsigned int ui;
	unsigned int mx;
	size_t loc;

	_numEntry = 0;

	mx = 0;
	memset(_nlist, false, _sz * sizeof(bool));
	// Check the string for errors and find the maximum value
	iss.str(_ns);
	getline(iss, junk, ',');
	while (!iss.fail()) {
		loc = junk.find('-');
		if (loc == std::string::npos) {
			iss2.str(junk);
			iss2.clear();
			iss2 >> x;
			if (iss2.fail() || !iss2.eof() || x == 0) {
				_errorString = "Invalid integer string ";
				_errorString += _ns;
				return 1;
			}
			if (x > mx)
				mx = x;
			else {
				_errorString = "Numbers in integer string must be in ascending order";
				return 1;
			}
		}
		else if (loc == 0 || junk.find('-', loc + 1) != std::string::npos) {
			_errorString = "Error reading integer string ";
			_errorString += _ns;
			return 1;
		}
		else {
			junk[loc] = ' ';
			iss2.str(junk);
			iss2.clear();
			iss2 >> x >> y;
			if (iss2.fail() || !iss2.eof() || x == 0 || y == 0) {
				_errorString = "Invalid integer string ";
				_errorString += _ns;
				return 1;
			}
			if (y <= x || x <= mx) {
				_errorString = "Numbers in integer string must be in ascending order";
				return 1;
			}
			if (y > mx) {
				mx = y;
			}
		}
		getline(iss, junk, ',');
	}
	if (mx > _sz) {
		_errorString = "Maximum value in integer list is larger than maximum value allowed";
		return 1;
	}
	// Read the string and store indicator array
	// No error check since that was already done
	iss.clear();
	iss.seekg(0);
	getline(iss, junk, ',');
	while (!iss.fail()) {
		loc = junk.find('-');
		if (loc == std::string::npos) {
			iss2.str(junk);
			iss2.clear();
			iss2.seekg(0);
			iss2 >> x;
			_nlist[x - 1] = true;
			++_numEntry;
		}
		else {
			junk[loc] = ' ';
			iss2.str(junk);
			iss2.clear();
			iss2.seekg(0);
			iss2 >> x >> y;
			for (ui = x - 1; ui < y; ++ui) {
				_nlist[ui] = true;
				++_numEntry;
			}
		}
		getline(iss, junk, ',');
	}
	return 0;
}

