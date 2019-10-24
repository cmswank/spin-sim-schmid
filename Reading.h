
#ifndef _READING_H_
#define _READING_H_

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
using namespace std;

/**
 * Helper class to interface with BOOST spirit regexp class
 */
class Reading {

 private:

  ifstream *file;

  int varLength;

  char varName[100]; // max length 100 chars

  double varVector[16];

  char expression[2048];

  char tempLine[8192];

  bool parse_assignment(char const* str, std::string &varname, vector<double>& v);

  bool parse_string_assignment(char const* str, std::string &varname, 
			       std::string &expression);

  
  bool parse_comment(char const* str);


  // store all parameters in vectors:
  int nParams; // contain the total number of assignments
 
  vector<std::string> names;  // the names of the parameters
  vector<int> nVars; // the number of arguments
  vector<vector<double> > vars; // the args (N.B. the 2D structure)
  


 public:

 
  Reading();
  
  int open (const char *filename);

  int getNextAssignment();
  
  char *getVarName() {return this->varName;};
  char *getExpression() {return this->expression;};
  const char *getVarName(int i) {return this->names[i].c_str();};

  int getVarLength() {return this->varLength;};
  int getVarLength(int i) {return this->nVars[i];};
  
  double *getVarVector() {return this->varVector;};
  vector<double> getVarVector(int i) {return this->vars[i];};

  char *getLineInput() {return this->tempLine;};

  int getEntries () {return this->nParams;};


  
  int parseFile();

};



#endif
