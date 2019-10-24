#include "Reading.h"

#ifndef __CINT__

#include <boost/version.hpp>

#if BOOST_VERSION >= 103800
 #define BOOST_SPIRIT_USE_OLD_NAMESPACE
 #include "boost/spirit/include/classic_core.hpp"
 #include "boost/spirit/include/classic_push_back_actor.hpp"
 #include "boost/spirit/include/classic_assign_actor.hpp"
#else
 #include "boost/spirit/core.hpp"
 #include "boost/spirit/actor/push_back_actor.hpp"
 #include "boost/spirit/actor/assign_actor.hpp"
#endif

#endif


///////////////////////////////////////////////////////////////////////////////


/**
 * Constructor
 */
Reading::Reading() {
  this->nParams = 0;
  this->file = 0;
}


///////////////////////////////////////////////////////////////////////////////
//
//  Parser for assignments with optional comments
//
///////////////////////////////////////////////////////////////////////////////

/**
 * parse simple assignment
 */
bool
Reading::parse_assignment(char const* str, std::string &varname, 
			  vector<double>& v) {
  
  using namespace std;
  using namespace boost::spirit;
  
    
    return parse(str,

        //  Begin grammar
        (
	 (+alnum_p)[assign_a(varname)] >> ch_p('=') >>
	 real_p[push_back_a(v)] >> *(',' >> real_p[push_back_a(v)]) >> ch_p(';') >> 
	 !(ch_p('/') >> ch_p('/') >> *(anychar_p))      // optional comment
	 >> end_p
        )
        ,
        //  End grammar

        space_p).full;
}

/**
 * parse string assignment
 */
bool
Reading::parse_string_assignment(char const* str, std::string &varname, 
				 std::string &expression) {
  
  using namespace std;
  using namespace boost::spirit;
  
    
    return parse(str,

        //  Begin grammar
        (
	 (+alnum_p)[assign_a(varname)] >> ch_p('=') >>
	 ch_p('"') >> (*(~ch_p('"')))[assign_a(expression)] >> ch_p('"') >> ch_p(';') >> 
	 !(ch_p('/') >> ch_p('/') >> *(anychar_p))      // optional comment
	 >> end_p
        )
        ,
        //  End grammar

        space_p).full;
}


/**
 * parse comment or empty line
 */
bool
Reading::parse_comment(char const* str)
{
  using namespace std;
    using namespace boost::spirit;

    return parse(str,

        //  Begin grammar
        (
	 !(ch_p('/') >> ch_p('/') >> *(anychar_p))
	 >> end_p
        )
        //  End grammar
		 
        ,space_p).full;
}

////////////////////////////////////////////////////////////////////////////
//
//  Main program
//
////////////////////////////////////////////////////////////////////////////

/**
 * Open file to be read
 */
int Reading::open (const char *filename)
{

  file = new ifstream();

  file->open(filename);
  if (!file->is_open()) {
    cout << "Error opening file " << filename <<  "\n";
    return -1;
  }

  return 0;
}

/**
 * Get next assignment from par file
 */
int Reading::getNextAssignment() {

  int i;

  std::string varName;
  std::string expression;

  vector<double> v;

  using namespace std;
  using namespace boost::spirit;

  if (!file->good())
    return 0;

  while (file->good()) {
    file->getline(tempLine, 8192);

    if (parse_assignment(tempLine, varName, v)) {
      /* cout << "ASSIGN : " << varName.c_str() << " --> ";
      for (i = 0; i < v.size(); i++) {
	cout << v[i] << " ";
      }
      cout << endl;
      */
      sprintf(this->varName, "%s", varName.c_str());
      this->varLength = v.size();
      for (i=0; i < (int)v.size() && i < 16; i++)
	this->varVector[i] = v[i];
      
      if (this->varLength > 16) {
	this->varLength = 16;
	cout << "Warning: trying to assign too big of a vector, part of it lost\n";
      }
      
      // save in vectors;
      this->names.push_back(varName);
      this->nVars.push_back(this->varLength);
      this->vars.push_back(v);
      nParams++;

      return 1;
    }
    else if (parse_string_assignment(tempLine, varName, expression)) {
      sprintf(this->varName, "%s", varName.c_str());
      sprintf(this->expression, "%s", expression.c_str());

      return 2;
    }
    // comments and empty lines
    else if (parse_comment(tempLine));
    else {
      cout << "ERROR : " << tempLine << endl;
      cout << "\t...ignored\n";
    }

  }

  return 0;
}

/**
 * Parse the parameter file
 */
int Reading::parseFile() {

  int nParams = 0;

  if (!file->good())
    return 0;

  while(this->getNextAssignment()) {
    nParams++;
  }
  return nParams;
}
