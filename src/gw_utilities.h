//!\file gw_utilities.h
//!\brief useful utilities
// Copyright 2010 2011 Gao Wang

/*  *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *
 *                                                                          *
 *  This program is free software: you can redistribute it and/or modify    *
 *  it under the terms of the GNU General Public License as published by    *
 *  the Free Software Foundation, either version 3 of the License, or       *
 *  (at your option) any later version.                                     *
 *                                                                          *
 *  This program is distributed in the hope that it will be useful,         *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *  GNU General Public License for more details.                            *
 *                                                                          *
 *  You should have received a copy of the GNU General Public License       *
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *                                                                          *
 *  *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   *   */

#ifndef GWUTILITY_H
#define GWUTILITY_H

#include <limits>
#include <iomanip>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <ctime>
#include <cstdlib>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <map>
#include <ctype.h>
#include <sys/types.h>
#include <dirent.h>
#include <stdlib.h>
#include "gsl/gsl_rng.h"
#include "gsl/gsl_cdf.h"
#include <sys/stat.h>
#include <json/json.h>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/join.hpp>
#include <boost/math/special_functions/fpclassify.hpp>
#include <boost/foreach.hpp>

typedef std::vector<double> vectorF;
typedef std::vector< std::vector<double> > vector2F;
typedef std::vector< std::vector< std::vector <double> > > vector3F;

typedef std::vector <char> vectorC;
typedef std::vector< std::vector <char> > vector2C;

typedef std::vector<int> vectorI;
typedef std::vector< std::vector <int> > vector2I;

typedef unsigned int UINT;
typedef std::vector<unsigned int> vectorUI;
typedef std::vector< std::vector<unsigned int> > vector2UI;
typedef std::vector< std::vector< std::vector <unsigned int> > > vector3UI;

typedef std::vector<bool> vectorL;
typedef std::vector< std::vector<bool> > vector2L;

typedef std::pair <std::string, double> pairSD;
typedef std::vector< std::pair <std::string, double> > vectorPairSD;

//!\brief setup GNU/RNG
class RNG {

  public:
  //!- the second-generation ranlux generators have the strongest proof of randomness.
  //!- However turns out gsl_rng_alloc(gsl_rng_ranlxs2) is buggy !!
  //gsl_rng_mt19937
    RNG() : rng( gsl_rng_alloc(gsl_rng_mt19937) ) {};
    ~RNG() { gsl_rng_free( rng ); }
    gsl_rng* get()
    {
      // time(NULL): number of seconds since 00:00:00 GMT Jan. 1, 1970
      __seed = static_cast<unsigned long>(time (NULL) + getpid());
      gsl_rng_set(rng, __seed);
      return rng;
    }
    gsl_rng* get(const unsigned long seed)
    {
      __seed = seed;
      gsl_rng_set(rng, __seed);
      return rng;
    }

  private:
    gsl_rng* rng;
    unsigned long __seed;
};


namespace std {
//!- Dump a vector to screen
template<class T> ostream & operator<<(ostream & out, const vector<T> & vec)
  {
    if (!vec.empty()) {
      typename vector<T>::const_iterator it = vec.begin();
      out << setiosflags(ios::fixed) << setprecision(4) << setw(7) << *it;
      for (++it; it != vec.end(); ++it)
        out << " " << setiosflags(ios::fixed) << setprecision(4) << setw(7) <<  *it ;
    }
    return out;
  }
}

template<class T> Json::Value buildJsonArray (const std::vector<T> & vec)
{
    Json::Value jVec(Json::arrayValue);
    if (!vec.empty()) {
        typename std::vector<T>::const_iterator it = vec.begin();
        for( ; it != vec.end(); ++it){
            jVec.append(Json::Value(*it));
        }
    }
    return jVec;
}

template<class T> inline bool scan_2dvector(std::string filename, std::vector< std::vector<T> > &vec)
{
	std::ifstream file;
	file.open(filename.c_str());
	if(!file) return false;

    vec.resize(0); //clear up
    std::string line;
	while(getline(file, line)){
        boost::trim(line); //get rid of whitespace on ends
        if(boost::starts_with(line, "#") || line.empty()) continue; //comments or blank

		std::vector<T> data;
		T value;
		std::istringstream iss(line);
		while(iss >> value) data.push_back(value);
		vec.push_back(data);
	}
	return true;
}

//make a directory, if not exists, wait a second for other potential process; if exists, return
inline bool makedir(std::string dirname)
{
	struct stat st;
	unsigned int count(0); // at most wait 10 cycles
	RNG rng; gsl_rng* gslr; gslr = rng.get();
	if(stat(dirname.c_str(), &st)!=0) usleep(gsl_rng_uniform_int(gslr, 1000)*1000); // waiting
    else return true;

	// then need to create a folder
	if(stat(dirname.c_str(),&st)!=0){
		while(mkdir(dirname.c_str(),0755)==-1 && count<10){
			usleep(gsl_rng_uniform_int(gslr, 1000)*1000); // waiting again
			count++;
		}
	}
	return count>=10?false:true;
}

inline bool is_file_empty(std::string filename)
{
	std::ifstream file;
	file.open(filename.c_str());
    if(!file) return true;

	std::string line;
	UINT length = 0;
	while(getline(file, line)) ++length;
	if(length>0) return false;
	else return true;
}

#endif ///:~
