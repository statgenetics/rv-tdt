//!\file gw_utilities.cpp
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

#include "gw_utilities.h"
bool is_file_empty (std::string filename) 
{
  std::ifstream file;
  file.open(filename.c_str());	
  if (!file) 
    return true;

  std::string line;
  UINT length = 0;  
  while ( getline(file, line) ) 
    ++length;
  if (length > 0) 
    return false;
  else
    return true;
}

void scan_vector2F(std::string filename, vector2F& info) 
{
  //!- scan file into vector2F
  //std::cout << "\tSource data: "<< filename << std::endl;
  std::ifstream file;
  file.open(filename.c_str());	
  if (!file) {
    std::cerr << "\tUnable to source [ " << filename << " ]. Quit." << std::endl;
    exit(-1); // terminate with error
  }
  std::string line;

  while ( getline(file, line) ) {
    vectorF data;
    double value;
    std::istringstream iss(line);
    while (iss >> value) data.push_back(value);
    info.push_back(data);
  }
  return;
}

void scan_vector2UI(std::string filename, vector2UI& info) 
{
  //!- scan file into vector2F
  //std::cout << "\tSource data: "<< filename << std::endl;
  std::ifstream file;
  file.open(filename.c_str());	
  if (!file) {
    std::cerr << "\tUnable to source [ " << filename << " ]. Quit." << std::endl;
    exit(-1); // terminate with error
  }
  std::string line;

  while ( getline(file, line) ) {
    vectorUI data;
    UINT value;
    std::istringstream iss(line);
    while (iss >> value) data.push_back(value);
    info.push_back(data);
  }
  return;
}

//make a directory
std::string makedir(std::string dirname)
{
  //frist check if exists
  struct stat st;
  if(stat(dirname.c_str(),&st) == 0){
    //wait for a second
    RNG rng; gsl_rng* gslr; gslr = rng.get();
    unsigned msecond = gsl_rng_uniform_int(gslr, 1000);
    usleep(msecond*1000);
  }
  if(stat(dirname.c_str(),&st) == 0){;}
  //then mkdir
  else{
    if (mkdir(dirname.c_str(),0755)==-1){
      //struct tm *current; time_t now;
      //time(&now);
      //current = localtime(&now);
      //
      //printf("At time %i:%i:%i: ", current->tm_hour, current->tm_min, current->tm_sec);
      //std::cout << "Unable to make directory named as " << dirname << std::endl << "Quit Now" << std::endl;
      //exit(-1);
      dirname += "_bak";
      int dummy = mkdir(dirname.c_str(),0755);
      dummy++;
    }
  }
  // return final dirname
  return dirname;
}
