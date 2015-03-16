#include "gw_utilities.h"
#include "tdt_tools.h"
#include <unistd.h>

//For unphased data
bool scan_unphasedTrioFile(std::string filename, vector2F& info)
{
  //!- scan file into vector2F
  //std::cout << "\tSource data: "<< filename << std::endl;
  std::ifstream file;
  file.open(filename.c_str());
  if (!file) {
    return false; // return error
  }
  std::string line;

  //erease first row (head line) and first two columns
  getline(file, line);
  while ( getline(file, line) ) {
    vectorF data;
    double value;
    std::istringstream iss(line);
    while (iss >> value) data.push_back(value);
    data.erase(data.begin(), data.begin()+2);
    info.push_back(data);
  }

  return true;
}

bool scan_unphasedAnnFile(std::string filename,std::vector<std::string>& position, std::vector<std::string>& mutType, vectorF& samMaf, vectorF& popMaf, double deMafs)
{
  int chrom(0), pos(1), colSam(4), colPop(5), mutOne(7), mutTwo(8);
  //!- scan file into vector2F
  //std::cout << "\tSource data: "<< filename << std::endl;
  std::ifstream file;
  file.open(filename.c_str());
  if (!file) {
    return false; // return error
  }
  std::string line;
  std::string chr;
  std::vector<std::string> location;
  // ignore first row (headline)
  getline(file, line);
  while ( getline(file, line) ) {
    std::vector<std::string> data;
    std::string value;
    std::istringstream iss(line);
    while (iss >> value) data.push_back(value);
    chr = data[chrom].c_str();
    location.push_back(data[pos].c_str());
    samMaf.push_back(atof(data[colSam].c_str()));

    std::string pop = data[colPop];
    // replace "NA" with 1/SampleSize
    if(pop == "NA") {popMaf.push_back(deMafs);}
    else
    {
      if(atof(pop.c_str()) == 0) {popMaf.push_back(deMafs);}
      else popMaf.push_back(atof(pop.c_str()));
    }

    mutType.push_back(data[mutOne] + "_" + data[mutTwo]);
  }
  // vector basic includes: chrom, location, mutType
  position.push_back(chr);
  position.push_back(location[0] + "-" + *(--location.end()));
  return true;
}


//For phased data
bool scan_phasedTrioFile(std::string filename, vector2F& info, std::vector<std::string>& makers)
{
  //!- scan file into vector2F
  //std::cout << "\tSource data: "<< filename << std::endl;
  std::ifstream file;
  file.open(filename.c_str());
  if (!file) {
    return false; // return error
  }
  std::string line;
  //second line is maker information
  getline(file, line);
  std::string pos;
  std::istringstream iss(line);
  while (iss >> pos) makers.push_back(pos);
  while ( getline(file, line) ) {
    vectorF data;
    double value;
    std::istringstream iss(line);
    while (iss >> value) data.push_back(value);
    info.push_back(data);
  }
  return true;
}

bool scan_phasedAnnFile(std::string filename, std::vector<std::string>& position, std::vector<std::string>& mutType, vectorF& samMaf, vectorF& popMaf, double deMafs)
{
  int chrom(0), pos(1), colSam(5), colPop(6), mutOne(2), mutTwo(3);
  //!- scan file into vector2F
  //std::cout << "\tSource data: "<< filename << std::endl;
  std::ifstream file;
  file.open(filename.c_str());
  if (!file) {
    return false; // return error
  }
  std::string line;
  std::string chr;
  std::vector<std::string> location;
  // ignore first row (headline)
  getline(file, line);
  while ( getline(file, line) ) {
    std::vector<std::string> data;
    std::string value;
    std::istringstream iss(line);
    while (iss >> value) data.push_back(value);
    chr = data[chrom].c_str();
    location.push_back(data[pos].c_str());
    samMaf.push_back(atof(data[colSam].c_str()));

    std::string pop = data[colPop];
    // replace "NA" with 1/SampleSize
    if(pop == "NA") {popMaf.push_back(deMafs);}
    else
    {
      if(atof(pop.c_str()) == 0) {popMaf.push_back(deMafs);}
      else popMaf.push_back(atof(pop.c_str()));
    }

    mutType.push_back(data[mutOne] + "_" + data[mutTwo]);
  }
  // vector basic includes: chrom, location, mutType
  position.push_back(chr);
  position.push_back(location[0] + "-" + *(--location.end()));
  for(UINT i=0; i<location.size(); i++) position.push_back(location[i]);
  return true;
}

//For indel data
bool scan_IndelTrioFile(std::string filename, vector2F& final, std::vector<std::string>& makers)
{
  vector2F info;
  //!- scan file into vector2F
  //std::cout << "\tSource data: "<< filename << std::endl;
  std::ifstream file;
  file.open(filename.c_str());
  //file.seekg(0, std::ios::end); // put the "cursor" at the end of the file
  //int length = file.tellg(); //find the position of the cursor
  if (!file) {
    return false; // return error
  }
  std::string line;
  //erease first row (head line)
  getline(file, line);
  //second line is maker information
  getline(file, line);
  std::string pos;
  std::istringstream iss(line);
  while (iss >> pos) makers.push_back(pos);
  while ( getline(file, line) ) {
    vectorF data;
    double value;
    std::istringstream iss(line);
    while (iss >> value) data.push_back(value);
    info.push_back(data);
  }

  //make it unphased
  for(UINT i=0; i < info.size()/2; i++)
  {
    vectorF oneind;
    for(UINT j=0; j<info[2*i].size(); j++)
    {
      if(info[2*i][j] == -9 || info[2*i+1][j] == -9) oneind.push_back(-9);
      else oneind.push_back(info[2*i][j] + info[2*i + 1][j]);
    }
    final.push_back(oneind);
  }
  return true;
}

bool scan_IndelAnnFile(std::string filename, std::vector<std::string>& position, vectorF& samMaf, double deMafs)
{
  int chrom(0), pos(1), colSam(4);
  //!- scan file into vector2F
  //std::cout << "\tSource data: "<< filename << std::endl;
  std::ifstream file;
  file.open(filename.c_str());
  if (!file) {
    return false; // return error
  }
  std::string line;
  std::string chr;
  std::vector<std::string> location;
  // ignore first row (headline)
  getline(file, line);
  while ( getline(file, line) ) {
    std::vector<std::string> data;
    std::string value;
    std::istringstream iss(line);
    while (iss >> value) data.push_back(value);
    chr = data[chrom].c_str();
    location.push_back(data[pos].c_str());

    std::string sam = data[colSam];
    // replace "NA" with 1/SampleSize
    if(sam == "NA") {samMaf.push_back(deMafs);}
    else
    {
      if(atof(sam.c_str()) < deMafs) {samMaf.push_back(deMafs);}
      else samMaf.push_back(atof(sam.c_str()));
    }
  }
  // vector basic includes: chrom, location, mutType
  position.push_back(chr);
  position.push_back(location[0] + "-" + *(--location.end()));
  for(UINT i=0; i<location.size(); i++) position.push_back(location[i]);
  return true;
}


// get the list of the files in folder dir, ignoring the files which has a '.' at the first char of file name
void getdir(std::string dir, std::vector<std::string> &files)
{
    DIR *dp;
    struct dirent *dirp;
    if((dp  = opendir(dir.c_str())) == NULL) {
        std::cout << "Error: opening " << dir << " . Quit Now!" << std::endl;
        exit(-1);
    }

    while ((dirp = readdir(dp)) != NULL) {
	//get ride of . and ..
	//only find *.dat file
	if( int(std::string(dirp->d_name).find(".dat")) > 0 )
	{
	    std::string geneDat = std::string(dirp->d_name);
	    //only keep gene name, erase .dat
	    geneDat.erase(geneDat.find(".dat"), 4);
            files.push_back(geneDat);
	}
    }
    closedir(dp);
    return;
}

////make a directory
//std::string makedir(std::string dirname)
//{
  ////frist check if exists
  //struct stat st;
  //if(stat(dirname.c_str(),&st) == 0){
    ////std::cout << "directory " << dirname << "already exists" << std::endl;
    ////std::cout << "trying to check if a directory named as " << dirname << "exists" << std::endl;
    ////wait for a second
    //RNG rng; gsl_rng* gslr; gslr = rng.get();
    //unsigned second = gsl_rng_uniform_int(gslr, 100000);
    //usleep(second*1000);
  //}
  ////
  //if(stat(dirname.c_str(),&st) == 0){;}
  ////then mkdir
  //else{
    //if (mkdir(dirname.c_str(),0755)==-1){
      //struct tm *current;
      //time_t now;
      //time(&now);
      //current = localtime(&now);

      //printf("At time %i:%i:%i: ", current->tm_hour, current->tm_min, current->tm_sec);
      //std::cout << "Unable to make directory named as " << dirname << std::endl << "Quit Now" << std::endl;
      //exit(-1);
    //}
    ////std::cout << "trying to make a directory named as " << dirname << std::endl;
  //}
  //// return final dirname
  //return dirname;
//}

//calculate sample minor allele frequence
vectorF SamMaf(vector2F& genos, std::string task)
{
  vectorI NonZero(genos[0].size(),0), NonMiss(genos[0].size(),0);
  for (UINT i=0; i<genos.size(); i++){
      for (UINT j=0; j < genos[i].size(); j++){
        if(genos[i][j] > 0) NonZero[j] += genos[i][j];
        if(genos[i][j] >= 0) NonMiss[j] ++;
      }
  }

  vectorF samMafs;
  for (UINT i = 0; i < NonZero.size(); i++){
    if(task == "2") samMafs.push_back((NonMiss[i] == 0)?0:float(NonZero[i])/float(NonMiss[i]));
    else samMafs.push_back((NonMiss[i] == 0)?0:float(NonZero[i])/float(2*NonMiss[i]));
  }
  return samMafs;
}

//for phased data only
void flipVector2F(vector2F& genos)
{
  for (UINT i=0; i<genos.size(); i++){
    for (UINT j=0; j < genos[i].size(); j++){
	genos[i][j] = mm_flip(genos[i][j]);
      }
    }
  return;
}

inline float mm_flip(float before)
{
  float after;
  if(before >= 0 && before <= 1) after = 1-before;
  else after = before;
  return after;
}
