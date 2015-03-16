#ifndef TDTTOOLS_H
#define TDTTOOLS_H
#include <sys/stat.h>
#include <sys/types.h>

//real data//
////////////
bool scan_unphasedTrioFile(std::string filename, vector2F& info);
bool scan_unphasedAnnFile(std::string filename, std::vector<std::string>& position, std::vector<std::string>& mutType, vectorF& samMaf, vectorF& popMaf, double deMafs);

bool scan_phasedTrioFile(std::string filename, vector2F& info, std::vector<std::string>& makers);
bool scan_phasedAnnFile(std::string filename, std::vector<std::string>& position, std::vector<std::string>& mutType, vectorF& samMaf, vectorF& popMaf, double deMafs);

bool scan_IndelTrioFile(std::string filename, vector2F& info, std::vector<std::string>& makers);
bool scan_IndelAnnFile(std::string filename, std::vector<std::string>& position, vectorF& samMaf, double deMafs);

// get the list of the files in folder dir, ignoring the files which has a '.' at the first char of file name
void getdir(std::string dir, std::vector<std::string> &files);
//make a directory
/*std::string makedir(std::string dirname);*/

vectorF SamMaf(vector2F& genos, std::string task);

//for phased data only
void flipVector2F(vector2F& genos);

inline float mm_flip(float before);

#endif
