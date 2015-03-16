
#ifndef MAINREAL_H
#define MAINREAL_H

#include "gw_utilities.h"
#include "tped.h" 

inline bool fileExists(std::string filename)
{
    std::ifstream file;
    file.open(filename.c_str());
    if (!file) return false;
    else return true;
}

void readDataForRareTdt(vector3F& genos, vectorSS& snps, vectorS& geneNames, vector2F& popMafs, bool phased, std::string genoFile, std::string phenoFile, std::string mapFile);


/* vectorF recalculate_maf(vector2F& genos, vectorF& popmaf); */

#endif


