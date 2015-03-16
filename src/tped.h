#ifndef TPED_H
#define TPED_H

#include "gw_utilities.h"
#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>

typedef std::vector< std::string > vectorS;
typedef std::vector< std::vector< std::string > > vectorSS;

typedef std::map<std::string, vectorF> mapSvF;
typedef std::map<std::string, double> mapSF;
typedef std::map<std::string, std::string> mapSS;
typedef std::map<std::string, std::vector< std::string > > mapSvS;
typedef std::pair<std::string, vectorF> pairSvF;
typedef std::pair<std::string, double> pairSF;
typedef std::pair<std::string, std::string> pairSS;
typedef std::pair<std::string, std::vector< std::string > > pairSvS;

inline bool loadFile(std::string filename, vectorSS& output)
{
    std::ifstream file;
    file.open(filename.c_str());
    if (!file) {return false;}

    output.resize(0); //clear up
    std::string line;
    while(getline(file, line)){
        boost::trim(line); //get rid of whitespace on ends
        if(boost::starts_with(line, "#") || line.empty()) continue; //comments or blank

        std::vector<std::string> data;
        std::string value;
        std::istringstream iss(line);
        while(iss >> value) data.push_back(value);
        output.push_back(data);
    }
    return true;
}

struct ngsGene
{
    ngsGene() {;};
    ngsGene(std::string _name, mapSvF _genoP, mapSF _maf){
        geneName = _name;
        genoP = _genoP;
        maf = _maf;
    }
    std::string geneName;
    mapSvF genoP;
    mapSF maf;
};

class pedFileReader
{
    public:
        pedFileReader();
        pedFileReader(std::string genoFile, bool phased);
        void loadMapFile(std::string mapFile);
        void loadPhenoFile(std::string phenoFile);
        void loadMapPheno(std::string mapFile, std::string phenoFile);
        ~pedFileReader();

        void getGenotypes(vector3F& genos, vectorSS& snps, vectorS& geneNames, vector2F& mafs);
        void getPhenotype(vectorS& inds, mapSvS& phenos);
        void dumpOutGenos(std::string folder);

    private:
        mapSvF __rawGenos;
        bool __phased;
        std::vector<ngsGene> __genes;
        std::string __errorMessage;
        mapSvS __phenotypes;
        vectorS __indOrder;
        UINT __indNum;

        bool m_scanTped(std::string file, mapSvF& geno);
};

#endif
