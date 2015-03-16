#ifndef RARETDT_H
#define RARETDT_H

#include "gw_utilities.h"
#include "tdtTrio.h"
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/shared_ptr.hpp>
#include "fisher2.h"

typedef boost::shared_ptr<trio> trioPtr;
typedef boost::shared_ptr<phasedTrio> phasedTrioPtr;
typedef boost::shared_ptr<unPhasedTrio> unPhasedTrioPtr;
typedef std::vector<trioPtr>::iterator trioIterator;


struct tdtResult
{
    tdtResult() {}
    tdtResult(double C_MZ, double C_CMC, double P_MZ, double P_CMC){
        mz_chi2=C_MZ;
        cmc_chi2=C_CMC;
        mz_pval=P_MZ;
        cmc_pval=P_CMC;
    }
	void print() {
		printf(">> %s=%+-4.2f %s=%+-4.2f %s=%+-4.2f %s=%+-4.2f\n", "mz_chi2", mz_chi2,"cmc_chi2", cmc_chi2,"mz_pval", mz_pval,"cmc_pval", cmc_pval);
	}

    double mz_chi2, cmc_chi2, mz_pval, cmc_pval;
};

struct rvTDTshuffleOutput
{
		rvTDTshuffleOutput() {}
		rvTDTshuffleOutput(const vector2F& __unTrans, const vector2F& __tdtBtable, const vector2F& __tdtCtable){
				unTrans=__unTrans;
				tdtBtable=__tdtBtable;
				tdtCtable=__tdtCtable;
		}
		vector2F unTrans, tdtBtable, tdtCtable;
};




class rareTdt
{
		public:
				rareTdt(bool phased, vector2F& genotypes, const vectorF& popMafs, bool skipMiss, double missCutoff, UINT minVarNum);
				rareTdt(bool phased, vector2F& genotypes, bool skipMiss, double missCutoff, UINT minVarNum);
				~rareTdt(){}
				void tdtTest(double lowerMaf, double upperMaf, UINT adaptive, UINT PermutateTimes, double alpha, std::map<std::string, double>& pvalues, std::vector<std::string>& tests, std::string mode);
				Json::Value getJsonLog() const {return logRoot;}
				bool checkInformGene() {return __informGene;}

				void noPermut(double lowerMaf, double upperMaf, std::map<std::string, double>& pvalues, std::vector<std::string>& tests);
				void updateWeights(const vectorF& weights);

				void ParentalBiasTest(double lowerMaf, double upperMaf, std::map<std::string, double>& pvalues, std::vector<std::string>& tests);

		private:
				Json::Value logRoot;
				std::vector<trioPtr> __trios;
				bool __phased;
				vector2F __genotypes;
				vectorF __popMafs;
				bool __hasPopMaf;
				bool __skipMiss;
				double __missCutoff;
				bool __informGene;
				vectorF __samMafs;
				UINT __nInformTrios;
				UINT __rowsPerTrio;
				vectorF __missingRatio;
				vectorL __varToBeAnalyzed;
				vectorUI __denovoCount;
				vectorUI __effTriosBySite;
				vector2F __tdtBtable;
				vector2F __tdtCtable;
				vectorF __wssWeight;
				bool __external_weights;

				std::string __test_mode;

				void m_preprocess();
				UINT m_getInformVarNum() const;
				void m_load();
				void m_tdtTableLoad();
				vectorF m_wssWeight(const vector2F& unTrans);
				vectorF m_singleSiteP() const;

				bool __first_call;

				rvTDTshuffleOutput m_shuffle(std::string shuffleMethod, gsl_rng* gslr);
				tdtResult m_rvAggreate(const vector2F& tdtBtable, const vector2F& tdtCtable, const vectorL& tobeAnalyzed, const vectorF& Weight, bool allowDenovo) ;
				tdtResult m_VT(const vector2F& tdtBtable, const vector2F& tdtCtable, const vectorL& tobeAnalyzed, const vectorF& weight) ;

				vectorF m_noPermut();
				vectorF m_Permut(UINT adaptive, UINT PermutateTimes, double alpha, gsl_rng* gslr);
				vectorF m_MzCmcPermut(vectorL& tobeAnalyzed, UINT adaptive, UINT PermutateTimes, double alpha, gsl_rng* gslr, double upperMaf);
				vectorF m_CmpPermut(std::string GenoHapo, const vectorL& tobeAnalyzed, UINT adaptive, UINT PermutateTimes, double alpha, gsl_rng* gslr, double lowerMaf, double upperMaf);

				vectorF  m_parentalTransmissionCount(const vector2F& tdtBtable, const vector2F& tdtCtable, const vectorL& tobeAnalyzed, bool allowDenovo) const;
};

inline vectorF m_colTotal(const vector2F& table)
{
		vectorF count(table[0].size(),0);
		for(UINT i = 0; i < table.size(); i++){
				for(UINT j = 0; j < table[i].size(); j++){
						count[j] += (table[i][j]>=100)?(table[i][j]-100):table[i][j];
				}
		}
		return count;
}

inline void m_updateCount(vectorUI& permCount, vectorF& oriRes, vectorF& permRes, gsl_rng* gslr)
{
		if(permCount.size() != oriRes.size() || permCount.size() != permRes.size()) {return;}
		for(UINT i=0; i<permCount.size(); i++){
				double diff=permRes[i] - oriRes[i];
				if(diff>1e-6) permCount[i]++;
				else if(diff<=1e-6 && diff>=-1e-6) {if(gsl_rng_uniform(gslr) > 0.5) permCount[i]++;}
				else {;}
		}
		return;
}

inline void m_updatePvalue(vectorF& permPval, vectorUI& permCount, UINT iPermut, double alpha)
{
		for(UINT i=0; i<permCount.size(); i++){
				double pval = (permCount[i]+1)*1.0/((iPermut+1)*1.0);
				double sigma = sqrt(pval*(1-pval)/iPermut);
				double beta = 0.05;
				double gs = gsl_cdf_gaussian_Pinv(1.0-beta/2.0, sigma);

				permPval[i] = (pval-gs > alpha)? pval:9.0;
		}
		return;
}

inline void m_ifTrueThenAddOne(const vectorL& logic, vectorUI& count)
{
		if(count.size()==0) count.resize(logic.size(), 0);
		if(count.size() != logic.size()) return;
		for(UINT i=0; i<logic.size(); i++){
				if(logic[i]) count[i]++;
		}
		return;
}

template<class T> std::vector<T> eraseTrue (std::vector<T>& vec, const vectorL stdV)
{
		std::vector<T> output;
		for(UINT i=0; i<vec.size(); i++){
				if(!stdV[i]) output.push_back(vec[i]); // remove element that is true in stdV
		}
		return output;
}

// borrowed from http://pngu.mgh.harvard.edu/~purcell/plink/dist/scratch/plink-1.07-src/fisher.cpp
// 2x2 only
inline double myfisher(vectorUI res)
{

	int nrow = 2;
	int ncol = 2;
	double table[4];

	int c=0;
    for (unsigned int i=0; i<res.size(); i++)    
	{
		table[c++] = res[i];
	}

	double expect = -1.0;
	double percnt = 100.0;
	double emin = 0;
	double pre = 0, prt = 0;
	int ws = 300000;
  
	fexact(&nrow, &ncol, table, &nrow, &expect, &percnt, &emin, &prt, &pre, &ws);

	return pre;
  
}

#endif
