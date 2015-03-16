/*!\file unphased.h
   \brief

   Copyright: 2011 Zongxiao He @ Baylor College of Medicine
*/

#ifndef UNPHASEDTDT_H
#define UNPHASEDTDT_H

#include <algorithm>
#include <cmath>
#include "gw_utilities.h"

// only for MZ and VT test
class unphasedTDT
{
		public:
				unphasedTDT(const vector2F& genotypes, const vectorF& popMafs, const vectorF& samMafs, std::string& chrom, double siteMissRatio, bool SkipMiss);

				unphasedTDT(const vector3F& InputData, const vectorF& popMafs, const vectorF& samMafs, std::string& chrom, bool SkipMiss);

				~unphasedTDT();

				/*!\brief Get data information and analysis result 
				*/
				UINT getVariNum();
				UINT getFamNum();
				UINT getXvariNum();
				vectorL getFlip();
				vectorUI getErrFam();
				vectorF getTdtCountB();
				vectorF getTdtCountC();
				vectorF getSingleP();
				vectorF getSortedMafs();
				vectorF getSortedCQ();
				UINT getPermTimes();
				vectorF getMissRatio();
				vectorL getBeAnalyzed(); 
				vectorUI getDenovo();
				vectorF getWSSweight();

				//!\brief Morris & Zeggni test, return one-sided test result
				vectorF TdtMZ(double mafLower, double mafUpper, bool returnChi) const;

				//!\brief Variable threshold (VT)  test, return one-sided test
				vectorF TdtPermutWSS(double mafLower, double mafUpper, UINT adaptive, UINT PermutateTimes, double alpha, bool shuffleThree, std::string method);
				vectorF TdtPermutVT(double mafLower, double mafUpper, UINT adaptive, UINT PermutateTimes, double alpha, bool shuffleThree, std::string method);
				vectorF TdtPermutMZ(double mafLower, double mafUpper, UINT adaptive, UINT PermutateTimes, double alpha, bool shuffleThree, std::string method);
				vectorF TdtPermut(double mafLower, double mafUpper, UINT adaptive, UINT PermutateTimes, double alpha,  bool shuffleThree, std::string method);
				// single site test, return [[site, pvalue], [site, pvalue]...]
				vector2F unphasedSST(double sstMaf);

				// one sided VT test, with and without boundary
				vectorF m_tdtVT(double mafLower, double mafUpper);
		private:
				//
				vector2F __twoDGenos;
				//!\brief 3D object, trio genotypes
				vector3F __trioGenos;
				//!\brief Minor allele frequencies of listed variant sites, from trio sample
				vectorF __samMafs;
				//!\brief Minor allele frequencies of listed variant sites, from ESP population
				vectorF __popMafs;
				//which sites will not be analyzed, because missing
				vectorL __lowMissing;
				//whether it is a informative site
				vectorL __unInformSites;
				//!\brief 2D object, tdt B&C count for each family
				vector2F __TdtCountB;
				vector2F __TdtCountC;
				//!brief variant site number
				UINT __varNum;
				UINT __varNumToBeAnalyzed;
				//!brief whether skip trio has missing
				bool __SkipMiss;
				//!brief exclued variants number
				UINT __Xvariants;
				//!brief record denovo events
				vectorUI __denovo;
				//!\brief 2D object, non-transmitted chromosomes, calculated weight for WSS
				vector2F __NonTransChrom;
				//!biref The weight for WSS method
				vectorF __weight;

				//!\brief chromosome, to difference X and autosomal chromosomes
				std::string __chr;

				//!brief family number with wrong genetype data
				UINT __wrongFam;
				//!brief Need to flip the code?
				vectorL __flip;
				//!brief count how many families are un-analyzed due to missing or genotype error in dataset;
				vectorUI __discardFam;

				//!brief sorted popMafs, without replicates
				vectorF __sortedpopMafs;
				//!brief chi square based on sorted popMafs;
				vectorF __sortedChiSqu;
				//! brief permutation times
				UINT __permTimes;
				//! brief missing proportion
				vectorF __MissingRatio;

				void m_markMiss();
				/*!\brief Generate trio genotypes
				 * \param affectedOffsprings The genotypes of affected child
				 * \return 3D object, __trioGenos if passed the FamilyCheck
				 */
				vector3F m_geTrio3D(const vector2F& genos, UINT famSize);

				// flip the 0/1/2 code if samMafs > 0.5
				float m_flip(float before);
				void flipVector2F();

				// void m_maf(): tdt count on 3-D genotype data, call m_tdtCount
				void m_Count();

				// count B&C for one trio (both fat and mother) on one site, no missing
				// return father counts (element 0, =b,c,n) and mother counts (element 1) 
				vectorC m_tdtCount(double fat, double mot, double kid, bool denovo);

				/*!\brief Check the genotypes of each family
				*/
				bool m_FamilyCheck(const vector2F& OneFam) const;

				vectorF m_tdtWSS(double mafLower, double mafUpper);

				//adaptive permutation
				double m_check(UINT permcount, UINT iPermut, UINT adaptive, double alpha);

				/*!\brief basic function for genotype shuffle
				*/    
				vector3F m_GenoShuffle(bool shuffleThree) const;
				vector3F m_HapoShuffle() const;

				vectorF m_untransmitted(double fat, double mot, double kid);

				void m_SiteMissRatio();
				void m_TrimTdtTable();
				void TrimUninformSites();
				bool m_zeroFamily(double P, double M, double C);
};

//Calculate p values
inline double Pvalue(double& oriResult, vectorF& PermutedResult)
{
    //vectorF GreaterThanTest(2, 0);
    double GreaterThanTest(0);
    
    for (UINT Test = 0; Test < PermutedResult.size(); Test++)
    {
        // one-sided
        if (PermutedResult[Test] >= oriResult) {GreaterThanTest++;}
        // two-sided
        //if (abs(PermutedResult[Test][1]) >= abs(oriResult[1])) {GreaterThanTest[1]++;}
    }

    //for (UINT i = 0; i < GreaterThanTest.size(); i++)  GreaterThanTest[i] = float(GreaterThanTest[i])/PermutedResult.size();
    GreaterThanTest = float(GreaterThanTest)/PermutedResult.size();
    //GreaterThanTest.push_back(gsl_cdf_gaussian_Q(oriResult[0], 1));
    //GreaterThanTest.push_back(gsl_cdf_chisq_Q(oriResult[1], 1));
    
    return GreaterThanTest;
}

inline double IsIt(char a, char b)
{
    return (a == b)? 1:0;
}

#endif
