/*!\file unphased.h
   \brief

   Copyright: 2011 Zongxiao He @ Baylor College of Medicine
*/

#ifndef PHASEDTDT_H
#define PHASEDTDT_H

#include <algorithm>
#include <cmath>
#include "gw_utilities.h"

extern bool DEBUG;
extern bool VTDEBUG;

// only for MZ and VT test
class phasedTDT
{
    public:

        /*!\brief Constructor
         * \param affectedChildren The genotypes of affected offsprings, 0/1/-9 coding
         */
        //phasedTDT(const vector2F& affectedChildren, const vectorF& popMafs, const vectorF& samMafs, std::string& chrom, double siteMissRatio, bool SkipMiss);
	phasedTDT(const vector2F& affectedChildren, const vectorF& popMafs, const vectorF& samMafs, std::string& chrom, double siteMissRatio, bool SkipMiss);

        /*!\brief Another constructor, used for permutation
         * \param InputData The genotypes of trios, result of GenoShuffle or HypoShuffle
         */
        //phasedTDT(const vector3F& InputData, const vectorF& popMafs, const vectorF& samMafs, std::string& chrom, bool SkipMiss);
	phasedTDT(const vector3F& InputData, const vectorF& popMafs, const vectorF& samMafs, std::string& chrom, bool SkipMiss);

        /*!\brief Deconstructor
         */
        ~phasedTDT();

        /*!\brief Get data information and analysis result
         */
	UINT getVariNum();
	UINT getFamNum();
	UINT getXvariNum();
	vectorL getFlip();
        //vectorUI getErrFam();
	vectorF getTdtCountB();
	vectorF getTdtCountC();
	vectorF getSingleP();
	vectorF getSortedMafs();
	vectorF getSortedCQ();
	UINT getPermTimes();
	vectorF getMissRatio();
	vectorF getpopMafs();
	vectorF getsamMafs();
	vectorUI getDenovo();
	vectorF getWSSweight();
	vectorL getBeAnalyzed();

        //!\brief Morris & Zeggni test, return one-sided test result, with and without boundary
        vectorF TdtMZ(double mafLower, double mafUpper) const;
	//!brief return one-sided test result, with and without boundary
	vectorF TdtCMC(double mafLower, double mafUpper) const;
        //!brief return one-sided test result, with and without boundary * 2
        vectorF TdtMZCMC(double mafLower, double mafUpper) const;
	//!brief return one-sided test result, with boundary, for WSS and VT
        vectorF TdtPermutWSS(double mafLower, double mafUpper, UINT adaptive, UINT PermutateTimes, double alpha, std::string method);
	vectorF TdtPermutMZ(double mafLower, double mafUpper, UINT adaptive, UINT PermutateTimes, double alpha, std::string method);
	vectorF TdtPermutVT(double mafLower, double mafUpper, UINT adaptive, UINT PermutateTimes, double alpha, std::string VTaggregate, std::string method);
	vectorF TdtPermut(double mafLower, double mafUpper, UINT adaptive, UINT PermutateTimes, double alpha, std::string method);
	vector2F phasedSST(double sstMaf);

    private:
	//
	vector2F __twoDGenos;
        //!\brief 3D object, trio genotypes
        vector3F __trioGenos;
        //!\brief Minor allele frequencies of listed variant sites, from trio sample
        vectorF __samMafs;
        //!\brief Minor allele frequencies of listed variant sites, from ESP population
        vectorF __popMafs;
        //!\brief chromosome, to difference X and autosomal chromosomes
        std::string __chr;
        //!brief family number with wrong genetype data
        UINT __wrongFam;
	//!brief whether skip trio has missing
	bool __SkipMiss;
	//!brief exclued variants number
	UINT __Xvariants;
        //!brief variant site number
        UINT __varNum;
	UINT __varNumToBeAnalyzed;
	//which sites will not be analyzed, because missing
	vectorL __lowMissing;
	//whether it is a informative site
	vectorL __unInformSites;
	//! brief missing proportion
	vectorF __MissingRatio;
	vectorL __flip;

        //!brief count how many families are un-analyzed due to missing or genotype error in dataset;
        //vectorUI __discardFam;
	//!brief record denovo events
	vectorUI __denovo;
        //!\brief 2D object, tdt B&C count for each family
        vector2F __TdtCountB;
        vector2F __TdtCountC;
	//!brief Chromosome Origin information
	vector2UI __chromOri;
        //!\brief 2D object, non-transmitted chromosomes, calculated weight for WSS
        vector2F __NonTransChrom;
	//!biref The weight for WSS method
	vectorF __weight;

	//!brief WSS and VT tests
	vectorF m_tdtWSS(double mafLower, double mafUpper);
	vectorF m_tdtVT(double mafLower, double mafUpper, std::string aggregate);

	//!brief sorted popMafs, without replicates
	vectorF __sortedpopMafs;
	//!brief chi square based on sorted popMafs;
	vectorF __sortedChiSqu;
	//! brief permutation times that run in GenoPermut function
	UINT __permTimes;
	
	float m_flip(float before);
	void m_markMiss();
        /*!\brief Generate trio genotypes
         * \param affectedOffsprings The genotypes of affected child
         * \return 3D object, __trioGenos if passed the FamilyCheck
         */
        vector3F m_geTrio3D(const vector2F& genos, UINT famSize);

        // void m_maf(): tdt count on 3-D genotype data, call m_tdtCount
        void m_TdtCountForAllFamilies();

        // count B&C for one trio (both fat and mother) on one site, no missing
        // return father counts (element 0, =b,c,n) and mother counts (element 1)
        vectorC m_tdtCount(double fat, double mot, double kid, bool denovo);

        /*!\brief Check the genotypes of each family
         */
        bool m_FamilyCheck(const vector2F& OneFam) const;

        void m_TdtCalForOneParent(vector2F& Parent, vectorF& Child, int whichChrom, vectorL& shouldBeAnalyzed);

        //adaptive permutation
        double m_check(UINT permcount, UINT iPermut, UINT adaptive, double alpha);

        /*!\brief basic function for genotype shuffle
         */
        vector3F m_GenoShuffle() const;
	vector3F m_HapoShuffle() const;
	
	//void m_SiteMissRatio();
	//vector2F m_TrimSites(vector2F& before, vectorL& should);
	bool m_zeroFamily(double P1, double P2, double M1, double M2, double C1, double C2);
	
	void m_TrimTdtTable();
	void TrimUninformSites();
	void m_SiteMissRatio();
};

bool m_ChromOrigin(const vector2F& Family, vectorUI& kidAlleleOri, bool Missing, bool Denovo);

bool m_ChromEqual(const vectorF& ChromOne, const vectorF& ChromTwo, bool Missing, bool Denovo);

inline bool m_IsItInform(const vectorF& parent)
{
    for (unsigned int locus = 0; locus < parent.size(); locus++)
    {
        if(parent[locus] != 0) return true;
    }
    return false;
}

//Calculate weight for WSS
inline double m_WssWeightCal(int n, int N, int M)
{
    double q = float(M+1)/(2*N+3);
    return (n==0)?0:1/sqrt(n*q*(1-q));
}


//inline bool m_zeroFamily(double Pone, double Ptwo, double Mone, double Mtwo, double Cone, double Ctwo)
//{
//    if((Pone == -9.0 || Ptwo == -9.0) && (Mone == -9.0 || Mtwo == -9.0)) return true;
//    else if(Cone == -9.0 || Ctwo == -9.0) return true;
//    else if(Pone == Ptwo && Pone == Mone && Pone == Mtwo && Pone == Cone && Pone == Ctwo) return true;
//    else return false;
//}
#endif
