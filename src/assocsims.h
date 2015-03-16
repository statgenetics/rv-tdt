// =====================================================================================
// 
//       Filename:  assocsims.h
// 
//    Description:  simulator for disease association studies
// 
//        Version:  1.0
//        Created:  01/04/2011 10:45:27 PM
//       Revision:  none
//       Compiler:  g++
// 
//         Author:  Gao Wang (gw), wangow@gmail.com
//                  Baylor College of Medicine, Texas, USA
//        License:  GNU General Public License <http://www.gnu.org/licenses/>
//                  Copyright (c) 2011, Gao Wang
// 
// =====================================================================================

#ifndef ASSOCSIMS_H
#define ASSOCSIMS_H

#include "gw_utilities.h"
//!\brief simulation of data for association tests

class gwSimulator {

  public:
    gwSimulator();

    gwSimulator(double boundary, double neutral_cutoff, const vectorF& pedInfos, const vectorF& mafs, const vector2F& genoFreqs, 
        const vectorF& fnctAnnotations, const vectorUI& positions);

    //!\fn createGenotypeComplexTraitsAssociations(...)
    //!\fn create_genotype_phenotype_mendelian(...)
    //!\fn mimic_genotyping(...)
    //!\fn create_ped_data_matrix(...)

    /*!\brief Generate complex traits genotype and phenotype data 
     * \param persons i.e. object to create: std::vector<gwPerson>& persons
     * \param pedInfos
     * \param mafs
     * \param genoFreqs
     * \param fnctAnnotations
     * \param positions
     * \param propFunctionalRv 
     * \param simulationTask
     * \param oddsRatios
     * \param baselinef
     * \param qtcoefs
     * \param qtcuts
     * \param shouldMarkCaseCtrl
     * \param pars 
     * \param isParConst 
     * \param moi
     * \param nPopulation
     * \param nCases
     * \param nCtrls
     * \param nUnphenotyped
     * \param gslr gsl_rng* object, initialize GNU/GSL rng
     * */
    ~gwSimulator();
    void createGenotypeComplexTraitsAssociations (const vectorF& propFunctionalRv, double baselinef, const vectorF& oddsRatios,  
        const char moi, unsigned nPopulation, gsl_rng* gslr, bool summary, std::string logFileName, const vector2UI& poolDat);

    void createGenotypeComplexTraitsAssociations (const vectorF& propFunctionalRv, double baselinef, const vectorF& oddsRatios,  
        const char moi, unsigned nCases, unsigned nCtrls, unsigned nUnphenotyped, gsl_rng* gslr, bool summary, std::string logFileName, const vector2UI& poolDat);

    void createGenotypeComplexTraitsAssociations (const vectorF& propFunctionalRv, const vectorF& pars, bool isParConst, const char moi,
        unsigned nCases, unsigned nCtrls, unsigned nUnphenotyped, gsl_rng* gslr);
    
    void createGenotypeComplexTraitsAssociations (const vectorF& propFunctionalRv, double baselinef, const vectorF& pars, bool isParConst, 
        const char moi, unsigned nCases, unsigned nCtrls, unsigned nUnphenotyped, gsl_rng* gslr, bool summary, std::string logFileName);

    void createGenotypeComplexTraitsAssociations (const vectorF& propFunctionalRv, const vectorF& qtcoefs, unsigned nPopulation, gsl_rng* gslr, const vector2UI& poolDat);

    void createGenotypeComplexTraitsAssociations (const vectorF& propFunctionalRv, const vectorF& qtcoefs, 
        const vectorF& qtcuts, bool shouldMarkCaseCtrl, unsigned nPopulation, unsigned nCases, unsigned nCtrls, 
        unsigned nUnphenotyped, gsl_rng* gslr, const vector2UI& poolDat);

    /*!\brief Generate Mendelian traits genotype and phenotype data*/
    void createGenotypeMendelianTraitsAssociations (double percentageCausal, bool isAllelicHeterogeneous,
        const char moi, unsigned nCases, double propHeterCases, unsigned nCtrls, gsl_rng* gslr);


    /*!\brief Add noise to genotypes and/or edit genotypes. Involves trimming genotypes, marking missing values, etc. 
     * \param proportionsMissingData	Portion of sites that has missing data, a vector of length 4: {deleterious%, protective%, non-causal%, synonymous%}
     * \param missingLowMaf data having maf below missingLowMaf will be marked missing (fail to genotype) 
     * \param shouldMarkMissing	Logical, if == true then "MARK_MISSING" (multiply the attribute by 1000) for future coding genotype as "MISSING_GENO(-9)" 
     * (otherwise will be coded "MAJOR_ALLELE(0)")
     * \param gslr gsl_rng* object, initialize GNU/GSL rng
     * */
    void mimicGenotyping ( const vectorF& proportionsMissingData, const double missingLowMaf,
        bool shouldMarkMissing, gsl_rng* gslr );


    /*!\brief Create the ped file as output and/or generate data matrix for assoc. analysis. 
     * \param isSynoTrimmed = true
     * \param isCvTrimmed = false
     * \param isPedWritten Whether or not to produce a ped file
     * \param projectName ped filename = projectName.ped
     * */
    void createPedfileMatrix(bool isSynoTrimmed, bool isCvTrimmed, 
        bool isPedWritten, std::string projectName, std::string simulationTask, vectorF& keepedMafs, vectorL& keepedSites); 

//    void calcVariantsPars(std::string filename, double mafLower, double mafUpper) const;
    
    vector2F getGenotypes() const;
    vectorF getPhenotypes() const;
    vectorF getObservedMafs() const;


  private:
    double __boundary;
    double __neutral_cutoff;
    //!\brief persons i.e. object to edit. std::vector<gwPerson>& persons
    std::vector<gwPerson> __persons;

    //!\brief 2D object, underlying genotype frequencies    
    vector2F __genoFreqs;
    //!\brief The first 6 columns of a *.ped file 
    vectorF __pedInfos;    
    //!\brief Additional phenotypes, if any
    vectorF __phenos;
    //!\brief Simulated current population MAF
    vectorF __mafs;
    //!\brief functional annotations (selection coef, synonymous nonsynonymous, etc) on each locus of the haplotype of simulated current generation
    vectorF __fnctAnnotations;    
    //!\brief positions in the gene region
    vectorUI __positions;
    //!\brief genotypes 2D vector of genotype matrix (to be modified)
    vector2F __mGenotypes;
    //!\brief phenotypes 1D vector of phenotypes (to be modified)
    vectorF __mPhenotypes; 
    //!\brief observedMafs 1D vector of observed MAF for each locus (to be modified)
    vectorF __mObservedMafs;
};
#endif///:~
