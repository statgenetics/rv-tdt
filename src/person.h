// =====================================================================================
// 
//       Filename:  person.h
// 
//    Description:  "person" object implementation
// 
//        Version:  1.0
//        Created:  01/04/2011 10:54:58 PM
//       Revision:  none
//       Compiler:  g++
// 
//         Author:  Gao Wang (gw), wangow@gmail.com
//                  Baylor College of Medicine, Texas, USA
//        License:  GNU General Public License <http://www.gnu.org/licenses/>
//                  Copyright (c) 2011, Gao Wang
// 
// =====================================================================================


#ifndef PERSON_H
#define PERSON_H

#include <vector>
#include "gsl/gsl_rng.h"
#include "gw_utilities.h"

/*!\brief Simulate genotype and pheontype of a person under various models. Notes:
 * <br>
 * - 1) haplotype pool = current population simulated via simuPop etc. 
 * - 2) (TODO) current implementation does not have confounders yet ... 
 */

class gwPerson {

  public:
    /*!\brief Constructor, default
    */
    gwPerson();

    /*!\brief Constructor, verbose, w/o confounders
     * \param pedInfos The first 6 columns of a *.ped file 
     * \param mafs Simulated current population MAF
     * \param fnctAnnotations Funcational annotations on each locus of the haplotype of simulated current generation
     * \param snvPositions
     */
    gwPerson( double boundary, double neutral_cutoff, const vectorF& pedInfos, const vectorF& mafs, 
        const vectorF& fnctAnnotations, const vectorUI& snvPositions );

    /*!\brief Constructor, verbose, w/ confounders
     * \param pedInfos The first 6 columns of a *.ped file 
     * \param mafs Simulated current population MAF
     * \param fnctAnnotations Functional annotations on each locus of the haplotype of simulated current generation
     * \param snvPositions
     * \param confounders 2D object like {0, 1, 2}(values for confounder 1), {0, 1}(values for confounder 2), etc; can be either categorical or continuous
     * \param confounderProbs 2D object like {.3,.3, .4}(probabilities for values of confounder 1), {.8, .2}(probabilities for values of confounder 2), etc
     * \param confounderEffects 2d object like {1, 1.3, 2}(effects of confounder 1), {1, 1.5}(effects of confounder 2), etc; effects can be Odds ratio for binary traits or unit mean shift for QT 
     * \param areEpistatic A logical vector, whether a confounder acts as a modifier to a genetic variant
     */   
    gwPerson( double boundary, double neutral_cutoff, const vectorF& pedInfos, const vectorF& mafs, const vectorF& fnctAnnotations, 
        const vectorUI& snvPositions, const vector2F& confounders, const vector2F& confounderProbs, 
        const vector2F& confounderEffects, const vectorL& areEpistatic );

    /*!\brief Deconstrutor, nothing to do
    */
    ~gwPerson();

    /*!\brief Edit locus arributes. Can be edited via 1) (vector) causal or non-causal probabilities 2) (vector) missing data probabilities
     * \param proportionsDeleteriousProtective causal%, protective%
     * \param gslr gsl_rng* object, initialize GNU/GSL rng
     * \return 
     */
    void updateLocusAttributes(const vectorF& proportionsDeleteriousProtective, gsl_rng* gslr);
         

    /*!\brief Compute group genotype frequency by PAR model
     * \param parsInput Population attributible risk, a vector of length 2 (deleterious PAR and protective PAR)
     * \param isParConstant Logical, if true then marginal PARs for variants are constant, otherwise will be weighted by 1/MAF
     * \param moi Mode of inheritance
     * \return 
     */  
    void updateGenotypeFreqs(const vectorF& parsInput, bool isParConstant, const char moi);


    /*!\brief Calculate underlying locus odds ratios given group population attributable risk
     * \param parsInput Population attributible risk, a vector of length 2 (deleterious PAR and protective PAR)
     * \param isParConstant Logical, if true then marginal PARs for variants are constant, otherwise will be weighted by 1/MAF
     * \param moi Mode of inheritance
     * \return
     * */
    void updateLocusUnderlyingOddsRatios(const vectorF& parsInput, bool isParConstant, const char moi);


    /*!\brief Compute group genotype frequency by a given genotype frequencies
     * \param genotypeFreqs a given genotype frequencies vector2F object
     * \return 
     */  
    void updateGenotypeFreqs(const vector2F& genotypeFreqs);


    /*!\brief Generate one genotype (diplotype)
     * \param byGenoFreqs if == 0 then sample by population mafs, else sample by population (or group) genotype frequencies (1, 2 or 3 where 2 and 3 takes into consideration the phases). == 1 random phase, == 2 heterozygote on the first haplotype, == 3 heterozygote on the second haplotype 
     * \param gslr gsl_rng* object, initialize GNU/GSL rng
     * \return (updated genotype can be obtained via getGenotype() method)
     */     
    void generateGenotype(int byGenoFreqs, gsl_rng* gslr);

    void generateGenotype(const vector2UI& poolDat, gsl_rng* gslr);

    /*!\brief Compute genotypic effect to disease (Odds ratio in case-control traits)
     * \param oddsRatios Odds ratios arguements {deleterious OR lower limit, deleterious OR upper limit, protective OR lower limit, protective OR upper limit, OR for CV} 
     * OR lower limit == 0 means constant OR model, when ORc == upper limit) 
     * \param baselinef Baseline penetrance, (appox. = prevalence) 
     * \param moi mode of inheritance
     * \return Odds of disease given the genotype 
     */     
    double computeGenotypicEffect(const vectorF& oddsRatios, double baselinef, const char moi);
    double computeGenotypicEffect(double baselinef, const char moi) const;

    /*!\brief Compute genotypic effect to disease (mean-shift per variant in QT)
     * \param meanShifts QT mean shift (beta) per variant arguments {RV beta lower limit, RV beta upper limit, CV beta} Note: protective variants simply take an opposite sign
     * \return mean QT shift due to this genotype 
     */     
    double computeGenotypicEffect(const vectorF& meanShifts) const;


    /*!\brief get mafs indexes under the pth percentile of sample mafs
     * \param percentageCausal Percentage p of causal variants. Causal = the ones having MAF \in {X(1) ... X(floor(n*p))} X(i) is order statistic of MAFs
     * \return mafs indexes under the pth percentile of sample mafs
     */ 
    vectorUI getMendelianMafIdxes(double percentageCausal) const;


    /*!\brief Further edit locus attributes: by marking the ones with index NOT in mendelianMafIdxes as having attribute "neutral = 17 or 67". 
     * \param mendelianMafIdxes mafs indexes under the pth percentile of sample mafs
     * \param isAllelicHeterogeneous If is false then the one Not to mark is the one that ranks pth percentile in MAFs
     * \return 
     */
    void updateLocusAttributes(const vectorUI& mendelianMafIdxes, bool isAllelicHeterogeneous);


 
     /*!\brief Compute genotypic effect for Mendelian disease (Odds = 0 or Odds = DBL_MAX)
      *\param mendelianMafIdxes mafs indexes under the pth percentile of sample mafs 
      *\param isAllelicHeterogeneous If is false then the fixed causal variant is the one that ranks pth percentile
      *\param moi mode of inheritance
      *\return Odds of diease given the genotype (either 0 or DBL_MAX)
      */
   // double computeGenotypicEffect(const vectorUI& mendelianMafIdxes, bool isAllelicHeterogeneous, const char moi) const;


    /*!\Update genotype frequency based on if the person is a Mendelian Case
     * \param mendelianMafIdxes mafs indexes under the pth percentile of sample mafs 
     * \param isAllelicHeterogeneous If is false then the fixed causal variant is the one that ranks pth percentile
     * \param moi mode of inheritance
     * \param gslr gsl_rng* object, initialize GNU/GSL rng
     */
    void updateGenotypeFreqs(const vectorUI& mendelianMafIdxes, bool isAllelicHeterogeneous, const char moi, gsl_rng* gslr);
    

    /*!\Update genotype frequency based on if the person is a Mendelian Case
     * \param mendelianMafIdxes mafs indexes under the pth percentile of sample mafs 
     * \param isAllelicHeterogeneous If is false then the fixed causal variant is the one that ranks pth percentile
     * \param moi mode of inheritance
     */
    void updateGenotypeFreqs(const vectorUI& mendelianMafIdxes, bool isAllelicHeterogeneous, const char moi);


    /*!\brief Generate one phenotype (diplotype)
     * \param genotypicEffect Will be either Odds of Affected or mean shift of QT
     * \param isQt if == true then genotypicEffect is the mean shift in QT, otherwise is odds in case-control
     * \param gslr gsl_rng* object, initialize GNU/GSL rng
     * \return (updated genotype can be obtained via getGenotype() method)
     */ 
    void generatePhenotype(double genotypicEffect, bool isQt, gsl_rng* gslr);
   
      
    /*!\brief Update primary phenotype of the person 
     *\param phenotype pedInfos[5]
     * */
    void updatePhenotype(double phenotype);
     

    /*!\brief Get Phenotype of the person 
     *\return Phenotypes (more than one) of the person excluding pedInfos[5]
     * */
    vectorF getPhenotypes() const;
    
    
    /*!\brief Get Phenotype of the person 
     *\return Phenotype of the person, pedInfos[5]
     * */
    double getPhenotype() const;


    /*!\brief Further edit locus attributes: set missing data via 1) proportion of missing deleterious sites 2) proportion of missing protective sites 3) proportion of missing non-causal sites, 4) proportion of missing synonymous sites
     * \param proportionsMissingData Portion of sites that has missing data, a vector of length 4: {deleterious%, protective%, non-causal%, synonymous%} 
     * \param shouldMarkMissing Logical, if == true then "MARK_MISSING" (multiply the attribute by 1000) for future coding genotype as "MISSING_GENO(-9)" (otherwise will be coded "MAJOR_ALLELE(0)")
     * \param gslr gsl_rng* object, initialize GNU/GSL rng
     * \return 
     */
    void updateLocusAttributes(const vectorF& proportionsMissingData, bool shouldMarkMissing, gsl_rng* gslr);


    /*!\brief Further edit locus attributes: by taking in a pre-set attribute vector 
     * \param locusAttributes Same type as __locusAttributes; pre-set
     * \return 
     */
    void updateLocusAttributes(const vectorI& locusAttributes);


    /*!\brief Further edit locus attributes: by taking a logical vector 
     * \param shouldTrim if entry i is true then for the ith attribute "MARK_MISSING" (multiply the attribute by 9999) or "MARK_WILD" (multiply the attribute by 1000)
     * \return 
     */
    void updateLocusAttributes(const vectorL& shouldTrim, bool shouldMarkMissing);


    /*!\brief Update genotype by __locusAttributes. May recode person's genotype as "MISSING_ALLELE(-9)" or "MAJOR_ALLELE(0)"
     *\return
     * */
    void updateGenotype();
    

    /*!\brief Update genotype by external genotypes
     *\return
     * */
    void updateGenotype(const vector2F& genotype);


    /*!\brief Get genotype of the person 
     *\param isMissingTrimmed Logical, if == true trim genotype: by checking for and knocking site out
     *\param isSynoTrimmed Logical, if == true trim genotype: by checking for and knocking site out
     *\param isCvTrimmed Logical, if == true trim genotype: by checking for and knocking site out
     *\param untrimmedPositions SNV positions that are kept
     *\return genotype of the person
     * */
    vector2F getGenotype(bool isMissingTrimmed, bool isSynoTrimmed, bool isCvTrimmed, vectorUI& untrimmedPositions, vectorF& keepedMafs, vectorL& keepedSites) const;
    vector2F getGenotype() const;
    vector2F getGenotypeFreqs() const; 
    vectorF getLocusEffects() const;

    /*!\brief Print locus attributes summary. To translate the codings back into words 
     *\return 
     * */
    void summarizeLocusAttributes(std::string logname) const;


    /*!\brief Get "locusAttributes" of the person. Can be used to keep track of simulations (same as person.debug(4)), or as input of updateLocusAttributes(const vectorI& locusAttributes)   
     *\return locusAttributes of the person
     * */
    vectorI getLocusAttributes() const;
    vectorF getLocusMafs() const;


    /*!\brief Get positions of variants 
     *\return Positions of variants in a gene region
     * */
    vectorUI getSnvPositions() const;


    /*!\brief Verbose screen out-put, for debug mainly
    */
    void setVerbose(int verbose);
    void debug(int showWhat) const;


  private:
    double __boundary; 
    double __neutral_cutoff;
    //!\brief 2D object, genotypes of a person
    vector2F __genos;
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
    //!\brief Variant positions in the gene region
    vectorUI __snvPositions;
    //!\brief D_RV = 1, D_CV = 6, P_RV = -1, P_CV = -6, SYNO_RV = 15, SYNO_CV = 65, N_RV = 17, N_CV = 67 (\times MARK_MISSING = 9999 or MARK_WILD = 1000) 
    vectorI __locusAttributes;    
    //!\brief Effects per locus, if applicable 
    vectorF __locusEffects;
    //!\brief locus population attributable risk, if applicable
    vectorF __locusPars;

    //!\brief 2D object like {0, 1, 2}(values for confounder 1), {0, 1}(values for confounder 2), etc; can be either categorical or continuous
    vector2F __confounders;
    //!\brief 2D object like {.3,.3, .4}(probabilities for values of confounder 1), {.8, .2}(probabilities for values of confounder 2), etc
    vector2F __confounderProbs;
    //!\brief A logical vector, whether a confounder acts as a modifier to a genetic variant
    vectorL __areEpistatic;
    //!\brief 2d object like {1, 1.3, 2}(effects of confounder 1), {1, 1.5}(effects of confounder 2), etc; effects can be Odds ratio for binary traits or unit mean shift for QT 
    vector2F __confounderEffects;
    //!\brief Logical variable, flag for debug mode
    bool __isDebug;


    /*!\brief Calculate odds from penetrance
     * \param pen penetrance
     * \return odds for disease
     */ 
    inline double m_calcOdds(double pen) const 
    {
      return (pen / (1 - pen));
    }  


    /*!\brief Calculate odds from penetrance
     * \param odds 
     * \return penetrance
     */ 
    inline double m_calcPenetrance(double odds) const 
    {
      return (odds / (1 + odds));
    } 


    /*!\brief Calculate effect of a variant under variable effects model
     * \param pi 
     * \param effectMin
     * \param effectMax
     * \param mafMin
     * \param mafMax
     * \param direction
     * \return Effect for the locus
     */ 
    inline double m_calcVariantEffect(double pi, double effectMin, double effectMax, 
        double mafMin, double mafMax, double direction) const
    { 
      if (mafMax == mafMin) 
        return (effectMin + effectMax) / 2.0;
      if (effectMin >= direction) 
        return effectMin + (effectMax - effectMin) * (mafMax - pi) / (mafMax - mafMin);
      else 
        return effectMin + (effectMax - effectMin) * (pi - mafMin) / (mafMax - mafMin);
    }


    /*!\brief Calculate actual locus odds ratio based on its genotype, variant attribute, odds ratio model and mode of inheritance
     * \param index locus index
     * \param oddsRatios Odds ratios arguements {deleterious OR lower limit, deleterious OR upper limit, protective OR lower limit, protective OR upper limit, OR for CV}
     * \param baselinef Baseline penetrance, (appox. = prevalence) 
     * \param moi mode of inheritance
     * */
    double m_calcLocusOddsRatio(UINT index, const vectorF& oddsRatios, double baselinef, const char moi);
    // calculation based on __locusEffects
    double m_calcLocusJointOddsRatio(double baselinef, const char moi) const;

    inline void m_printPunches(int n) const
    {
      for (int i = 0; i != n; ++i)
        std::cout << "#";
      std::cout << "\n" << std::endl;
      return;
    }
};
#endif ///:~
