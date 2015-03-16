// =====================================================================================
// 
//       Filename:  assocsims.cpp
// 
//    Description:  simulator for disease association studies
// 
//        Version:  1.0
//        Created:  01/04/2011 10:40:38 PM
//       Revision:  none
//       Compiler:  g++
// 
//         Author:  Gao Wang (gw), wangow@gmail.com
//                  Baylor College of Medicine, Texas, USA
//        License:  GNU General Public License <http://www.gnu.org/licenses/>
//                  Copyright (c) 2011, Gao Wang
// 
// =====================================================================================

#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <cstdlib>
#include <algorithm>
#include "gsl/gsl_cdf.h"
#include "gsl/gsl_randist.h"

#include "gw_utilities.h"
#include "person.h"
#include "assocsims.h"

namespace {
  const double AFFECTED = 2.0, UNAFFECTED = 1.0, UNPHENOTYPED = 0.0, MISSING_ALLELE = -9.0, MAJOR_ALLELE = 0.0; 
}
  
//!\brief This program is the core simulator function that simulates complex traits via combined OR model, PAR model, QT model or Extreme QT models.
//use the "gwPerson" class

gwSimulator::gwSimulator(double boundary, double neutral_cutoff, const vectorF& pedInfos, const vectorF& mafs, const vector2F& genoFreqs, 
        const vectorF& fnctAnnotations, const vectorUI& positions)
{
  __boundary = boundary;
  __neutral_cutoff = neutral_cutoff;
  __pedInfos = pedInfos;
  __mafs = mafs;
  __genoFreqs = genoFreqs;
  __fnctAnnotations = fnctAnnotations;
  __positions = positions;
  __persons.resize(0);
  __mObservedMafs.resize(0);
}
gwSimulator::gwSimulator() {}
gwSimulator::~gwSimulator() {}

//!-# "pop-odds"
void gwSimulator::createGenotypeComplexTraitsAssociations (const vectorF& propFunctionalRv, double baselinef, 
    const vectorF& oddsRatios, const char moi, unsigned nPopulation, gsl_rng* gslr, bool summary, std::string logFileName, const vector2UI& poolDat) 
{
  //!- Initialize a "template" person object
  gwPerson* personTemplate = new gwPerson (__boundary, __neutral_cutoff, __pedInfos, __mafs, __fnctAnnotations, __positions);
  personTemplate->updateLocusAttributes(propFunctionalRv, gslr);

  bool isInputOk = (nPopulation > 0);
  if (isInputOk);
  else {
    std::cerr << "Population size not valid. Nothing to do. Quit now." << std::endl;
    exit(-1);
  }

  //!- actual prevalence, will be computed based on updated penetrance of each person in the population
  double prevalence = 0.0; 

  //!- Generate population with disease status, accept all generated samples
  for (unsigned i = 0; i != nPopulation; ++i) {
    if (poolDat.size() == 0) personTemplate->generateGenotype(0, gslr);
    else personTemplate->generateGenotype(poolDat, gslr);
    double odds = personTemplate->computeGenotypicEffect(oddsRatios, baselinef, moi);
    personTemplate->generatePhenotype(odds, 0, gslr);
    __persons.push_back(*personTemplate);
    prevalence += odds / (1 + odds);
  }

  if (summary) 
  {
    //!- Output the disease prevalence for the binary trait
    std::string outname = logFileName + ".prevalence";
    std::ofstream myout;
    myout.open(outname.c_str(), std::ios::app);
    myout << "prevalence: " << prevalence / (1.0 * nPopulation) << std::endl;
    myout.close();
  }
  //  std::cout << personTemplate->getPhenotype() << std::endl;
  delete personTemplate; 
  return;
}


//!-# "dichot-odds"
void gwSimulator::createGenotypeComplexTraitsAssociations (const vectorF& propFunctionalRv, double baselinef, 
    const vectorF& oddsRatios, const char moi, unsigned nCases, unsigned nCtrls, unsigned nUnphenotyped, gsl_rng* gslr, bool summary, 
    std::string logFileName, const vector2UI& poolDat)
{
  //!- Initialize a "template" person object
  gwPerson* personTemplate = new gwPerson (__boundary, __neutral_cutoff, __pedInfos, __mafs, __fnctAnnotations, __positions);
  personTemplate->updateLocusAttributes(propFunctionalRv, gslr);

  bool isInputOk = (nCases > 0 || nCtrls > 0 || nUnphenotyped > 0);
  if (isInputOk);
  else {
    std::cerr << "Case/Control/Cohort size not valid. Nothing to do. Quit now." << std::endl;
    exit(-1);
  }

  //!- Generate case-ctrl samples (resample, reject or accept)
  unsigned iCase = 0, iCtrl = 0, iUnphenotyped = 0;

  while (iCase != nCases || iCtrl != nCtrls) {
    if (poolDat.size() == 0) personTemplate->generateGenotype(0, gslr);
    else personTemplate->generateGenotype(poolDat, gslr);

    double odds = personTemplate->computeGenotypicEffect(oddsRatios, baselinef, moi);
    personTemplate->generatePhenotype(odds, 0, gslr);
    double trait = personTemplate->getPhenotype();

    if (trait == AFFECTED) {
      // get a case
      if (iCase != nCases) {  
        // if case is still in need, collect it.
        __persons.push_back(*personTemplate);
        ++iCase;
      }
      else;  
      // if case is enough, do nothing.
    }
    else {
      // get a control
      if (iCtrl != nCtrls) {  
        // if control is still in need, collect it.
        __persons.push_back(*personTemplate);
        ++iCtrl;
      }
      else;  
      // if control is enough, do nothing.          
    }
  }

  personTemplate->updatePhenotype(UNPHENOTYPED); 

  //!- Generate unphenotyped cohorts
  while (iUnphenotyped != nUnphenotyped) {
    if (poolDat.size() == 0) personTemplate->generateGenotype(0, gslr);
    else personTemplate->generateGenotype(poolDat, gslr);
    __persons.push_back(*personTemplate);
    ++iUnphenotyped; 
  }
  //  std::cout << personTemplate->getPhenotype() << std::endl;
  if (summary) 
    personTemplate->summarizeLocusAttributes(logFileName);
  delete personTemplate;
  return;
}

 
//!-# "dichot-par"
void gwSimulator::createGenotypeComplexTraitsAssociations (const vectorF& propFunctionalRv, const vectorF& pars, bool isParConst, const char moi,
        unsigned nCases, unsigned nCtrls, unsigned nUnphenotyped, gsl_rng* gslr)
{  
  //!- Initialize a "template" person object
  gwPerson* personTemplate = new gwPerson (__boundary, __neutral_cutoff, __pedInfos, __mafs, __fnctAnnotations, __positions);
  personTemplate->updateLocusAttributes(propFunctionalRv, gslr);
  
  bool isInputOk = (nCases > 0 || nCtrls > 0 || nUnphenotyped > 0);
  if (isInputOk);
  else {
    std::cerr << "Case/Control/Cohort size not valid. Nothing to do. Quit now." << std::endl;
    exit(-1);
  }

  //!- Generate case-ctrl samples. No rejecting because case/ctrl use different underlying MAF
  vector2F genoFreqsRecovery = __genoFreqs;
  personTemplate->updatePhenotype(AFFECTED); 
  for (unsigned i = 0; i != nCases; ++i) {
    personTemplate->updateGenotypeFreqs(pars, isParConst, moi); 
    personTemplate->generateGenotype(1, gslr);
    __persons.push_back(*personTemplate);
    personTemplate->updateGenotypeFreqs(genoFreqsRecovery);
  }

  personTemplate->updatePhenotype(UNAFFECTED); 
  for (unsigned i = 0; i != nCtrls; ++i) {
    personTemplate->updateGenotypeFreqs(pars, isParConst, moi); 
    personTemplate->generateGenotype(1, gslr);
    __persons.push_back(*personTemplate);
    personTemplate->updateGenotypeFreqs(genoFreqsRecovery);
  }

  //!- Generate unphenotyped cohorts. Use the MAF directly from SRV_batch

  personTemplate->updatePhenotype(UNPHENOTYPED); 
  for (unsigned i = 0; i != nUnphenotyped; ++i) {
    personTemplate->generateGenotype(0, gslr);
    __persons.push_back(*personTemplate);
  }
  //  std::cout << personTemplate->getPhenotype() << std::endl;
  delete personTemplate;
  return;
}


//!-# "dichot-par-odds"
void gwSimulator::createGenotypeComplexTraitsAssociations (const vectorF& propFunctionalRv, double baselinef, 
    const vectorF& pars, bool isParConst, const char moi, unsigned nCases, unsigned nCtrls, unsigned nUnphenotyped, gsl_rng* gslr, bool summary, std::string logFileName)
{
  //!- Initialize a "template" person object
  gwPerson* personTemplate = new gwPerson (__boundary, __neutral_cutoff, __pedInfos, __mafs, __fnctAnnotations, __positions);
  personTemplate->updateLocusAttributes(propFunctionalRv, gslr);

  bool isInputOk = (nCases > 0 || nCtrls > 0 || nUnphenotyped > 0);
  if (isInputOk);
  else {
    std::cerr << "Case/Control/Cohort size not valid. Nothing to do. Quit now." << std::endl;
    exit(-1);
  }
  
  personTemplate->updateLocusUnderlyingOddsRatios(pars, isParConst, moi);

  //!- Generate case-ctrl samples (resample, reject or accept)
  unsigned iCase = 0, iCtrl = 0, iUnphenotyped = 0;

  while (iCase != nCases || iCtrl != nCtrls) {
    personTemplate->generateGenotype(0, gslr);
    double odds = personTemplate->computeGenotypicEffect(baselinef, moi);
    personTemplate->generatePhenotype(odds, 0, gslr);
    double trait = personTemplate->getPhenotype();

    if (trait == AFFECTED) {
      // get a case
      if (iCase != nCases) {  
        // if case is still in need, collect it.
        __persons.push_back(*personTemplate);
        ++iCase;
      }
      else;  
      // if case is enough, do nothing.
    }
    else {
      // get a control
      if (iCtrl != nCtrls) {  
        // if control is still in need, collect it.
        __persons.push_back(*personTemplate);
        ++iCtrl;
      }
      else;  
      // if control is enough, do nothing.          
    }
  }

  personTemplate->updatePhenotype(UNPHENOTYPED); 

  //!- Generate unphenotyped cohorts
  while (iUnphenotyped != nUnphenotyped) {
    personTemplate->generateGenotype(0, gslr);
    __persons.push_back(*personTemplate);
    ++iUnphenotyped; 
  }
  //  std::cout << personTemplate->getPhenotype() << std::endl;
  if (summary) 
    personTemplate->summarizeLocusAttributes(logFileName);
  
  delete personTemplate;
  return;
}


//!-# "qt"
void gwSimulator::createGenotypeComplexTraitsAssociations (const vectorF& propFunctionalRv, 
    const vectorF& qtcoefs, unsigned nPopulation, gsl_rng* gslr, const vector2UI& poolDat)
{
  //!- Initialize a "template" person object
  gwPerson* personTemplate = new gwPerson (__boundary, __neutral_cutoff, __pedInfos, __mafs, __fnctAnnotations, __positions);
  personTemplate->updateLocusAttributes(propFunctionalRv, gslr);
  bool isInputOk = (nPopulation > 0);
  if (isInputOk);
  else {
    std::cerr << "Population size not valid. Nothing to do. Quit now." << std::endl;
    exit(-1);
  }

  //!- Generate population with QT values, accept all generated samples
  for (unsigned i = 0; i != nPopulation; ++i) {
    if (poolDat.size() == 0) personTemplate->generateGenotype(0, gslr);
    else personTemplate->generateGenotype(poolDat, gslr);
    double meanShift = personTemplate->computeGenotypicEffect(qtcoefs);
    personTemplate->generatePhenotype(meanShift, 1, gslr);
    __persons.push_back(*personTemplate);
  }
  //  std::cout << personTemplate->getPhenotype() << std::endl;
  delete personTemplate;
  return;
}


//!-# "dichot-qt" (finite and infinite)
void gwSimulator::createGenotypeComplexTraitsAssociations (const vectorF& propFunctionalRv, const vectorF& qtcoefs, 
        const vectorF& qtcuts, bool shouldMarkCaseCtrl, unsigned nPopulation, unsigned nCases, unsigned nCtrls, 
        unsigned nUnphenotyped, gsl_rng* gslr, const vector2UI& poolDat) 
{
  //!- Initialize a "template" person object
  gwPerson* personTemplate = new gwPerson (__boundary, __neutral_cutoff, __pedInfos, __mafs, __fnctAnnotations, __positions);
  personTemplate->updateLocusAttributes(propFunctionalRv, gslr);

  //!- most tidious sampling scheme to program ... 
  //!- #0 = sample from infinite population, #!0 = sample extremes from finite cohort
  bool isInputOk = (qtcuts.size() == 2 && qtcuts[0] > 0.0 
      && qtcuts[1] < 1.0 && qtcuts[1] >= qtcuts[0]);
  if (isInputOk);
  else {
    std::cerr << "Extreme Qt lower/upper percentiles not valid. Quit now." << std::endl;
    exit(-1);
  }

  isInputOk = (((nCases > 0 || nCtrls > 0) && nPopulation == 0) || 
      (nPopulation > 0 && nCases == 0 && nCtrls == 0));
  if (isInputOk);
  else {
    std::cerr << "Case/Control/Cohort or population size not valid. Nothing to do. Quit now." << std::endl;
    exit(-1);
  }

  switch (nPopulation) {
    //!- extreme QT from unlimited population, using "qtcuts[0]" as percentile cut-off for ctrls, "qtcuts[1]" for cases
    //!- Will/Wont mark the phoentype as binary, determined by "shouldMarkCaseCtrl"
    case 0 :
      {
        unsigned iCase = 0, iCtrl = 0;

        double lower = gsl_cdf_ugaussian_Pinv(qtcuts[0]);
        double upper = gsl_cdf_ugaussian_Pinv(qtcuts[1]);

        while (iCase != nCases || iCtrl != nCtrls) {
          if (poolDat.size() == 0) personTemplate->generateGenotype(0, gslr);
          else personTemplate->generateGenotype(poolDat, gslr);
          double meanShift = personTemplate->computeGenotypicEffect(qtcoefs);
          personTemplate->generatePhenotype(meanShift, 1, gslr);
          double trait = personTemplate->getPhenotype();

          if (trait >= upper) {
            // get a case
            if (iCase != nCases) {  
              // if case is still in need, collect it.
              if (shouldMarkCaseCtrl) 
                personTemplate->updatePhenotype(AFFECTED);

              __persons.push_back(*personTemplate);
              ++iCase;
            }
            else;  
            // if case is enough, do nothing.
          }

          else if (trait <= lower) {
            // get a control
            if (iCtrl != nCtrls) {  
              // if control is still in need, collect it.
              if (shouldMarkCaseCtrl) 
                personTemplate->updatePhenotype(UNAFFECTED);

              __persons.push_back(*personTemplate);
              ++iCtrl;
            }
            else;  
            // if control is enough, do nothing.          
          }

          else;
        }
      }
      break;

    default :
      {
        //!- extreme QT from given population, using "qtcuts[0]" as percentile cut-off for ctrls, "qtcuts[1]" for cases
        //!- Needs to convert bounds into index (floored by type conversion -- no need to worry about edge problem)
        unsigned lower = (unsigned) (nPopulation * qtcuts[0]);
        unsigned upper = (unsigned) (nPopulation * qtcuts[1]);
        if (upper == lower) ++upper;

        vectorF sortedPhenotypes(nPopulation);

        //!- First generate a cohort with QT known, as in "qt" model
        for (unsigned i = 0; i != nPopulation; ++i) {
          if (poolDat.size() == 0) personTemplate->generateGenotype(0, gslr);
          else personTemplate->generateGenotype(poolDat, gslr);
          double meanShift = personTemplate->computeGenotypicEffect(qtcoefs);
          personTemplate->generatePhenotype(meanShift, 1, gslr);
          __persons.push_back(*personTemplate);
          sortedPhenotypes[i] = personTemplate->getPhenotype();
        }

        //!- Need to sort __mPhenotypes to figure out extreme values
        sort (sortedPhenotypes.begin(), sortedPhenotypes.end());

        //!- Collect only extreme samples
        //!- Will/Wont mark the phoentype as binary, determined by "shouldMarkCaseCtrl"
        std::vector<gwPerson>* __personsExtreme = new std::vector<gwPerson>();
        for (unsigned i = 0; i != __persons.size(); ++i) {
          if (__persons[i].getPhenotype() >= sortedPhenotypes[upper]) {
            if (shouldMarkCaseCtrl) { 
              __persons[i].updatePhenotype(AFFECTED);
            }
            __personsExtreme->push_back(__persons[i]);
          }
          else if (__persons[i].getPhenotype() <= sortedPhenotypes[lower]) {
            if (shouldMarkCaseCtrl) {
              __persons[i].updatePhenotype(UNAFFECTED);
            }
            __personsExtreme->push_back(__persons[i]);
          }
          else continue;
        }

        //!- Return to the extreme samples, a subset of the cohort.           
        __persons = *__personsExtreme;
        delete __personsExtreme;
      }
      break;   
  }
  
  //  std::cout << personTemplate->getPhenotype() << std::endl;
  delete personTemplate;

  return;
}


/*!\brief Generate Mendelian traits genotype and phenotype data*/
/*
void gwSimulator::createGenotypeMendelianTraitsAssociations (double percentageCausal, bool isAllelicHeterogeneous,
        const char moi, unsigned nCases, unsigned nHeterCases, unsigned nCtrls, gsl_rng* gslr)
{
  //!\brief his program is the core simulator function that simulates Mendelian traits by person's information (not pedigrees)



  //!- Initialize a "template" person object
  gwPerson* personTemplate = new gwPerson (__boundary, __neutral_cutoff, __pedInfos, __mafs, __fnctAnnotations, __positions);
  vectorUI mendelianMafIdxes = personTemplate->getMendelianMafIdxes(percentageCausal);
  personTemplate->updateLocusAttributes(mendelianMafIdxes, isAllelicHeterogeneous);

  //!- Generate case-ctrl samples (resample, reject or accept)
  unsigned iCase = 0, iCtrl = 0, iHeterCases = 0;

  while (iCase != nCases || iCtrl != nCtrls || iHeterCases != nHeterCases) {
    personTemplate->generateGenotype(0, gslr);
    double odds = personTemplate->computeGenotypicEffect(mendelianMafIdxes, isAllelicHeterogeneous, moi);
    personTemplate->generatePhenotype(odds, 0, gslr);
    double trait = personTemplate->getPhenotype();

    if (trait == AFFECTED) {
      // get a case
      if (iCase != nCases) {  
        // if case is still in need, collect it.
        __persons.push_back(*personTemplate);
        ++iCase;
      }
      else;  
      // if case is enough, do nothing.
    }
    else {
      // get a control
      if (iCtrl != nCtrls) {  
        // if control is still in need, collect it.
        __persons.push_back(*personTemplate);
        ++iCtrl;
      }
      else if (iHeterCases != nHeterCases) {
        personTemplate->updatePhenotype(AFFECTED);
        __persons.push_back(*personTemplate);
        ++iHeterCases;
      }  
      // if HeterCases is needed collect it.          
      else;
    }
  }
  delete personTemplate;
  return;
}
*/
void gwSimulator::createGenotypeMendelianTraitsAssociations (double percentageCausal, bool isAllelicHeterogeneous,
        const char moi, unsigned nCases, double propHeterCases, unsigned nCtrls, gsl_rng* gslr)
{  
  //!- Initialize a "template" person object
  gwPerson* personTemplate = new gwPerson (__boundary, __neutral_cutoff, __pedInfos, __mafs, __fnctAnnotations, __positions);
  vectorUI mendelianMafIdxes = personTemplate->getMendelianMafIdxes(percentageCausal);
  personTemplate->updateLocusAttributes(mendelianMafIdxes, isAllelicHeterogeneous);  

  bool isInputOk = (nCtrls > 0 && nCases > 0 && propHeterCases >= 0.0 && propHeterCases <= 1.0);
  if (isInputOk);
  else {
    std::cerr << "Case/Control/Prop-heterogenious input not valid. Quit now." << std::endl;
    exit(-1);
  }

  personTemplate->updateGenotypeFreqs(mendelianMafIdxes, isAllelicHeterogeneous, moi);
  vector2F genoFreqsRecovery = personTemplate->getGenotypeFreqs();


  personTemplate->updatePhenotype(AFFECTED); 
  
  for (unsigned i = 0; i != nCases; ++i) {
    double runif = gsl_rng_uniform(gslr);
    //!- Generate case-ctrl samples. No rejecting because case/ctrl use different underlying MAF
    if (runif >= propHeterCases) {
      if (moi != 'C') {
        personTemplate->updateGenotypeFreqs(mendelianMafIdxes, isAllelicHeterogeneous, moi, gslr); 
        personTemplate->generateGenotype(1, gslr);
        __persons.push_back(*personTemplate);
        personTemplate->updateGenotypeFreqs(genoFreqsRecovery);
      }
      else {
        vector2F tmpGenotype(2);
        personTemplate->updateGenotypeFreqs(mendelianMafIdxes, isAllelicHeterogeneous, moi, gslr); 
        personTemplate->generateGenotype(2, gslr);
        tmpGenotype[0] = personTemplate->getGenotype()[0];
        personTemplate->updateGenotypeFreqs(genoFreqsRecovery);
        personTemplate->updateGenotypeFreqs(mendelianMafIdxes, isAllelicHeterogeneous, moi, gslr); 
        personTemplate->generateGenotype(3, gslr);
        tmpGenotype[1] = personTemplate->getGenotype()[1];
        personTemplate->updateGenotype(tmpGenotype);
        __persons.push_back(*personTemplate);
        personTemplate->updateGenotypeFreqs(genoFreqsRecovery);
      }
    }
    //!- Generate unphenotyped cohorts. Use the MAF directly from SRV_batch
    else {
      personTemplate->generateGenotype(1, gslr);
      __persons.push_back(*personTemplate);
    } 
  }

  personTemplate->updatePhenotype(UNAFFECTED); 
  for (unsigned i = 0; i != nCtrls; ++i) {
    personTemplate->generateGenotype(1, gslr);
    __persons.push_back(*personTemplate);
  }
  
  //  std::cout << personTemplate->getPhenotype() << std::endl;
  delete personTemplate;
  return;
}


void gwSimulator::mimicGenotyping ( const vectorF& proportionsMissingData, const double missingLowMaf, 
    bool shouldMarkMissing, gsl_rng* gslr )
{
  //!\brief This program further edits the simulated sample by masking out some functional variants, etc.
  
  bool isInputOk = (__persons.size() != 0 && missingLowMaf < 1.0 && missingLowMaf >= 0.0);
  //!- Check input: make sure the matrix is square
  //!- (otherwise cannot compute MAF because locus information incomplete)
  for (unsigned i = 0; i != __persons.size(); ++i) {
    if (__persons[0].getLocusAttributes().size() != __persons[i].getGenotype()[0].size()) {
      isInputOk = false;
      break;
    }
  }

  if (isInputOk);
  else {
    std::cerr << "Input data-set not square or maf cut-offs not valid. Quit now." << std::endl;
    exit(-1);
  }

  //!- 1) Update locus attributes for each "person" object 2) Trim __mGenotypes 3) make data-set
  
  unsigned nVariantSites = __mafs.size();
  vectorL shouldTrimSites(nVariantSites, false);
  for (unsigned i = 0; i != shouldTrimSites.size(); ++i) {
    if (__mafs[i] < missingLowMaf) shouldTrimSites[i] = true;
  }

  
  for (unsigned i = 0; i != __persons.size(); ++i) {
    __persons[i].updateLocusAttributes(proportionsMissingData, shouldMarkMissing, gslr);
    __persons[i].updateLocusAttributes(shouldTrimSites, shouldMarkMissing);
    __persons[i].updateGenotype();
  }

  return;
}


void gwSimulator::createPedfileMatrix(bool isSynoTrimmed, bool isCvTrimmed, 
        bool isPedWritten, std::string projectName, std::string simulationTask, vectorF& keepedMafs, vectorL& keepedSites)
{
  //!\brief This program writes a ped file and make data matrix for assoc. analysis
  
  bool isInputOk = (__persons.size() != 0);
  //!- Check input: make sure the matrix is square. 
  //!- Even with missing data, to write a ped file the __mGenotypes should still be square
  for (unsigned i = 0; i != __persons.size(); ++i) {
    if (__persons[0].getLocusAttributes().size() != __persons[i].getGenotype()[0].size()) {
      isInputOk = false;
      break;
    }
  }

  if (isInputOk);
  else {
    std::cerr << "Input data-set not square / require genotype and phenotype vectors be pre-allocated / MAFs size should be 0. Quit now." << std::endl;
    exit(-1);
  }

  __mGenotypes.resize(__persons.size()); 
  __mPhenotypes.resize(__persons.size());

  vectorUI snvPositions(0);
  vector2F tmpGeno = __persons[0].getGenotype(false, isSynoTrimmed, isCvTrimmed, snvPositions, keepedMafs, keepedSites);
  unsigned  nVariantSites = tmpGeno[0].size();
  vectorUI tmpSnvPositions(0);
  __mObservedMafs.resize(nVariantSites, 0.0);

  //!- Make the data-matrix 
  //!- Calculate observed MAFs while making the matrix
  for (unsigned i = 0; i != __persons.size(); ++i) {
    
    __mPhenotypes[i] = __persons[i].getPhenotype();
    if (simulationTask != "5" && simulationTask != "4" && __mPhenotypes[i] == UNPHENOTYPED)
      __mPhenotypes[i] = UNAFFECTED;

    vectorF remainedMafs;
    vectorL nouse(__persons[i].getGenotype()[0].size(), false);
    tmpGeno = __persons[i].getGenotype(false, isSynoTrimmed, isCvTrimmed, tmpSnvPositions, remainedMafs, nouse);

    if (tmpGeno[0].size() != nVariantSites || tmpSnvPositions.size() != snvPositions.size()) {
      std::cerr << "Data-set should be square after trimmings of CV and synonymous sites. Quit now." << std::endl;
      exit(-1);
    }

    for (unsigned j = 0; j != nVariantSites; ++j) {
      if (tmpGeno[0][j] == MISSING_ALLELE || tmpGeno[1][j] == MISSING_ALLELE)
        __mGenotypes[i].push_back(MISSING_ALLELE);
      else {
        __mGenotypes[i].push_back( tmpGeno[0][j] + tmpGeno[1][j] );
        __mObservedMafs[j] += __mGenotypes[i][j];
      }
    }
  }
  
  for (unsigned j = 0; j != nVariantSites; ++j)
    __mObservedMafs[j] = __mObservedMafs[j] / (2.0 * __persons.size());


  //!- Write the data-matrix
  //!- 3 files to write: PED, MAP and LOG
  if (isPedWritten) {
    std::string pedFileName = projectName;
    pedFileName.append(".ped");
    std::string pedFileNameLog = projectName;
    pedFileNameLog.append(".log");    
    std::string mapFileName = projectName;
    mapFileName.append(".pos");
    std::ofstream fout(pedFileName.c_str());
    std::ofstream lout(pedFileNameLog.c_str());
    std::ofstream mout(mapFileName.c_str());
    fout.precision(5);
    lout.precision(5);

    for (unsigned i = 0; i != __persons.size(); ++i) {
      for (unsigned j = 0; j != 5; ++j) fout << 0 << " ";
      fout << __mPhenotypes[i] << " ";
      fout << __mGenotypes[i] << std::endl;
    }
    fout.close();
//    std::clog << "\tPED file written [ " << pedFileName <<  " ]." << std::endl;
    
    for (unsigned i = 0; i != snvPositions.size(); ++i) 
      mout << "M" << i + 1 << " " << snvPositions[i] << std::endl;
    mout.close();
//    std::clog << "\tMAP file written [ " << mapFileName <<  " ]." << std::endl;
    
    lout << "# " << pedFileName << " log file" << std::endl;
    lout << "Sample size:\n" << __persons.size() << std::endl;
    lout << "Length of variant sites:\n" << nVariantSites << std::endl;
    lout << "Are synonymous sites included in the data? (1 = Yes, 0 = No):\n" << (!isSynoTrimmed) << std::endl;
    lout << "Are (underlying) common variant sites included in the data? (1 = Yes, 0 = No):\n" << (!isCvTrimmed) << std::endl;
    lout << "Observed Locus MAF in the sample data:\n" << __mObservedMafs << std::endl;
    lout.close();
//    std::clog << "\tLog file written [ " << pedFileNameLog <<  " ]." << std::endl;
  }
  
  return;
}

/*
void gwSimulator::calcVariantsPars(std::string filename, double mafLower, double mafUpper) const 
{
  bool isInputOk = (__mObservedMafs.size() != 0 && __mGenotypes.size() != 0 && __mPhenotypes.size() == __mGenotypes.size());
  if (isInputOk);
  else {
    std::cerr << "Input data-set not valid. Cannot calculate PAR. Now Quit." << std::endl;
    exit(-1);
  }

  vectorF pars(__mObservedMafs.size());
  for (unsigned j = 0; j != pars.size(); ++j) {
    double a1 = 0.0, a2 = 0.0, u1 = 0.0, u2 = 0.0;
    for (unsigned i = 0; i != __mPhenotypes.size(); ++i) {
      if (__mPhenotypes[i] == AFFECTED) {
        if (__mGenotypes[i][j] != MAJOR_ALLELE && __mGenotypes[i][j] != MISSING_ALLELE)
          a1 += __mGenotypes[i][j];
        else 
          a2 += 2.0;
      }
      else if (__mPhenotypes[i] == UNAFFECTED) { 
        if (__mGenotypes[i][j] != MAJOR_ALLELE && __mGenotypes[i][j] != MISSING_ALLELE)
          u1 += __mGenotypes[i][j];
        else 
          u2 += 2.0;
      }
    }
    pars[j] = a1 / (a1 + u1) - a2 / (a2 + u2);
  }
  
  double par = 0.0;
  for (unsigned j = 0; j != pars.size(); ++j) {
    if (__mObservedMafs[j] <= mafLower || __mObservedMafs[j] > mafUpper)
      continue;
    else
      par += pars[j];
  }

  std::string parFileName = filename + "_PARS_allReplicates.log";
  std::ofstream outPar;
  outPar.open(parFileName.c_str(), std::ios::app);
  outPar << par << " ";
  return;
}
*/

vector2F gwSimulator::getGenotypes() const 
{
  return  __mGenotypes;
}
vectorF gwSimulator::getPhenotypes() const
{
  return __mPhenotypes;
}
vectorF gwSimulator::getObservedMafs() const
{
  return __mObservedMafs;
}
