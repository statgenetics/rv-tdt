#include "gw_utilities.h"
#include "person.h"
#include "assocsims.h"
#include "Argument_helper.h"
#include "mainPower.h"
#include "tdtTrio.h"
#include "rareTdt.h"
#include "phasedTDT.h"
#include "unphasedTDT.h"
using namespace std;

int main(int argc, const char* argv[])
{
  // Parameters
  unsigned seed = 0;
  string gFile("./example/kyrukov_sfscode");
  double boundary = 0.01;
  double neutral_cutoff = 0.0001;
  
  vector<double> propFunctionalRv(2);
  propFunctionalRv[0] = 1.0;
  //!- proportion of effective deleterious variant (vs. non-causal)
  propFunctionalRv[1] = 0.0;
  //!- proportion of effective protective variant (vs. non-causal)
  char moi = 'A';

  string projectName = "CompareExample";
  string simulationTask = "1";
  //!- options: 1.phased, 2.unphased

  vector<double> oddsRatios(5);
  oddsRatios[0] = 0.0;
  oddsRatios[1] = 2.0;
  //!- odds ratio for deleterious variants
  oddsRatios[2] = 0.0;
  oddsRatios[3] = 1.0;
  //!- odds ratio for protective variants
  oddsRatios[4] = 1.0;
  //!- odds ratio for common variants
  double baselinef = 0.01;
  //!- base-line penetrance ~ prevalence
  
  vector<double> propMissingData(4);
  propMissingData[0] = 0.0;
  propMissingData[1] = 0.0;
  propMissingData[2] = 0.0;
  propMissingData[3] = 0.0;
  double missingLowMaf = 0.0;
  bool shouldMarkMissing = false;

  unsigned nCases = 500;
  unsigned nCtrls = 0;

  /*** Write simulated data ***/
  bool isPedWritten = false;
  bool isSynoKept = false;
  bool isCvTrimmed = false;
  vector<double> pedInfos(6, 0.0);

  /*** settings for analysis ***/
  //!-Trim data parameters
  double mafLower = 0;
  double mafUpper = 1.0;
  //double alpha = 0.05;
  unsigned nPermutations = 2000;
  unsigned nReplicates = 50;
  // to simulate sequencing missing
  double indMiss(0), siteMiss(0);
  // analyze missing
  double missingratio(0);
  bool skipMiss(true);

  // Mz permut for rareTdt class
  bool MzPermut(false);
    
  // argsparser
  dsr::Argument_helper ah;
  // the only required argument
  ah.new_string("task", "\n\tAnalysis task. \n\t| Type values \"1~2\"\n\t| 1: phased \n\t| 2: unphased \n\t", simulationTask);
  ah.new_optional_string("pname", "\n\tProject name. \n\t| STRING \n\t| set output file names. \n\t", projectName);
  ah.new_optional_string("gdata", "\n\tGenetic data files for the simulation to be based on. \n\t| STRING\n\t| Proper gdata.maf, gdata.sel and gdata.pos files need to be provided to the program\n\t", gFile);
  // named arguments
  ah.new_named_double('A', "define_rare", "<frequency>", "\n\tDefinition of rare variants. \n\t| FLOAT\n\t| variants having MAF <= frequency will be considered as \"rare\" variants\n\t", boundary);
  ah.new_named_double('B', "prop_func_deleterious", "<fraction>", "\n\tProportion of FUNCTIONAL deleterious variants. \n\t| FLOAT\n\t| 0.0 <= fraction <= 1.0\n\t| (1-proportion)x100\% is the proportion of non-causal deleterious variants (noise) \n\t| This is NOT the proportion of deleterious variants, which should have been defined in <gdata.sel> file\n\t", propFunctionalRv[0]);
  ah.new_named_char('D', "mode_of_inheritance", "<moi>", "\n\tMode of inheritance under which the phenotype data is simulated.\n\t| CHAR\n\t| \"A\" (additive), \"D\" (dominant), \"R\" (recessive), \"M\" (multiplicative), \n\t| \"C\" (compound dominant for non-mendelian traits, or compound recessive for mendelian traits)\n\t", moi);
  //
  ah.new_named_double('E', "OR_deleterious_min", "<effect_size>", "\n\t[task=1,2] Minimum odds ratio for deleterious variants. \n\t| FLOAT\n\t| 0.0 for fixed effect size model; >=1.0 for variable effect sizes model\n\t", oddsRatios[0]);
  ah.new_named_double('F', "OR_deleterious_max", "<effect_size>", "\n\t[task=1,2] Maximum odds ratio for deleterious variants. \n\t| FLOAT\n\t| >=1.0 AND > $OR_deleterious_min\n\t| will be the odds ratio for deleterious variants in fixed effect size model, or the maximum odds ratio in variable effect sizes model\n\t", oddsRatios[1]);
  ah.new_named_double('I', "OR_common", "<effect_size>", "\n\t[task=1,2] Odds ratio for common variants. \n\t| FLOAT\n\t| =1.0 for neutral, >1.0 for deleterious, <1.0 for protective \n\t| this is the odds ratio for all variants having MAF > $define_rare, i.e., the \"common\" variants \n\t| assuming all common variants have fixed effect size\n\t", oddsRatios[4]);
  ah.new_named_double('J', "prevalence", "<fraction>", "\n\t[task=1,2] Disease prevalence. \n\t| FLOAT \n\t| will be used as baseline penetrance of a gene (baseline penetrance ~= disease prevalence)\n\t", baselinef);
  //
  ah.new_flag('Q', "keep_synonymous", "\n\tKeep synonymous variants from analysis. \n\t| BOOLEAN\n\t", isSynoKept);
  ah.new_flag('R', "remove_common_loci", "\n\tRemove common variant sites from analysis. \n\t| BOOLEAN\n\t| the \"common\" loci refers to variants in the haplotype pool having MAF > $define_rare \n\t", isCvTrimmed);
  ah.new_named_unsigned_int('S', "num_cases", "<#cases>", "\n\t[task=1,3,5,6] Number of cases. \n\t| INT (>0)\n\t| for extreme QT simulations (task=5) it will be #samples having high QT values from the population when $num_all_samples is set to 0\n\t", nCases);

  ah.new_named_double('b', "define_neutral", "<annotation_cutoff>", "\n\tAnnotation value cut-off that defines a variant to be \"neutral\" in evolution (either synonymous or non-coding). \n\t| FLOAT\n\t| loci having annotation value (in gdata.ann) between (-$annotation_cutoff, +$annotation_cutoff) will be regarded \"neutral\" and will not contribute to phenotype \n\t", neutral_cutoff);
  ah.new_named_double('c', "individ_miss", "<fraction>", "\n\tThe proportion of individual that has entire missing, to mimic sample collecting process. \n\t| FLOAT\n\t", indMiss);
  ah.new_named_double('d', "site_miss", "<fraction>", "\n\tThe proportion of missing on each site, to mimic sequencing process. \n\t| FLOAT\n\t", siteMiss);
  ah.new_flag('e', "trim_miss", "\n\tHow to analyze missing site. \n\t| BOOLEAN\n\t| Skip incomplete trio by default, if envoked, will analyze those sites that have missing in a trio\n\t", skipMiss);
  ah.new_named_double('m', "missing_ratio", "<fraction>", "\n\tThe max proportion of missing on each site, to mimic data cleanning process. \n\t| FLOAT\n\t", missingratio);
  //
  ah.new_named_double('f', "maf_lower", "<frequency>", "\n\tLower bound of observed sample minor allele frequency. \n\t| FLOAT\n\t| loci having observed MAF < $maf_lower will not be analyzed \n\t", mafLower);
  ah.new_named_double('g', "maf_upper", "<frequency>", "\n\tUpper bound of observed sample minor allele frequency. \n\t| FLOAT\n\t| loci having observed MAF > $maf_upper will not be analyzed \n\t", mafUpper);
  //ah.new_named_double('i', "significance", "<alpha_level>", "\n\tSignificance level at which power will be evaluated. \n\t| FLOAT\n\t", alpha);
  ah.new_named_unsigned_int('j', "replicates", "<#replicates>", "\n\tNumber of replicates for power evaluation. \n\t| INT (>0)\n\t", nReplicates);
  ah.new_named_unsigned_int('k', "permutations", "<#permutations>", "\n\tNumber of permutations, only applicable to permutation based methods. \n\t| INT (>0)\n\t", nPermutations);
  ah.new_named_unsigned_int('l', "rng_seed", "<long_integer>", "\n\tSeed for random number generator. \n\t| INT (>=0) \n\t| =0 is to use a random seed (seed = system time + process ID)\n\t", seed);
  ah.new_flag('p', "MzPermut", "\n\tprint out MZ statistic using genotype and hapotype permutation in rareTdt class. \n\t| BOOLEAN\n\t", MzPermut);

  // program information
  string program_name = "compare";
  double version = 0.0;
  string banner = "\n\t:--------------------------------------------------------:\n\t: Power Calculator for Transmission disequilibrium Tests :\n\t:--------------------------------------------------------:\n\t: (c) 2012 Zongxiao He  |  http://bcm.edu/genetics/leal  :\n\t:--------------------------------------------------------:\n";
  ah.set_name(program_name.c_str());
  ah.set_description(banner.c_str());
  ah.set_version(version);
  ah.set_author("Zongxiao He <zongxiah@bcm.edu>");
  time_t rawtime;
  struct tm * timeinfo;
  time ( &rawtime );
  timeinfo = localtime ( &rawtime );
  ah.set_build_date(asctime (timeinfo));

  ah.process(argc, argv);
  std::vector<std::string> arguments(argv + 1, argv + argc);
  std::string cmdcurrent = boost::algorithm::join(arguments, "_");
  
  bool isSynoTrimmed = !isSynoKept;
  string mafFile = gFile + ".maf";
  string selFile = gFile + ".ann";
  string posFile = gFile + ".pos";

  // Power calculations
  RNG rng;
  gsl_rng* gslr;
  if (seed == 0) {
    gslr = rng.get();
  }
  else {
    gslr = rng.get(seed);
  }

  vector2F mafDat;
  vector2F selDat;
  vector2UI posDat;
  scan_vector2F(selFile, selDat);
  scan_vector2F(mafFile, mafDat);
  scan_vector2UI(posFile, posDat);
  
  // For MzPermut: set random seed, and only one replicates
  if(MzPermut) {gslr = rng.get(); nReplicates = 1;}

  // all output goes to std::cout
  unsigned iReplicate = 0;
  while (iReplicate != nReplicates) {
    /*** Choose MAF data ***/
    unsigned dataIdx = gsl_rng_uniform_int(gslr, mafDat.size());
    vector<double>& mafs =  mafDat[dataIdx];
    vector<double>& selcoefs = selDat[dataIdx];
    vector<unsigned>& positions = posDat[dataIdx];
    vector2UI poolDat;
    poolDat.resize(0);
    
    vector2F genoFreqs(mafs.size());
    for (unsigned i = 0; i != mafs.size(); ++i) {
      genoFreqs[i].push_back( (1 - mafs[i]) * (1 - mafs[i]) );
      genoFreqs[i].push_back( mafs[i] * mafs[i] );
    }

    /*** Generate population data ***/
    gwSimulator* simulator = new gwSimulator(boundary, neutral_cutoff, pedInfos, mafs, genoFreqs, selcoefs, positions);
    //simulator->createGenotypeComplexTraitsAssociations (propFunctionalRv, baselinef, oddsRatios, moi, nCases, nCtrls, nUnphenotyped, gslr, dsr::verbose, projectName); // "dichot-odds"
    simulator->createGenotypeComplexTraitsAssociations (propFunctionalRv, baselinef, oddsRatios, moi, nCases, nCtrls, 0, gslr, false, projectName, poolDat); // "dichot-odds"
    /*** Exclude some variants (a mimic genotyping procedure) ***/
    simulator->mimicGenotyping( propMissingData, missingLowMaf, shouldMarkMissing, gslr);
    vectorF keepedMafs;
    vectorL keepedSites(mafs.size(), false);
    //default: isSynoTrimmed = true; isCvTrimmed = false; isPedWritten = false;
    //modify: person.cpp line 1564, person.h 236; associms.cpp 585; assocsims.h 102
    simulator->createPedfileMatrix(isSynoTrimmed, isCvTrimmed, isPedWritten, projectName, "dichot-odds", keepedMafs, keepedSites);
    
    // need to make sure mafs has same variant sites as genotypes
    vector2F genotypes = simulator->getGenotypes();
    if(keepedMafs.size() != genotypes[0].size()){
      std::cerr << "ERROR: unmatched maf data and genotype data. after trim. Now Quit." << endl;
      exit(-1);
    }
    delete simulator;
    // simulate sequencing process, introduce missing
    vector2F genotypesTDT = zx_genGenos(genotypes, keepedMafs, gslr);
    // sequencing --> missing
    if(indMiss > 0 || siteMiss > 0) zx_createMissing(genotypesTDT, indMiss, siteMiss, gslr);
      
    std::string geneName = "Gene_" + n2s(iReplicate);
    std::string projName = projectName + "_TDT";
    bool phased = (simulationTask == "1")?true:false;
    if(simulationTask == "2") genotypesTDT = zx_phased2unphased(genotypesTDT);
    rareTdt tdtObj(projName, geneName, phased, genotypesTDT, keepedMafs, skipMiss, missingratio);

    /*** FOR ME_PERMUT OPTION ***/
    if(MzPermut){
      vectorF pval = tdtObj.tdtTest("MzPermut", mafLower, mafUpper, nPermutations*2, nPermutations, 0.05, false); pval.push_back(0.0);
      iReplicate++;
      continue;
    }

    vectorF pval_New = tdtObj.tdtTest("All", mafLower, mafUpper, nPermutations*2, nPermutations, 0.05, false);
    //order: mz, cmc, VT-MZ-Geno, VT-CMC-Geno (if phased), WSS-Geno, VT-MZ-Hapo, VT-CMC-Hapo (if phased), WSS-Hapo
    rareTdt tdtObj2(projName, geneName, phased, genotypesTDT, keepedMafs, skipMiss, missingratio);
    vectorF pval_Old = tdtObj2.tdtTest("All", mafLower, mafUpper, nPermutations*2, nPermutations, 0.05, false);
    
    //std::string NotX = "NotX";
    //vectorF pval_Old1, pval_Old2, pval_Old3, pval_Old;
    /*** Phased TDT TESTs ***/
    //if(simulationTask == "1") { 
	    //phasedTDT tdtTest(genotypesTDT, keepedMafs, keepedMafs, NotX, missingratio, skipMiss);
	    //pval_Old1 = tdtTest.TdtMZCMC(mafLower, mafUpper);
	    //pval_Old2 = tdtTest.TdtPermut(0, 1, nPermutations*2, nPermutations, 0.05, "GenoShuffle"); // no cutoff for WSS and VT
	    //phasedTDT tdtTest2(genotypesTDT, keepedMafs, keepedMafs, NotX, missingratio, skipMiss);
	    //pval_Old3 = tdtTest2.TdtPermut(0, 1, nPermutations*2, nPermutations, 0.05, "HapoShuffle"); // return WSS VT-MZ VT-CMC
	    //pval_Old.push_back(pval_Old1[0]); pval_Old.push_back(pval_Old1[2]);
	    //pval_Old.push_back(pval_Old2[1]); pval_Old.push_back(pval_Old2[2]); pval_Old.push_back(pval_Old2[0]);
	    //pval_Old.push_back(pval_Old3[1]); pval_Old.push_back(pval_Old3[2]); pval_Old.push_back(pval_Old3[0]);
    //}
    //else if(simulationTask == "2") {
	    ///*** unPhased TDT TESTs ***/
	    //unphasedTDT tdtTest(genotypesTDT, keepedMafs, keepedMafs, NotX, missingratio, skipMiss);
	    //pval_Old1 = tdtTest.TdtMZ(mafLower, mafUpper, false);
	    //pval_Old2 = tdtTest.TdtPermut(0, 1, nPermutations*2, nPermutations, 0.05, true, "GenoShuffle"); // no cutoff for WSS and VT
	    //unphasedTDT tdtTest2(genotypesTDT, keepedMafs, keepedMafs, NotX, missingratio, skipMiss);
	    //pval_Old3 = tdtTest2.TdtPermut(0, 1, nPermutations*2, nPermutations, 0.05, true, "HapoShuffle"); // return WSS and VT-MZ
	    //pval_Old.push_back(pval_Old1[0]);
	    //pval_Old.push_back(pval_Old2[1]); pval_Old.push_back(pval_Old2[0]);
	    //pval_Old.push_back(pval_Old3[1]); pval_Old.push_back(pval_Old2[0]);
    //}
    
    //in case wrong p value
    bool informative(true);
    for(UINT i=0; i<pval_New.size(); i++){if(pval_New[i] == std::numeric_limits<double>::infinity()) informative = false;}
    for(UINT i=0; i<pval_Old.size(); i++){if(pval_Old[i] == std::numeric_limits<double>::infinity()) informative = false;}
    if(!informative){
      std::cout << "\tuninformative dataset\n--------\n";
      continue;
    }
    
    if(iReplicate == 0) {
      if(phased) std::cout << "========\nMZ CMC VT-MZ-Geno VT-CMC-Geno WSS-Geno VT-MZ-Hapo VT-CMC-Hapo WSS-Hapo\n========\n";
      else std::cout << "========\nMZ VT-MZ-Geno WSS-Geno VT-MZ-Hapo WSS-Hapo\n========\n";
    }

    std::cout << "New: " << pval_New << "\nOld: " << pval_Old << std::endl << "========\n";      
    ++iReplicate;
  }
  return 0;
}
