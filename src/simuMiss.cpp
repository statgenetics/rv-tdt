#include "gw_utilities.h"
#include "person.h"
#include "assocsims.h"
#include "Argument_helper.h"
#include "mainPower.h"
#include "tdtTrio.h"
#include "rareTdt.h"
using namespace std;

int main(int argc, const char* argv[])
{
    // Parameters
    unsigned seed = 0;
    string gFile("Boyko2008European1p5k");
    //bool shouldUseGenPool = false;
    double boundary = 0.01;
    double neutral_cutoff = 0.0001;

    vector<double> propFunctionalRv(2);
    propFunctionalRv[0] = 1.0;
    //!- proportion of effective deleterious variant (vs. non-causal)
    propFunctionalRv[1] = 0.0;
    //!- proportion of effective protective variant (vs. non-causal)
    char moi = 'A';

    string projectName = "JustForFun";
    string simulationTask = "1";
    //!- options: 1.simulated phased data, 2.simulated unphased data

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
    string test = "CMC-one";
    double mafLower = 0.0;
    double mafUpper = 1.0;
    double alpha = 0.05;
    unsigned nPermutations = 2000;
    unsigned adaptive = 1000;
    unsigned nReplicates = 1000;
    // to simulate sequencing missing
    double indMiss(0), siteMiss(0);
    // analyze missing
    double missingratio(0);
    //bool skipMiss(true);
    //bool usePOPmaf(true);

    // write information about every replicate
    bool recordGenos(false);
    bool recordInfos(true);
    bool recordStatic(true);
    // record informative replicates
    UINT unAnalyzedGens(0);

    // argsparser
    dsr::Argument_helper ah;
    // the only required argument
    ah.new_string("task", "\n\tAnalysis task. \n\t| Type values \"1~2\"\n\t| 1: Simulated phased data \n\t| 2: Simulated unphased data\n\t", simulationTask);
    ah.new_string("pname", "\n\tProject name. \n\t| STRING \n\t| set output file names. \n\t", projectName);
    ah.new_optional_string("gdata", "\n\tGenetic data files for the simulation to be based on. \n\t| STRING\n\t| Proper gdata.maf, gdata.sel and gdata.pos files need to be provided to the program\n\t", gFile);
    // named arguments
    ah.new_named_double('A', "define_rare", "<frequency>", "\n\tDefinition of rare variants. \n\t| FLOAT\n\t| variants having MAF <= frequency will be considered as \"rare\" variants\n\t", boundary);
    ah.new_named_double('B', "prop_func_deleterious", "<fraction>", "\n\tProportion of FUNCTIONAL deleterious variants. \n\t| FLOAT\n\t| 0.0 <= fraction <= 1.0\n\t| (1-proportion)x100\% is the proportion of non-causal deleterious variants (noise) \n\t| This is NOT the proportion of deleterious variants, which should have been defined in <gdata.sel> file\n\t", propFunctionalRv[0]);
    ah.new_named_double('C', "prop_func_protective", "<fraction>", "\n\tProportion of FUNCTIONAL protective variants. \n\t| FLOAT\n\t| 0.0 <= fraction <= 1.0\n\t| (1-proportion)x100\% is the proportion of non-causal protective variants (noise) \n\t| This is NOT the proportion of protective variants, which should have been defined in <gdata.sel> file\n\t", propFunctionalRv[1]);
    ah.new_named_char('D', "mode_of_inheritance", "<moi>", "\n\tMode of inheritance under which the phenotype data is simulated.\n\t| CHAR\n\t| \"A\" (additive), \"D\" (dominant), \"R\" (recessive), \"M\" (multiplicative), \n\t| \"C\" (compound dominant for non-mendelian traits, or compound recessive for mendelian traits)\n\t", moi);
    //
    ah.new_named_double('E', "OR_deleterious_min", "<effect_size>", "\n\t[task=1,2] Minimum odds ratio for deleterious variants. \n\t| FLOAT\n\t| 0.0 for fixed effect size model; >=1.0 for variable effect sizes model\n\t", oddsRatios[0]);
    ah.new_named_double('F', "OR_deleterious_max", "<effect_size>", "\n\t[task=1,2] Maximum odds ratio for deleterious variants. \n\t| FLOAT\n\t| >=1.0 AND > $OR_deleterious_min\n\t| will be the odds ratio for deleterious variants in fixed effect size model, or the maximum odds ratio in variable effect sizes model\n\t", oddsRatios[1]);
    ah.new_named_double('G', "OR_protective_min", "<effect_size>", "\n\t[task=1,2] Minimum odds ratio for protective variants. \n\t| FLOAT\n\t| 0.0 for fixed effect size model; <1.0 for variable effect sizes model\n\t", oddsRatios[2]);
    ah.new_named_double('H', "OR_protective_max", "<effect_size>", "\n\t[task=1,2] Maximum odds ratio for protective variants. \n\t| FLOAT\n\t| <=1.0 AND > $OR_protective_min\n\t| will be the odds ratio for protective variants in fixed effect size model, or the maximum odds ratio in variable effect sizes model\n\t", oddsRatios[3]);
    ah.new_named_double('I', "OR_common", "<effect_size>", "\n\t[task=1,2] Odds ratio for common variants. \n\t| FLOAT\n\t| =1.0 for neutral, >1.0 for deleterious, <1.0 for protective \n\t| this is the odds ratio for all variants having MAF > $define_rare, i.e., the \"common\" variants \n\t| assuming all common variants have fixed effect size\n\t", oddsRatios[4]);
    ah.new_named_double('J', "prevalence", "<fraction>", "\n\t[task=1,2] Disease prevalence. \n\t| FLOAT \n\t| will be used as baseline penetrance of a gene (baseline penetrance ~= disease prevalence)\n\t", baselinef);
    //
    ah.new_named_double('K', "prop_missing_deleterious", "<fraction>", "\n\tProportion of deleterious variants missing data. \n\t| FLOAT \n\t| missing genotypes will be coded as wildtype genotype by default\n\t", propMissingData[0]);
    ah.new_named_double('L', "prop_missing_protective", "<fraction>", "\n\tProportion of protective variants missing data. \n\t| FLOAT \n\t| missing genotypes will be coded as wildtype genotype by default\n\t", propMissingData[1]);
    ah.new_named_double('M', "prop_missing_non_causal", "<fraction>", "\n\tProportion of non-causal variants missing data. \n\t| FLOAT \n\t| missing genotypes will be coded as wildtype genotype by default\n\t", propMissingData[2]);
    ah.new_named_double('N', "prop_missing_synonymous", "<fraction>", "\n\tProportion of synonymous variants missing data. \n\t| FLOAT \n\t| missing genotypes will be coded as wildtype genotype by default\n\t", propMissingData[3]);
    ah.new_named_double('O', "missing_low_maf", "<frequency>", "\n\tVariants having MAF < $missing_low_maf will be marked as missing data. \n\t| FLOAT \n\t| note that $missing_low_maf is compared against the haplotype pool, not the sample\n\t", missingLowMaf);
    ah.new_flag('P', "recode_missing", "\n\tRe-code missing data. \n\t| BOOLEAN\n\t| if envoked, will re-code missing data from wildtype genotype to \"-9\", indicating missingness\n\t", shouldMarkMissing);
    ah.new_flag('Q', "keep_synonymous", "\n\tKeep synonymous variants from analysis. \n\t| BOOLEAN\n\t", isSynoKept);
    ah.new_flag('R', "remove_common_loci", "\n\tRemove common variant sites from analysis. \n\t| BOOLEAN\n\t| the \"common\" loci refers to variants in the haplotype pool having MAF > $define_rare \n\t", isCvTrimmed);
    ah.new_named_unsigned_int('S', "num_cases", "<#cases>", "\n\t[task=1,3,5,6] Number of cases. \n\t| INT (>0)\n\t| for extreme QT simulations (task=5) it will be #samples having high QT values from the population when $num_all_samples is set to 0\n\t", nCases);

    //ah.new_flag('a', "use_haplotype_pool", "\n\t[task=1,2,4,5] Randomly sample haplotypes from haplotype pool file $gdata.hap, rather than generating haplotypes on the fly. \n\t| BOOLEAN \n\t", shouldUseGenPool);
    ah.new_named_double('b', "define_neutral", "<annotation_cutoff>", "\n\tAnnotation value cut-off that defines a variant to be \"neutral\" in evolution (either synonymous or non-coding). \n\t| FLOAT\n\t| loci having annotation value (in gdata.ann) between (-$annotation_cutoff, +$annotation_cutoff) will be regarded \"neutral\" and will not contribute to phenotype \n\t", neutral_cutoff);
    ah.new_named_double('c', "individ_miss", "<fraction>", "\n\tThe proportion of individual that has entire missing, to mimic sample collecting process. \n\t| FLOAT\n\t", indMiss);
    ah.new_named_double('d', "site_miss", "<fraction>", "\n\tThe proportion of missing on each site, to mimic sequencing process. \n\t| FLOAT\n\t", siteMiss);
    //ah.new_flag('e', "skip_miss", "\n\tHow to analyze missing site. \n\t| BOOLEAN\n\t| Skip incomplete trio by default, if envoked, will analyze those sites that have missing in a trio\n\t", skipMiss);
    //ah.new_flag('s', "use_sMaf", "\n\tUse sample maf to infer miss \n\t| BOOLEAN\n\t| Use provided population maf by default, if envoked, will use calcualted sample maf to infer miss and permut incomplete trio\n\t", usePOPmaf);
    ah.new_named_double('m', "missing_ratio", "<fraction>", "\n\tThe max proportion of missing on each site, to mimic data cleanning process. \n\t| FLOAT\n\t", missingratio);
    //ah.new_named_double('n', "prop_non_causal", "<fraction>", "\n\tThe proporation of non-causal variants that are included in analysis. \n\t| FLOAT\n\t| They are added after gw_simulator generates genotype, the sites are randomly selected from keepedmafs \n\t", propNonCausal);
    //
    ah.new_named_double('f', "maf_lower", "<frequency>", "\n\tLower bound of observed sample minor allele frequency. \n\t| FLOAT\n\t| loci having observed MAF < $maf_lower will not be analyzed \n\t", mafLower);
    ah.new_named_double('g', "maf_upper", "<frequency>", "\n\tUpper bound of observed sample minor allele frequency. \n\t| FLOAT\n\t| loci having observed MAF > $maf_upper will not be analyzed \n\t", mafUpper);
    ah.new_named_string('t', "test", "<tdt_test>", "\n\tTDT test method. \n\t| STRING \n\t| All, GenoPermut, HapoPermut, NoPermut\n\t", test);
    ah.new_named_double('i', "significance", "<alpha_level>", "\n\tSignificance level at which power will be evaluated. \n\t| FLOAT\n\t", alpha);
    ah.new_named_unsigned_int('j', "replicates", "<#replicates>", "\n\tNumber of replicates for power evaluation. \n\t| INT (>0)\n\t", nReplicates);
    ah.new_named_unsigned_int('k', "permutations", "<#permutations>", "\n\tNumber of permutations, only applicable to permutation based methods. \n\t| INT (>0)\n\t", nPermutations);
    ah.new_named_unsigned_int('a', "adaptive", "<int>", "\n\tNumber of permutations for adaptive check,  only applicable to permutation based methods.\n\t", adaptive);
    ah.new_named_unsigned_int('l', "rng_seed", "<long_integer>", "\n\tSeed for random number generator. \n\t| INT (>=0) \n\t| =0 is to use a random seed (seed = system time + process ID)\n\t", seed);
    ah.new_flag('U', "record_genos", "\n\tIf or not record genotypes in each replicates. \n\t| BOOLEAN\n\t| if envoked, write down genotypes in each replicates\n\t", recordGenos);
    ah.new_flag('V', "record_infos", "\n\tIf or not record genotypes in each replicates. \n\t| BOOLEAN\n\t| if envoked, write down information in each replicates\n\t", recordInfos);
    ah.new_flag('Z', "record_static", "\n\tIf or not record rare variant static summary in each replicates. \n\t| BOOLEAN\n\t| if envoked, write down information for each replicates\n\t", recordStatic);

    //ah.new_flag('Z', "simulation_only", "\n\tWrite out simulated genotype-phenotype data to file $pname, rather than calculating power. \n\t| BOOLEAN \n\t", isPedWritten);

    // program information
    string program_name = "tdtSimu";
    double version = 2.1;
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

    std::string pvalFile = projectName + ".pval";
    FILE * pValout;
    pValout = fopen(pvalFile.c_str(), "a");

    std::string StatFile = projectName + ".stat";
    FILE * statout;
    statout = fopen(StatFile.c_str(), "a");

    std::string GenoDir = projectName + "_genos";
    std::string dirname = makedir(GenoDir);

    //return 8 pvalues at most
    vectorUI iPowerful(20,0);
    UINT testNum(0);
    unsigned iReplicate = 0;

    //record rare variant information//
    vector2F RaraVariantSummary;

    while (iReplicate != nReplicates) {
        /*** Choose MAF data ***/
        unsigned dataIdx = gsl_rng_uniform_int(gslr, mafDat.size());
        vector<double>& mafs =  mafDat[dataIdx];
        vector<double>& selcoefs = selDat[dataIdx];
        vector<unsigned>& positions = posDat[dataIdx];
        vector2UI poolDat;
        poolDat.resize(0);
        //if (shouldUseGenPool) {
        //  scan_vector2UI(gFile + "_hap/" + "hapo" + n2s(dataIdx+1), poolDat);
        //} else {
        //  poolDat.resize(0);
        //}

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
        vector2F genotypeTDT = zx_genGenos(genotypes, keepedMafs, gslr);
        vectorL unInformSites(keepedMafs.size(), true);
        for(UINT i=0; i<genotypeTDT[0].size(); i++){
            for(UINT j=0; j<genotypeTDT.size(); j++){
                // keep this site if any trio is informative
                if(genotypeTDT[j][i]>0) {unInformSites[i]=false; break;}
            }
        }
        vectorF newMafs; vector2F genotypesTDT;
        for(UINT i=0; i<unInformSites.size(); i++) {if(!unInformSites[i]) newMafs.push_back(keepedMafs[i]);}
        for(UINT i=0; i<genotypeTDT.size(); i++) {
            vectorF oneHapo;
            for(UINT j=0; j<genotypeTDT[i].size(); j++) {if(!unInformSites[j]) oneHapo.push_back(genotypeTDT[i][j]);}
            genotypesTDT.push_back(oneHapo);
        }
        keepedMafs = newMafs;
        // sequencing --> missing
        if(indMiss > 0 || siteMiss > 0) zx_createMissing(genotypesTDT, indMiss, siteMiss, gslr);
        /*** rare variant summary ***/
        RaraVariantSummary.push_back(zx_vStatSum(dataIdx, boundary, mafs, keepedMafs, genotypesTDT));
        /*** TDT TESTs ***/
        std::string geneName = "Gene_" + n2s(iReplicate);
        std::string projName = projectName + "_TDT";
        bool phased = (simulationTask == "1")?true:false;
        if(simulationTask == "2") genotypesTDT = zx_phased2unphased(genotypesTDT);

        rareTdt tdtObj(projName, geneName, phased, genotypesTDT, keepedMafs, true, missingratio);
        rareTdt tdtObjPopInfer(projName, geneName, phased, genotypesTDT, keepedMafs, false, missingratio);
        rareTdt tdtObjSamInfer(projName, geneName, phased, genotypesTDT, false, missingratio);

        vectorF pvalue = tdtObj.tdtTest(test, mafLower, mafUpper, adaptive, nPermutations, 0.05);
        vectorF pvaluePop = tdtObjPopInfer.tdtTest(test, mafLower, mafUpper, adaptive, nPermutations, 0.05);
        vectorF pvalueSam = tdtObjSamInfer.tdtTest(test, mafLower, mafUpper, adaptive, nPermutations, 0.05);
	for(UINT i=0; i<pvaluePop.size(); i++) pvalue.push_back(pvaluePop[i]);
	for(UINT i=0; i<pvalueSam.size(); i++) pvalue.push_back(pvalueSam[i]);

        //in case wrong p value
        bool informative(true);
        for(UINT i=0; i<pvalue.size(); i++){if(pvalue[i] == std::numeric_limits<double>::infinity()) informative = false;}
        if(!informative){
            unAnalyzedGens++; continue;
        }

        testNum = pvalue.size();
        for(UINT i=0; i<pvalue.size(); i++) {if (pvalue[i] <= alpha) iPowerful[i]++;}

        flockfile(pValout);
        fprintf(pValout, "%d\t", dataIdx+1);
        for(UINT i=0; i<pvalue.size(); i++) {fprintf(pValout, "%f\t", pvalue[i]);}
        fprintf(pValout, "\n");
        funlockfile(pValout);

        //record genotypes or informations//
        if(recordInfos) {
            std::string ldirname = makedir(projName + "_log");
            std::string logname = "./" + ldirname + "/" + geneName + ".log";
            std::ofstream log;
            log.open(logname.c_str());
            Json::StyledStreamWriter writer;
            writer.write(log, tdtObj.getJsonLog());
            log.close();
        }

        if(recordGenos){
            std::string genoFile = "./" + GenoDir + "/geno" + n2s(iReplicate);
            FILE * Genout;
            Genout = fopen(genoFile.c_str(), "w");
            for(UINT i=0; i<genotypesTDT.size(); i++){
                for(UINT j=0; j<genotypesTDT[i].size(); j++) fprintf(Genout, "%d\t", int(genotypesTDT[i][j]));
                fprintf(Genout, "\n");
            }
            fclose(Genout);
        }

        ++iReplicate;
    }

    vectorF power;
    for(UINT i=0; i<testNum; i++) power.push_back(1.0*iPowerful[i]/(1.0*nReplicates));

    vectorF RaraVariantStatic(RaraVariantSummary[0].size(), 0.0);
    for(UINT i=0; i<RaraVariantSummary.size(); i++){
        for(UINT j=0; j<RaraVariantSummary[i].size(); j++){
            RaraVariantStatic[j] += RaraVariantSummary[i][j];
        }
    }
    for(UINT i=0; i<RaraVariantStatic.size(); i++) RaraVariantStatic[i]=RaraVariantStatic[i]/(RaraVariantSummary.size()); 

    if(recordStatic){
        flockfile(statout);
        for(UINT i=0; i<RaraVariantSummary.size(); i++){
            fprintf(statout, "%d ",int(RaraVariantSummary[i][0]));
            for(UINT j=1; j<RaraVariantSummary[i].size(); j++){
                fprintf(statout, "%-6.4f ", RaraVariantSummary[i][j]);
            }
            fprintf(statout, "\n");
        }
        funlockfile(statout);
    }

    //if (!dsr::quiet) std::clog << std::endl << "INFO: output format header [ METHOD|POWER|PROJECT_ID|COMMAND ]" << std::endl << test << "|" << power << "|" << projectName << "|" << cmdcurrent << std::endl;
    std::string logFile = projectName + ".log";
    FILE * logout;
    logout = fopen(logFile.c_str(), "a");
    flockfile(logout);

    fprintf(logout, "#projName gFile nCases boundary neutral_cutoff ORlow ORhigh prop mafUpper skipMiss usePOPmaf siteMiss oriVarSites  AnalyzedSite  nVarInTrio  nVarInParent  nVarInKid  nCommonVariantInTrio  nTrioHasVariantInParent  nTrioHasVariantInKid  MissInParent  MissInSite test power1 power2 power3 (power4)\n");
    fprintf(logout, "%s %s %d %-4.2f %-7.6f %-3.1f %-3.1f %-4.2f %-4.2f %d %d %-4.2f ", projectName.c_str(), gFile.c_str(), nCases, boundary, neutral_cutoff, oddsRatios[0], oddsRatios[1], propFunctionalRv[0], mafUpper, true, true, siteMiss);
    for(UINT i=1; i<RaraVariantStatic.size(); i++){
        fprintf(logout, "%-6.4f ", RaraVariantStatic[i]);
    } //write down rare variant static summary
    fprintf(logout, "%s ", test.c_str());
    for(UINT i=0; i<power.size(); i++){
        fprintf(logout, "%-6.4f ", power[i]);
    }
    fprintf(logout, "| %s | %d\n", cmdcurrent.c_str(), unAnalyzedGens);
    funlockfile(logout);

    fclose(pValout);
    fclose(logout);
    fclose(statout);
    return 0;
}

