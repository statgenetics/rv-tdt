#include "gw_utilities.h"
#include "person.h"
#include "assocsims.h"
#include "Argument_helper.h"
#include "mainPower.h"
#include "tdtTrio.h"
#include "rareTdt.h"
#include "fb-skat/fb_skat.h"
#include <gsl/gsl_statistics.h>

using namespace std;

inline float skatPval(float v1, float v2){
	// y1=1-pchisq(x[,4],x[,5]);
	float res; 
	try {
		res = (v1 != v1 || v2 != v2)? NAN: (1 - gsl_cdf_chisq_P(v1, v2));
	} catch (char *e) {
		res = NAN;
	}
	return res;
}


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
	string test = "All";
	double mafLower = 0.0;
	double mafUpper = 0.05; //usually 5% cutoff for WSS and VT
	double alpha = 0.05;
	unsigned nPermutations = 2000;
	unsigned adaptive = 1000;
	unsigned nReplicates = 1000;
	// to simulate sequencing missing
	double indMiss(0), siteMiss(0);
	// analyze missing
	//double missingratio(0);
	bool skipMiss(true);
	bool usePOPmaf(true);

	// write information about every replicate
	bool recordGenos(false);
	bool recordInfos(false);
	bool recordStatic(false);

    // LD T1E
    double shuffleLD(0);

	// argsparser
	dsr::Argument_helper ah;
	// basic argument
	ah.new_string("task", "\n\tAnalysis task. \n\t| OPTION: phased; unphased; ld\n\t* phased: simulated phased data \n\t* unphased: simulated unphased data \n\t* ld-phased/unphased: fake LD data\n\t", simulationTask);
	ah.new_string("projName", "\n\tProject name. (STRING) \n\t* set output file names. \n\t", projectName);
	ah.new_optional_string("gdata", "\n\tPrefix of genetic data files. (STRING) \n\t* gdata.maf, gdata.sel and gdata.pos files are needed \n\t", gFile);
	// simulated genetic data
	ah.new_named_double('A', "rareCutOff", "<float>", "\n\tDefinition of rare variants. \n\t* variants with MAF<=$rareCutOff will be considered as \"rare\" variants\n\t", boundary);
	ah.new_named_double('B', "propFuncDel", "<float>", "\n\tProportion of FUNCTIONAL deleterious variants. \n\t* (1-proportion)x100\% is the proportion of non-causal deleterious variants (noise) \n\t* This is NOT the proportion of deleterious variants, which should have been defined in <gdata.sel> file\n\t", propFunctionalRv[0]);
	ah.new_named_double('C', "propFuncProtect", "<float>", "\n\tProportion of FUNCTIONAL protective variants. \n\t* (1-proportion)x100\% is the proportion of non-causal protective variants (noise) \n\t* This is NOT the proportion of protective variants, which should have been defined in <gdata.sel> file\n\t", propFunctionalRv[1]);
	ah.new_named_char('D', "inheritMode", "<char>", "\n\tMode of inheritance under which the phenotype data is simulated.\n\t| OPTION: A; D; R; M \n\t* \"A\" (additive), \"D\" (dominant), \"R\" (recessive), \"M\" (multiplicative), \n\t* \"C\" (compound dominant for non-mendelian traits, or compound recessive for mendelian traits)\n\t", moi);
	ah.new_named_double('E', "minDelOR", "<float>", "\n\t Minimum odds ratio for deleterious variants. \n\t* 0.0 for fixed effect size model; >=1.0 for variable effect sizes model\n\t", oddsRatios[0]);
	ah.new_named_double('F', "maxDelOR", "<float>", "\n\t Maximum odds ratio for deleterious variants. \n\t* >=1.0 AND > $minDelOR\n\t* will be the odds ratio for deleterious variants in fixed effect size model, or the maximum odds ratio in variable effect sizes model\n\t", oddsRatios[1]);
	ah.new_named_double('G', "minProtOR", "<float>", "\n\t Minimum odds ratio for protective variants. \n\t* 0.0 for fixed effect size model; <1.0 for variable effect sizes model\n\t", oddsRatios[2]);
	ah.new_named_double('H', "maxProtOF", "<float>", "\n\t Maximum odds ratio for protective variants. \n\t* <=1.0 AND > $minProtOR\n\t* will be the odds ratio for protective variants in fixed effect size model, or the maximum odds ratio in variable effect sizes model\n\t", oddsRatios[3]);
	ah.new_named_double('I', "comOR", "<float>", "\n\t Odds ratio for common variants. \n\t* =1.0 for neutral, >1.0 for deleterious, <1.0 for protective \n\t* this is the odds ratio for all variants having MAF > $rareCutOff, i.e., the \"common\" variants \n\t* assuming all common variants have fixed effect size\n\t", oddsRatios[4]);
	ah.new_named_double('J', "prevalence", "<float>", "\n\t Disease prevalence. \n\t* will be used as baseline penetrance of a gene (baseline penetrance ~= disease prevalence)\n\t", baselinef);
	ah.new_flag('Q', "keepSynonymous", "\n\tKeep synonymous variants from analysis. \n\t", isSynoKept);
	ah.new_flag('R', "removeComSites", "\n\tRemove common variant sites from analysis. \n\t* the \"common\" loci refers to variants in the haplotype pool having MAF > $rareCutOff \n\t", isCvTrimmed);
	ah.new_named_unsigned_int('S', "numCase", "<int>", "\n\t Number of cases (trios). \n\t", nCases);
	//ah.new_flag('a', "use_haplotype_pool", "\n\t[task=1,2,4,5] Randomly sample haplotypes from haplotype pool file $gdata.hap, rather than generating haplotypes on the fly. \n\t| BOOLEAN \n\t", shouldUseGenPool);
	ah.new_named_double('T', "neutralCutOff", "<float>", "\n\t Annotation value cut-off that defines a variant to be \"neutral\" in evolution (either synonymous or non-coding). \n\t* loci with annotation value (in gdata.ann) between (-$neutralCutOff, +$neutralCutOff) will be regarded \"neutral\" and will not contribute to phenotype \n\t", neutral_cutoff);
	ah.new_named_double('U', "indMiss", "<float>", "\n\t The proportion of individual that has entire missing, to mimic sample collecting process. \n\t| FLOAT\n\t", indMiss);
	ah.new_named_double('V', "siteMiss", "<float>", "\n\t The proportion of missing on each site, to mimic sequencing process. \n\t", siteMiss);

	// rvTDT arugment
	ah.new_flag('a', "skipMiss", "\n\t Skip missing trios \n\t* if envoked, the incomplete trios will be excluded from analysis\n\t", skipMiss);
	ah.new_flag('b', "useSamMaf", "\n\t Use sample maf to inference missing and maf threshold \n\t* will be envoked automatically if population maf aren't provided \n\t", usePOPmaf);
	//ah.new_named_double('c', "maxMissRatio", "<float>", "\n\t The max missing ratio allowed on each site \n\t* A variant site will be excluded from analysis, if its missing ratio > $maxMissRatio \n\t", missingratio);
	//
	ah.new_named_double('d', "lowerMaf", "<float>", "\n\t Lower bound of observed sample minor allele frequency \n\t* loci having observed MAF < $lowerMaf will not be analyzed \n\t", mafLower);
	ah.new_named_double('e', "upperMaf", "<float>", "\n\t Upper bound of observed sample minor allele frequency \n\t* loci having observed MAF > $upperMaf will not be analyzed \n\t", mafUpper);
	ah.new_named_string('t', "test", "<string>", "\n\t rvTDT test method. \n\t| OPTION: NoPermut; GenoPermut; HapoPermut; All \n\t* NoPermut: MZ, CMC^ \n\t* GenoPermut: VT-MZ-Geno, VT-CMC-Geno^, VT-WSS-Geno  \n\t* HapoPermut: VT-MZ-Hapo, VT-CMC-Hapo^, VT-WSS-Hapo  \n\t* All: all above \n\t* (^ for phased only) \n\t", test);
	ah.new_named_double('f', "alpha", "<float>", "\n\t Significance level at which power will be evaluated. \n\t", alpha);
	ah.new_named_unsigned_int('g', "replicates", "<int>", "\n\t Number of replicates for power evaluation\n\t", nReplicates);
	ah.new_named_unsigned_int('p', "permutations", "<int>", "\n\t Number of permutations, only applicable to permutation based methods \n\t", nPermutations);
	ah.new_named_unsigned_int('i', "adaptive", "<int>", "\n\t Number of permutations for adaptive check, only applicable to permutation based methods.\n\t", adaptive);
	ah.new_named_unsigned_int('j', "rng_seed", "<int>", "\n\t Seed for random number generator. \n\t* =0 is to use a random seed (seed = system time + process ID)\n\t", seed);
    ah.new_named_double('l', "shuffleLD", "<double>", "\n\t probability to shuffle fixed LD pattern. \n\t* =0 no shuffle; perfect LD\n\t", shuffleLD);


	// SKAT parameter
	UINT permut_max = 10000;
	double mend_th = 0.01;
	UINT min_site = 2;
	ah.new_named_unsigned_int('m', "permut_max", "<int>", "\n\t max permut for SKAT", permut_max);
    ah.new_named_double('n', "mend_th", "<double>", "\n\t mendilian theshold\n\t", mend_th);
	ah.new_named_unsigned_int('o', "min_size", "<int>", "\n\t min variant site", min_site);


	// output argument
	ah.new_flag('X', "recordGenos", "\n\t record genotypes in each replicates.\n\t* if envoked, write down genotypes in each replicates\n\t", recordGenos);
	ah.new_flag('Y', "recordInfos", "\n\t record genotypes in each replicates. \n\t* if envoked, write down information in each replicates\n\t", recordInfos);
	ah.new_flag('Z', "recordStatic", "\n\t record rare variant static summary in each replicates. \n\t* if envoked, write down information for each replicates\n\t", recordStatic);

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

    string haploFile = gFile + ".haplo";

	// Power calculations
	RNG rng;
	gsl_rng* gslr;
	if (seed == 0) {
		gslr = rng.get();
	}
	else {
		gslr = rng.get(seed);
	}

	vector2F mafDat, selDat;
	vector2UI posDat;
    vectorF ldMafs, ldProb;
    vector2F ldHaplo;
	
    bool phased = (simulationTask.find("phased")!=string::npos)?true:false;
    bool ld = (simulationTask.find("ld")!=string::npos)?true:false;
    //std::cout << phased << "|" << ld << std::endl;

    if(ld){
        // ld simulation
        // ld File: pobabilty haplo....
        // also, calculate mafs
        vector2F rawDat;
        if(!(scan_2dvector(haploFile, rawDat) )){
            std::cerr << "ERROR: Unable to open Haplotype data " << haploFile << ". Now Quit." << endl;
            exit(-1);
        }
        for(UINT i=0; i<rawDat.size(); i++){
            ldProb.push_back(rawDat[i][0]);
            vectorF oneHaplo;
            for(UINT j=1; j<rawDat[i].size(); j++) oneHaplo.push_back(rawDat[i][j]);
            ldHaplo.push_back(oneHaplo);
        }
        for(UINT i=0; i<ldHaplo[0].size(); i++){
            double one(0);
            for(UINT j=0; j<ldHaplo.size(); j++) one += ldProb[j]*ldHaplo[j][i]*0.1; // fake 0.1: BL01 and BL05
            ldMafs.push_back(one);
        }
        double sum(0);
        for(UINT i=0; i<ldProb.size(); i++) sum += ldProb[i];
        if(sum-1 != 0){
            std::cerr << "ERROR: sum of probability is not 1, can't use haplotype pool. Now Quit." << endl;
            exit(-1);
        }
        //std::cout << ldProb << std::endl << ldMafs << std::endl << "----\n";
        //for(UINT i=0; i<ldProb.size(); i++) std::cout << ldHaplo[i] << std::endl;
    } else {
        if(!(scan_2dvector(selFile, selDat) && scan_2dvector(mafFile, mafDat) && scan_2dvector(posFile, posDat))){
            std::cerr << "ERROR: Unable to open genetic data " << selFile << ". Now Quit." << endl;
            exit(-1);
        }
    }

	// output
	std::string pvalFile = projectName + ".pval";
	FILE * pValout;
	pValout = fopen(pvalFile.c_str(), "a");
	std::string StatFile = projectName + ".stat";
	std::string GenoDir = projectName + "_genos";

	//return 8 pvalues at most
	vectorUI iPowerful(20,0), iCount(20,0);
	UINT testNum(0), iReplicate(0);

	//record rare variant information//
	vector2F RaraVariantSummary;
    vector2F totPvals;

	while (iReplicate != nReplicates) {
        vectorF mafs, keepedMafs; 
        UINT dataIdx;
        vector2F genotypesTDT;
        if(ld){
            keepedMafs = ldMafs;
            mafs = ldMafs;
            dataIdx=iReplicate;
            for(UINT i=0; i<nCases; i++){
                UINT idx1=ldGenHaplo(ldProb, gslr); // if idx=0; wildtype
                UINT idx2=ldGenHaplo(ldProb, gslr);
                UINT idx3=ldGenHaplo(ldProb, gslr);
                UINT idx4=ldGenHaplo(ldProb, gslr);

                vectorF haplo1=ldHaplo[idx1]; // kids 1
                vectorF haplo2=ldHaplo[idx2]; 
                vectorF haplo3=ldHaplo[idx3]; // kids 2
                vectorF haplo4=ldHaplo[idx4]; 

                if(shuffleLD>0){
                    if(idx1) haplo1=breakLD(haplo1, shuffleLD, gslr);
                    if(idx2) haplo2=breakLD(haplo2, shuffleLD, gslr);
                    if(idx3) haplo3=breakLD(haplo3, shuffleLD, gslr);
                    if(idx4) haplo4=breakLD(haplo4, shuffleLD, gslr);
                }

                genotypesTDT.push_back(haplo1); genotypesTDT.push_back(haplo2); genotypesTDT.push_back(haplo3); genotypesTDT.push_back(haplo4); 
                genotypesTDT.push_back(haplo1); genotypesTDT.push_back(haplo3); // kids
            }
            //std::cout << genotypesTDT.size() << std::endl;
            //for(UINT i=0; i<genotypesTDT.size(); i++) std::cout << genotypesTDT[i] << std::endl;
        } else {
            /*** Choose MAF data ***/
            dataIdx = gsl_rng_uniform_int(gslr, mafDat.size());
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

            // remove all wildtype sites
            vector2F genotypeTDT = zx_genGenos(genotypes, keepedMafs, gslr);
            vectorL unInformSites(keepedMafs.size(), true);
            for(UINT i=0; i<genotypeTDT[0].size(); i++){
                for(UINT j=0; j<genotypeTDT.size(); j++){
                    // keep this site if any trio is informative
                    if(genotypeTDT[j][i]>0) {unInformSites[i]=false; break;}
                }
            }
            vectorF newMafs;
            for(UINT i=0; i<unInformSites.size(); i++) {if(!unInformSites[i]) newMafs.push_back(keepedMafs[i]);}
            for(UINT i=0; i<genotypeTDT.size(); i++) {
                vectorF oneHapo;
                for(UINT j=0; j<genotypeTDT[i].size(); j++) {if(!unInformSites[j]) oneHapo.push_back(genotypeTDT[i][j]);}
                genotypesTDT.push_back(oneHapo);
            }
            keepedMafs = newMafs;
            // simulate sequencing process, introduce missing
            // sequencing --> missing
            if(indMiss > 0 || siteMiss > 0) zx_createMissing(genotypesTDT, indMiss, siteMiss, gslr);
        }

        /*** rare variant summary ***/
        RaraVariantSummary.push_back(zx_vStatSum(dataIdx, boundary, mafs, keepedMafs, genotypesTDT));

		// minimum variant sites
		if(keepedMafs.size () < min_site) 
			continue;

		/*** TDT TESTs ***/
		std::string geneName = "Gene_" + n2s(iReplicate);
		std::string projName = projectName + "_TDT";
		if(!phased) genotypesTDT = zx_phased2unphased(genotypesTDT);

		rareTdt tdtObj(phased, genotypesTDT, keepedMafs, skipMiss, 1, 0);
		if(!usePOPmaf) tdtObj = rareTdt(phased, genotypesTDT, skipMiss, 1, 0);
		std::vector<std::string> tests;
		std::map<std::string, double> pvalues;
		tdtObj.tdtTest(test, mafLower, mafUpper, adaptive, nPermutations, 0.05, pvalues, tests);


		/*** SKAT test ***/
		vector2F genos = zx_phased2unphased(genotypesTDT);
		UINT num_sites = genos[0].size();
		UINT num_trios = genos.size()/3;
		vectorF Y(num_trios, 1), pass(num_sites, 1), weight(num_sites, 1);
			
		/* SKAT use sample maf as cutoff */
		for(UINT i=0; i<num_sites; i++){
			if( keepedMafs[i] > mafUpper) 
				pass[i] = 0;
		}

		std::string gene_name = "Rep_" + n2s(iReplicate);
		fbSkat obj(gene_name, genos, Y, pass, weight);
		obj.TestLoad(permut_max, 1, mend_th, 0, min_site); // skat
		vectorF vals = obj.SKAT();
		tests.push_back( "skat-1");
	    pvalues["skat-1"] = skatPval(vals[0], vals[1]);
		tests.push_back( "skat-2");
	    pvalues["skat-2"] = skatPval(vals[2], vals[3]);

		obj.TestLoad(permut_max, 1, mend_th, 1, min_site); // skat-burden
		vals = obj.SKAT();
		tests.push_back( "skat-burden-1");
	    pvalues["skat-burden-1"] = skatPval(vals[0], vals[1]);
		tests.push_back( "skat-burden-2");
	    pvalues["skat-burden-2"] = skatPval(vals[2], vals[3]);


		if(recordInfos) {
			Json::Value JsonLog=tdtObj.getJsonLog();
			Json::Value basic;
			basic["task"] = projName;
			basic["gene"] = geneName;
			basic["phased"] = Json::Value(phased);
			JsonLog["basic"] = basic;
			Json::Value tdtRes;
			for(UINT i=0; i<tests.size(); i++){
				if(pvalues[tests[i]] == std::numeric_limits<double>::infinity()) tdtRes[pvalues[tests[i]]] = "null";
				else tdtRes[pvalues[tests[i]]] = pvalues[tests[i]];
			}
			tdtRes["tests"] = buildJsonArray(tests);
			JsonLog["tdtResult"] = tdtRes;

			if(makedir(projName + "_log")){
				std::string logname = "./" + projName + "_log/" + geneName + ".log";
				std::ofstream log;
				log.open(logname.c_str());
				Json::StyledStreamWriter writer;
				writer.write(log, JsonLog);
				log.close();
			}
		}
	
		if(!testNum) testNum = tests.size();
		vectorF pvals;
		for(UINT i=0; i<tests.size(); i++) 
			pvals.push_back(pvalues[tests[i]]);
        totPvals.push_back(pvals);
        //std::cout << pvals << "|" << tests << std::endl;
		for(UINT i=0; i<pvals.size(); i++) {
			if(pvals[i] != std::numeric_limits<double>::infinity() && pvals[i] == pvals[i]) {
				iCount[i]++; // one test done
				if(pvals[i] <= alpha) 
					iPowerful[i]++; // and it's significant
			}
		}

		//flockfile(pValout);
        //fprintf(pValout, "%d\t", dataIdx+1);
		//for(UINT i=0; i<pvals.size(); i++) {fprintf(pValout, "%f\t", pvals[i]);}
		//fprintf(pValout, "\n");
		//funlockfile(pValout);

		if(recordGenos){
			if(makedir(GenoDir)){
				std::string genoFile = "./" + GenoDir + "/geno" + n2s(iReplicate);
				FILE * Genout;
				Genout = fopen(genoFile.c_str(), "w");
				for(UINT i=0; i<genotypesTDT.size(); i++){
					for(UINT j=0; j<genotypesTDT[i].size(); j++) fprintf(Genout, "%d\t", int(genotypesTDT[i][j]));
					fprintf(Genout, "\n");
				}
				fclose(Genout);
			}
		}
		++iReplicate;
	}

    flockfile(pValout);
    for(UINT i=0; i<totPvals.size(); i++){ 
        for(UINT j=0; j<totPvals[i].size(); j++){
            fprintf(pValout, "%f\t", totPvals[i][j]);
        }
        fprintf(pValout, "\n");
    }
    funlockfile(pValout);

	vectorF power;
	for(UINT i=0; i<testNum; i++) 
		power.push_back(1.0*iPowerful[i]/(1.0*iCount[i]));

	// get mean for each entry (column)
	vectorF RaraVariantStatic(RaraVariantSummary[0].size(), 0.0);
	for(UINT i=0; i<RaraVariantSummary.size(); i++){
		for(UINT j=0; j<RaraVariantSummary[i].size(); j++){
			RaraVariantStatic[j] += RaraVariantSummary[i][j];
		}
	}
	for(UINT i=0; i<RaraVariantStatic.size(); i++) RaraVariantStatic[i]=RaraVariantStatic[i]/(RaraVariantSummary.size()); 

	if(recordStatic){
		FILE * statout;
		statout = fopen(StatFile.c_str(), "a");
		flockfile(statout);
		for(UINT i=0; i<RaraVariantSummary.size(); i++){
			fprintf(statout, "%d ",int(RaraVariantSummary[i][0]));
			for(UINT j=1; j<RaraVariantSummary[i].size(); j++){
				fprintf(statout, "%-6.4f ", RaraVariantSummary[i][j]);
			}
			fprintf(statout, "\n");
		}
		funlockfile(statout);
		fclose(statout);
	}

	std::string logFile = projectName + ".log";
	FILE * logout;
	logout = fopen(logFile.c_str(), "a");
	flockfile(logout);

	fprintf(logout, "#projName gFile nCases boundary neutral_cutoff ORlow ORhigh prop mafUpper skipMiss usePOPmaf siteMiss oriVarSites AnalyzedSite  nVarNotInTrio  nVarInParent  nVarInKid  nComVariantInTrio  nTrioHasVarInParent  nTrioHasVarInKid  MissInParent  MissInSite test power1 power2 power3 (power4)\n");
	fprintf(logout, "%s %s %d %-4.2f %-7.6f %-3.1f %-3.1f %-4.2f %-4.2f %d %d %-4.2f ", projectName.c_str(), gFile.c_str(), nCases, boundary, neutral_cutoff, oddsRatios[0], oddsRatios[1], propFunctionalRv[0], mafUpper, skipMiss, usePOPmaf, siteMiss);
	// variant summary
	for(UINT i=1; i<RaraVariantStatic.size(); i++){
		fprintf(logout, "%-6.4f ", RaraVariantStatic[i]);
	} 
	fprintf(logout, "%s ", test.c_str());
	for(UINT i=0; i<power.size(); i++){
		fprintf(logout, "%-6.4f ", power[i]);
	}
	fprintf(logout, "| %s | ", cmdcurrent.c_str());
	for(UINT i=0; i<power.size(); i++){
		fprintf(logout, "%d ", iCount[i]);
	}
	fprintf(logout, "\n");
	funlockfile(logout);

	fclose(pValout);
	fclose(logout);
	return 0;
}
