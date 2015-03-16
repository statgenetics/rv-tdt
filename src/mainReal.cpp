#include <cstdlib>
#include <numeric>
#include <algorithm>
#include <cmath>
#include <ctype.h>
#include <sys/types.h>
#include <dirent.h>
#include <stdlib.h>
#include <stdio.h>

#include "mainReal.h"
#include "tped.h"
#include "rareTdt.h"
#include "gw_utilities.h"
#include "gsl/gsl_cdf.h"
#include "Argument_helper.h"

int main(int argc, char* argv[])
{
	// input files
	std::string task = "phased";
	std::string projectName("rareTdtTest");
	std::string genoFile("NoInputFile"), phenoFile("NoInputFile"), mapFile("NoInputFile");
	// test
	double mafLower(0), mafUpper(0.05);
	bool nopermut(false);
	//bool biasTest(true);


	double alpha(0.00001);
	UINT nPermutations(200000), adaptive(500);
	//bool useSamMaf(false);
	//bool samplemaf(false);
	// handle missing
	double maxMissRatio(1);
	bool skipMiss(true);
	UINT minSiteForGene(4);
	std::string test_mode("A");

	/* UINT mafindex(3); */
	/* UINT weightindex(0); // =0 means use maf in untransmitted haplotypes has weights */

	dsr::Argument_helper arg;
	// required arguments
	/* arg.new_string("task", "\n\tAnalysis task. \n\t| OPTION: phased; unphased\n\t* phased: tdt test for rare variants on phased data \n\t* unphased: tdt test for rare variants on unphased data \n\t", task); */
	arg.new_string("projName", "\n\tProject name. (STRING) \n\t* set output folder name \n", projectName);
	// input files
	arg.new_named_string('G', "genoFile", "<string>", "\n\tGenotype file, with folder path. \n\t| Format (no header): SNP_ID 0 1 1 0 .... \n\t* missing site marked as -9 \n\t* If phased, nearby two genotypes stand for one allele\n", genoFile);
	arg.new_named_string('P', "phenoFile", "<string>", "\n\tPhenotype file, with folder path. \n\t| Format (no header): Ind_ID Fam_ID Father_ID Mother_ID Sex Status ... \n\t* Fahter_ID or Mother_ID = 0, if not available \n\t* SEX = (male)?1:0 \n\t* Status = (affected)?1:0 \n", phenoFile);
	arg.new_named_string('M', "mapFile", "<string>", "\n\tmapping file, with folder path. \n\t| Format (no header): Gene SNP_ID Annotation \n", mapFile);

	
	// test arguments
//	arg.new_named_string('t', "test_mode", "<string>", "\n\tP: paternal transmission only; M: maternal transmission only; A: all transmission \n", test_mode);
	/* arg.new_named_string('t', "test", "<string>", "\n\trvTDT test method. All test \n\t| OPTION: NoPermut; GenoPermut; HaploPermut; All \n\t* NoPermut: MZ, CMC^ \n\t* GenoPermut: VT-MZ-Geno, VT-CMC-Geno^, VT-WSS-Geno  \n\t* HapoPermut: VT-MZ-Haplo, VT-CMC-Halpo^, VT-WSS-Haplo  \n\t* All: all above \n\t* (^ for phased only) \n\t", test); */
	arg.new_named_double('l', "lower_cutoff", "<float>", "\n\tThe lower bound of variants to be included in analysis \n\t* Will compared with the third column of map file, only loci having annotated value > $lower_cutoff will be analyzed \n\t", mafLower);
	arg.new_named_double('u', "upper_cutoff", "<float>", "\n\tThe upper bound of variants to be included in analysis \n\t* Will compared with the third column of map file, only loci having annotated value <= $upper_cutoff will be analyzed \n\t", mafUpper);
	/* arg.new_named_double('i', "samMafCutoff", "<float>", "\n\tThe upper maf bound of sample MAF for variants to be included in analysis \n\t* Only loci having sample MAF < $samMafCutoff will be analyzed \n\t", samMafCutoff); */

	arg.new_named_unsigned_int('p', "permut", "<int>", "\n\tNumber of permutations \n\t* only applicable to permutation based methods\n\t", nPermutations);
	arg.new_named_unsigned_int('a', "adapt", "<int>", "\n\tNumber of permutations for adaptive check \n\t* only applicable to permutation based methods\n\t", adaptive);
	arg.new_named_double('s', "alpha", "<float>", "\n\tAlpha level for adaptive check. \n\t", alpha);

	/* arg.new_named_unsigned_int('j', "weight", "<int>", "\n\tColumn index of weights in mapfile (for WSS)\n\t", weightindex); */

	// handle missing
	arg.new_named_double('m', "maxMissRatio", "<float>", "\n\tThe max missing ratio allowed on each site \n\t* A variant site will be excluded from analysis, if its missing ratio > $maxMissRatio \n\t", maxMissRatio);
	/* arg.new_flag('k', "skipMiss", "\n\tSkip missing trios \n\t* if envoked, the incomplete trios will be inferred in analysis\n\t", skipMiss); */
	arg.new_named_unsigned_int('e', "minVariants", "<int>", "\n\tThe minimum variant sites for a gene to be analyzed \n\t* genes with variant site No. < $minVariants will be excluded from analysis (after check missing).\n\t", minSiteForGene);

	// 
//	arg.new_flag('b', "parental_bias", "\n\tTest for the parental transmission bias \n\t", biasTest);
	arg.new_flag('n', "nopermut", "\n\tOnly run CMC-analytical test\n\t", nopermut);



	// program information
	std::string program_name = "rvTDT";
	double version = 2.0;
	std::string banner = "\n\t:------------------------------------------------------:\n\t:  Transmission Disequilibrium Test for Rare Variants  :\n\t:------------------------------------------------------:\n\t: (c) 2015 Zong-Xiao He | http://bcm.edu/genetics/leal :\n\t:------------------------------------------------------:\n";

	arg.set_name(program_name.c_str());
	arg.set_description(banner.c_str());
	arg.set_version(version);
	arg.set_author("Zong-Xiao He <zongxiah@bcm.edu> & Suzanne M. Leal");
	time_t rawtime;
	struct tm * timeinfo;
	time ( &rawtime );
	timeinfo = localtime ( &rawtime );
	arg.set_build_date(asctime (timeinfo));

	arg.process(argc, argv);
	// Arguments check
	std::string pdirname(projectName + "_pval"), ldirname(projectName + "_rvTDT");
	if((!makedir(projectName + "_pval")) || (!makedir(projectName + "_rvTDT"))){
		std::cerr << "Error: can't create folder for pval and rvTDT"<< std::endl; exit(-1);
	}

	vector3F genotypes;
	vectorSS snps;
	vectorS geneNames;
	vector2F popMafs;
	bool phased=(task=="unphased")?false:true;
	/* bool hasPopMaf=(mapFile=="NoInputFile")?false:true; */

	readDataForRareTdt(genotypes, snps, geneNames, popMafs, phased, genoFile, phenoFile, mapFile); //


	for(UINT i=0; i<genotypes.size(); i++){
		vectorS tests;
		std::map<std::string, double> pvalues;
		Json::Value root;
		std::string logname = "./" + ldirname + "/" + geneNames[i] + ".rvTDT";
		std::ofstream log;
		log.open(logname.c_str());

		vectorF popmaf = popMafs[i];

	
		/* if(hasPopMaf){ */
		rareTdt tdtObj(phased, genotypes[i], popmaf, skipMiss, maxMissRatio, minSiteForGene);

		if(tdtObj.checkInformGene()) {

//			if(biasTest)
//				tdtObj.ParentalBiasTest(mafLower, mafUpper, pvalues, tests);

			if(nopermut)
				tdtObj.noPermut(mafLower, mafUpper, pvalues, tests);
			else {
				tdtObj.tdtTest(mafLower, mafUpper, adaptive, nPermutations, alpha, pvalues, tests, test_mode);
			}
			root = tdtObj.getJsonLog();
		}
		else 
			continue;
		/* }  */
		/* else  */
		/* { */
		/* 	rareTdt tdtObj(phased, genotypes[i], skipMiss, maxMissRatio, minSiteForGene); */
		/* 	if(tdtObj.checkInformGene()) { */

		/* 		if(biasTest) */
		/* 			tdtObj.ParentalBiasTest(mafLower, mafUpper, pvalues, tests); */

		/* 		if(nopermut) */
		/* 			tdtObj.noPermut(mafLower, mafUpper, pvalues, tests); */
		/* 		else */
		/* 			tdtObj.tdtTest(mafLower, mafUpper, adaptive, nPermutations, alpha, pvalues, tests); */
		/* 		root = tdtObj.getJsonLog(); */
		/* 	} */
		/* 	else continue; */
		/* } */

		root["basicInformation"]["variants"] = buildJsonArray(snps[i]);
		// update json value with returned pvalues
		Json::Value testRoot, tdtRes;
		for(UINT j=0; j<tests.size(); j++){
			if(pvalues[tests[j]] == std::numeric_limits<double>::infinity()) 
				tdtRes[tests[j]] = "-9";
			else 
				tdtRes[tests[j]] = pvalues[tests[j]];
		}

		testRoot["tests"] = buildJsonArray(tests);
		testRoot["pvalues"] = tdtRes;		
		root["tdtTest"] = testRoot;
		Json::StyledStreamWriter writer;
		writer.write(log, root);
		log.close();

		std::string pValName = "./" + pdirname + "/" + geneNames[i] + ".pval";
		FILE * pValout;
    		pValout = fopen(pValName.c_str(), "w");
		// first line
		fprintf(pValout, "#%s\t", "gene");
		for(UINT j=0; j<tests.size(); j++) {
			fprintf(pValout, "%s\t", tests[j].c_str());
		}
		fprintf(pValout, "\n");
		fprintf(pValout, "%s\t", geneNames[i].c_str());

		for(UINT j=0; j<tests.size(); j++) {
			if(pvalues[tests[j]] == std::numeric_limits<double>::infinity()) 
				fprintf(pValout, "%s\t", "null");
			else 
				fprintf(pValout, "%f\t",pvalues[tests[j]]);

		}
		fprintf(pValout, "\n");
		fclose(pValout);		
	}
	return 0;
}

/* vectorF recalculate_maf(vector2F& genos, vectorF& popmaf) */
/* { */

/* 	UINT nTrios = genos.size()/6; */
/* 	UINT nSites = genos[0].size(); */
/* 	vectorF missing(nSites, 0); */
/* 	vectorF allele(nSites, 0); */

/* 	for (UINT i=0; i<genos.size(); i++){ */
/* 		for (UINT j=0; j<genos[i].size(); j++){ */
/* 			if(genos[i][j] < 0)  */
/* 				missing[j] += 1; */

/* 			if(i%6 < 4) { */
/* 				if(genos[i][j] > 0)  */
/* 					allele[j] += 1; */
/* 			} */
/* 		} */
/* 	} */

/* 	vectorF maf(nSites, 0); */
/* 	for (UINT j=0; j<nSites; j++){ */
/* 		missing[j] /= nTrios*6; */
/* 		allele[j] /=nTrios*4; */

/* 		if(missing[j] < 0.1)  */
/* 			maf[j] = allele[j]; */
/* 		else */
/* 			maf[j] = popmaf[j]; */
/* 	} */

/* 	// use sample maf if missing ratio < 10%, otherwise use population maf */
/* 	return maf; */
/* } */


void readDataForRareTdt(vector3F& genos, vectorSS& snps, vectorS& geneNames, vector2F& popMafs, bool phased, std::string genoFile, std::string phenoFile, std::string mapFile)
{
	pedFileReader obj(genoFile, phased);
	vectorSS data;
	/* if(hasPopMaf && fileExists(mapFile))  */
	obj.loadMapFile(mapFile);
	obj.loadPhenoFile(phenoFile);

	// get genotypes and individual order, then reorganize the order
	vector3F noOrderGenos;
	geneNames.resize(0);
	popMafs.resize(0);
	
	vectorS indOrder; mapSvS phenoMap;
	obj.getGenotypes(noOrderGenos,snps, geneNames, popMafs);
	obj.getPhenotype(indOrder, phenoMap);
	// For RareTdt: re-order the genotypes, based on phenotype information; by order: father, mother, kid
	vectorSS families;

	for(UINT i=0; i<indOrder.size(); i++){
		std::string person = indOrder[i];
		vectorS pheno = phenoMap[person]; // pheno: famID, fatherID, motherID, sex, disease
		if(pheno[1].compare("0") != 0 && pheno[2].compare("0") != 0){
			vectorS oneFam;
			oneFam.push_back(pheno[1]); oneFam.push_back(pheno[2]); oneFam.push_back(person);
			families.push_back(oneFam);
		}
	}


	UINT ph=(phased)?2:1;
	genos.resize(0);
	BOOST_FOREACH(vector2F oneGene, noOrderGenos){
		if(oneGene.size() != indOrder.size()*ph) {
			std::cout << "error: unmatched genotype and phenotype file, different sample size\n"; exit(-1);
		}
		std::map<std::string, vector2F> caseToGeno;
		for(UINT i=0; i<indOrder.size(); i++){
			vector2F genoForAcase(oneGene.begin()+ph*i, oneGene.begin()+ph*i+ph); // in case phased
			caseToGeno.insert(std::make_pair(indOrder[i], genoForAcase));
		}

		vector2F Genotypes;
		UINT geneLength = oneGene[0].size();
		BOOST_FOREACH(vectorS fam, families){
			// check if all individual in this family available
			bool thisFamOK(true);
			BOOST_FOREACH(std::string ind, fam){
				// check if exist genoetypes for this individual
				if(caseToGeno.find(ind) == caseToGeno.end()) {thisFamOK = false; break;} 
				BOOST_FOREACH(vectorF oneChr, caseToGeno[ind]){
					if(oneChr.size() != geneLength) {thisFamOK = false; break;}
				} // check if the genotype has right length
			}

			if(thisFamOK){
				BOOST_FOREACH(std::string ind, fam){
					BOOST_FOREACH(vectorF oneChr, caseToGeno[ind]){
						Genotypes.push_back(oneChr);
					}
				}
			}
		}
		genos.push_back(Genotypes);
	}
	return;
}

