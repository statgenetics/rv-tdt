/**
 * @file        rareTdt.cpp
 * @author      Zong-Xiao He (zongxiah@bcm.edu)
 * @copyright   Copyright 2014 Zong-Xiao He
 * @date        2014-06-14
 * @see		    rareTdt.hpp
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program. If not, see <http://www.gnu.org/licenses/>.
 *
 *
 * @brief
 *
 */


#include "gw_utilities.h"
#include "tdtTrio.h"
#include "rareTdt.h"
#include <math.h> 
//#include "unphasedTDT.h"

#define PR(a) std::cout << #a " = " << a << std::endl;

void rareTdt::m_preprocess()
{

	__test_mode = "A";
	__first_call = true;

	UINT vNum = __genotypes[0].size();
	__rowsPerTrio = (__phased)?6:3;
	if(__genotypes.size()%__rowsPerTrio != 0) {
		__informGene=false;
		return;
	}
	__nInformTrios = __genotypes.size()/__rowsPerTrio;
	__missingRatio.assign(vNum, 0);
	__samMafs.assign(vNum, 0);
	vectorUI parentNotMiss(vNum, 0); //need this to calculate sample maf
	double maxGenos = (__phased)?1:2;
	vectorL varInParent(vNum, false), varInKid(vNum, false);

	for (UINT i = 0; i < __genotypes.size(); i++){
		if(__genotypes[i].size() != vNum) {
			__informGene=false;
			return;
		}
		for (UINT j=0; j<__genotypes[i].size(); j++){
			if (__genotypes[i][j] < 0 || __genotypes[i][j] > maxGenos) {
				__missingRatio[j]++;
			} else {
				if(i%__rowsPerTrio <= (__rowsPerTrio*2)/3-1) {
					__samMafs[j] += __genotypes[i][j]; //only consider parents when calculating sample maf
					parentNotMiss[j]++;
					varInParent[j] = true;
				} else {
					varInKid[j] = true;
				}
			} // record missing and sample maf
		}
	}

	for(UINT i=0; i<vNum; i++){
		__samMafs[i] = (__samMafs[i]==0)?1:__samMafs[i]; // in case zero 
		__samMafs[i] = __samMafs[i]/float(parentNotMiss[i]*(6/__rowsPerTrio));
		__missingRatio[i] = __missingRatio[i]/ float(__genotypes.size());
	}

	// flip genotype if samMaf > 0.5
	bool flip(false);
	vectorL needFlip(vNum, false);
	for(UINT i=0; i<vNum; i++){
		if(__samMafs[i]>0.5) {
			flip=true; 
			needFlip[i]=true;
		}
	}
	if(flip){
		for(UINT i=0; i<vNum; i++){
			if(needFlip[i]){
				for(UINT x=0; x<__genotypes.size(); x++) {
					if(__genotypes[x][i] >= 0 && __genotypes[x][i] <= maxGenos) __genotypes[x][i]=maxGenos-__genotypes[x][i];
				}
			}
		}
	}

	// which variant sites to be analyzed
	__varToBeAnalyzed.assign(vNum, true);
	if(__missCutoff<1 && __missCutoff>0){
		for(UINT i = 0; i < __varToBeAnalyzed.size(); i++) {
			if(__missingRatio[i] > __missCutoff) 
				__varToBeAnalyzed[i]=false;
		}
	}


	if(!__hasPopMaf) __popMafs=__samMafs;
	// write log Json
	Json::Value vStatic;
	vStatic["varAnnotations"] = buildJsonArray(__popMafs);
	vStatic["sampleMafs"] = buildJsonArray(__samMafs);
	vStatic["missingRatio"] = buildJsonArray(__missingRatio);
	vStatic["varFoundInParent"] = buildJsonArray(varInParent);
	vStatic["varFoundInKid"] = buildJsonArray(varInKid);
	logRoot["variantStatic"] = vStatic;

	return;
}


rareTdt::rareTdt(bool phased, vector2F& genotypes, const vectorF& popMafs, bool skipMiss, double missCutoff, UINT minVarNum): __phased(phased), __genotypes(genotypes), __popMafs(popMafs), __hasPopMaf(true), __skipMiss(skipMiss), __missCutoff(missCutoff), __informGene(true) {
	if(__popMafs.size() != __genotypes[0].size()) {
		__informGene=false;
		return;
	}
	m_preprocess();
	if(m_getInformVarNum()<minVarNum){
		__informGene=false;
		return;
	}
	m_load();
	m_tdtTableLoad();
	return;
}

rareTdt::rareTdt(bool phased, vector2F& genotypes, bool skipMiss, double missCutoff, UINT minVarNum): __phased(phased), __genotypes(genotypes), __hasPopMaf(false), __skipMiss(skipMiss), __missCutoff(missCutoff), __informGene(true) {
	m_preprocess();
	if(m_getInformVarNum()<minVarNum){
		__informGene=false;
		return;
	}
	m_load();
	m_tdtTableLoad();
	return;
}

UINT rareTdt::m_getInformVarNum() const {
	UINT varNum(0);
	BOOST_FOREACH(bool yes, __varToBeAnalyzed) {
		if(yes) varNum++;
	}
	return varNum;
}

void rareTdt::m_load(){
	for(UINT i = 0; i < __genotypes.size()/__rowsPerTrio; i++){
		vector2F OneFamily; OneFamily.resize(__rowsPerTrio);
		copy(__genotypes.begin() + i*__rowsPerTrio, __genotypes.begin() + (i+1)*__rowsPerTrio, OneFamily.begin());

		if(__phased){
			phasedTrioPtr oneTrio = phasedTrioPtr(new phasedTrio);
			oneTrio->load(OneFamily, __popMafs, __varToBeAnalyzed, __skipMiss);
			if(oneTrio->checkInformTrio()){
				trioPtr baseTrio = oneTrio;
				__trios.push_back(baseTrio);
			}   
			else __nInformTrios--;
		} else {
			unPhasedTrioPtr oneTrio = unPhasedTrioPtr(new unPhasedTrio);
			oneTrio->load(OneFamily, __popMafs, __varToBeAnalyzed, __skipMiss);
			if(oneTrio->checkInformTrio()){
				trioPtr baseTrio = oneTrio;
				__trios.push_back(baseTrio);
			}
			else __nInformTrios--;
		}
	}
	if(__trios.size()==0) __informGene=false;
	return;
}

void rareTdt::m_tdtTableLoad(){
	if(!__informGene) return;
	__denovoCount.resize(0); __effTriosBySite.resize(0); 
	vector2F unTrans;
	for(trioIterator iter=__trios.begin(); iter!=__trios.end(); iter++){
		if(!(*iter)->checkInformTrio()) continue;
		vector2F tdtTable = (*iter)->getTdtTable();
		for(UINT i=0; i<tdtTable.size(); i++){
			if(i%2==0) __tdtBtable.push_back(tdtTable[i]);
			else __tdtCtable.push_back(tdtTable[i]);
		}
		vector2F unTransObj = (*iter)->getNuTransPatChrs();
		for(UINT i=0; i<unTransObj.size(); i++) unTrans.push_back(unTransObj[i]);
		m_ifTrueThenAddOne((*iter)->getDenovo(), __denovoCount);
		m_ifTrueThenAddOne((*iter)->getAnalyzed(), __effTriosBySite);
	}
	__wssWeight = m_wssWeight(unTrans); // need __effTriosBySite
	__external_weights = false;

    //for(UINT i=0; i<__tdtBtable.size(); i++) std::cout << __tdtBtable[i] << std::endl;
	// write log Json
	Json::Value tdtStat;
	tdtStat["effectiveTrioNum"] = __nInformTrios;
	tdtStat["siteBeAnalyzed"] = buildJsonArray(__varToBeAnalyzed);
	tdtStat["Transmitted"] = buildJsonArray(m_colTotal(__tdtCtable));
	tdtStat["unTransmitted"] = buildJsonArray(m_colTotal(__tdtBtable));
	tdtStat["denovoCount"] = buildJsonArray(__denovoCount);
	tdtStat["effectiveTrioNum"] = buildJsonArray(__effTriosBySite);
	tdtStat["wssWeight"] = buildJsonArray(__wssWeight);
	tdtStat["singleSitePval"] = buildJsonArray(m_singleSiteP());
	logRoot["tdtStatic"] = tdtStat;
	return;
}

vectorF rareTdt::m_singleSiteP() const {
	vectorF singleSiteP;
	vectorF Bcount(m_colTotal(__tdtBtable)), Ccount(m_colTotal(__tdtCtable));
	for(UINT i=0; i<Bcount.size(); i++){
		if((Ccount[i] + Bcount[i]) == 0) singleSiteP.push_back(std::numeric_limits<double>::infinity());
		else singleSiteP.push_back(gsl_cdf_gaussian_Q((Ccount[i]-Bcount[i])/sqrt(Ccount[i]+Bcount[i]), 1));
	}
	return singleSiteP;
}

vectorF rareTdt::m_wssWeight(const vector2F& unTrans)
{
	vectorUI NotZero(unTrans[0].size(),0), NotMiss(unTrans[0].size(),0);
	for(UINT i = 0; i < unTrans.size(); i++){
		for(UINT j = 0; j < unTrans[i].size(); j++){
			if(unTrans[i][j] >= 0.0) {NotMiss[j] ++; NotZero[j] += unTrans[i][j];}
		}
	}
	//std::cout << NotMiss << std::endl << NotZero << std::endl;
	vectorF weight(unTrans[0].size(),0);
	// w=sqrt(n*q*(1-q))
	for(UINT i = 0; i < NotMiss.size(); i++){
		if(NotMiss[i] == 0) weight[i]=0;
		else {
			double q = (NotZero[i]+1)/float(NotMiss[i]+1);
			if(q==1) q = (NotZero[i]+1)/float(NotMiss[i]+2);
			else if(q>1) {weight[i]=0; continue;}
			else {;}
			weight[i] = 1/sqrt(__effTriosBySite[i]*4*q*(1-q)); // only consider parents
		}
	}
	return weight;
}

/* set up for rvTDT test */
// shuffle trios, get tdt table and untransmitted haplotypes 
rvTDTshuffleOutput rareTdt::m_shuffle(std::string shuffleMethod, gsl_rng* gslr)
{
	vector2F unTrans, tdtBtable, tdtCtable;
	for(trioIterator iter=__trios.begin(); iter!=__trios.end(); iter++){
		if(!(*iter)->checkInformTrio()) continue;
		shuffleOutput shuffled = (*iter)->shuffle(shuffleMethod,gslr);

		if(shuffled.tdtTable.size()%2 != 0 ) {std::cerr << "Error: shuffle return a invalid structure" << std::endl; exit(-1);}
		for(UINT i=0; i<shuffled.tdtTable.size(); i++){
			if(i%2 == 0) tdtBtable.push_back(shuffled.tdtTable[i]);
			else tdtCtable.push_back(shuffled.tdtTable[i]);
		}
		unTrans.push_back(shuffled.untransmitted[0]); unTrans.push_back(shuffled.untransmitted[1]);
	}
	return rvTDTshuffleOutput(unTrans, tdtBtable, tdtCtable);
}

// aggreate events count
// TODO: how to deal with missing infer in CMC 
tdtResult rareTdt::m_rvAggreate(const vector2F& tdtBtable, const vector2F& tdtCtable, const vectorL& tobeAnalyzed, const vectorF& weight, bool allowDenovo) 
{
    //for(UINT i=0; i<__tdtBtable.size(); i++) std::cout << __tdtBtable[i] << std::endl;
    //for(UINT i=0; i<__tdtCtable.size(); i++) std::cout << __tdtCtable[i] << std::endl;
    /* std::cout << tobeAnalyzed << " | " << weight << " | " << allowDenovo << std::endl; */
	tdtResult errorResult(0, 0, std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity());
	if(tdtBtable.size() != tdtCtable.size() || tdtBtable[0].size() != tobeAnalyzed.size() || tdtBtable[0].size() != weight.size()) return errorResult;

	double MZ_B(0), MZ_C(0), CMC_B(0), CMC_C(0);

	// quick and dirty way to check trans/untrans per trio
	unsigned int num_trio = __genotypes.size()/__rowsPerTrio;
	vectorF father_trans(num_trio, 0);
	vectorF father_untrans(num_trio, 0);
	vectorF mother_trans(num_trio, 0);
	vectorF mother_untrans(num_trio, 0);

	for(UINT i=0; i<tdtBtable.size(); i++){

		unsigned int fam_id = i/2;
		double trans(0), untrans(0);

		// parental test
		// father
		if(__test_mode == "M" && i%2 == 0) 
			continue;

		// mother
		if(__test_mode == "P" && i%2 == 1) 
			continue;

		if(tdtBtable[i].size() != tdtCtable[i].size()) return errorResult;
		double oneB(0), oneC(0), floorB(0), floorC(0);
		for(UINT j=0; j<tdtBtable[i].size(); j++){

			if(tobeAnalyzed[j]) {

				oneB += tdtBtable[i][j]*weight[j]; floorB += floor(tdtBtable[i][j])*weight[j];

				untrans = untrans + tdtBtable[i][j];
				// denovo marked as 101
				if(tdtCtable[i][j]>=100){
					double offset=(allowDenovo)?100:101;
					oneC += (tdtCtable[i][j]-offset)*weight[j]; 
					floorC += floor(tdtCtable[i][j]-offset)*weight[j];

					trans = trans + tdtCtable[i][j] - offset;
				}else {
					oneC += tdtCtable[i][j]*weight[j];
					floorC += floor(tdtCtable[i][j])*weight[j];

					trans = trans + tdtCtable[i][j];
				}

			}
		}

		// father
		if(__first_call) {
			if(i%2 == 0) {
				father_trans[fam_id] = trans;
				father_untrans[fam_id] = untrans;
			} else {
				mother_trans[fam_id] = trans;
				mother_untrans[fam_id] = untrans;
			}
		}

		MZ_B += oneB; MZ_C += oneC;
		if(__skipMiss && (oneB + oneC !=0)) {CMC_B += oneB/(oneB+oneC); CMC_C += oneC/(oneB+oneC);} // skip miss
		if((!__skipMiss) && (floorB + floorC !=0)) {CMC_B += floorB/(floorB+floorC); CMC_C += floorC/(floorB+floorC);} // infer miss
	}

	double MZ_chi2(0), CMC_chi2(0);
	double MZ_P(std::numeric_limits<double>::infinity()), CMC_P(std::numeric_limits<double>::infinity());
	if(MZ_B+MZ_C !=0 ){
		/* MZ_chi2 = (MZ_B - MZ_C)/sqrt(MZ_B+MZ_C); */
		MZ_chi2 = (MZ_C - MZ_B)/sqrt(MZ_B+MZ_C);
		MZ_P = gsl_cdf_gaussian_Q(MZ_chi2, 1);
	}
	if(CMC_B+CMC_C !=0 ){
		// check the overtransmission of wildtype in unaffected trios, 09/05/2014
		/* CMC_chi2 = (CMC_B - CMC_C)/sqrt(CMC_B+CMC_C); */
		CMC_chi2 = (CMC_C - CMC_B)/sqrt(CMC_B+CMC_C);
		CMC_P = gsl_cdf_gaussian_Q(CMC_chi2, 1);
	}

	// quick and dirty way to check trans/untrans per trio
	if(__first_call) {

		Json::Value transStat;
		transStat["father_trans"] = buildJsonArray(father_trans);
		transStat["father_untrans"] = buildJsonArray(father_untrans);
		transStat["mother_trans"] = buildJsonArray(mother_trans);
		transStat["mother_untrans"] = buildJsonArray(mother_untrans);
		logRoot["transStat"] = transStat;
	
		__first_call = false;
	}


	/* std::cout << __test_mode << std::endl; */
	/* std::cout << MZ_C << " " << MZ_B << " " << CMC_C << " " << CMC_B << std::endl; */
    /* std::cout << MZ_chi2 <<" " <<  CMC_chi2 <<" " <<  MZ_P <<" " <<  CMC_P << std::endl; */
	return tdtResult(MZ_chi2, CMC_chi2, MZ_P, CMC_P);
}

// variable threshold
tdtResult rareTdt::m_VT(const vector2F& tdtBtable, const vector2F& tdtCtable, const vectorL& tobeAnalyzed, const vectorF& weight) 
{
	tdtResult errorResult(0, 0, std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity());
	if(tdtBtable.size() != tdtCtable.size() ||  tdtBtable[0].size() != tobeAnalyzed.size() || tdtBtable[0].size() != weight.size() || tdtBtable[0].size() != __popMafs.size()) return errorResult;

	vectorF validMafs;
	for(UINT i=0; i<tobeAnalyzed.size(); i++){ if(tobeAnalyzed[i]) validMafs.push_back(__popMafs[i]);}
	sort(validMafs.begin(), validMafs.end());
	validMafs.erase(unique(validMafs.begin(), validMafs.end()), validMafs.end());
	
	double maxMZ(-999.0), maxCMC(-999.0);
	for(UINT i=0; i < validMafs.size(); i++){
		vectorL oneAnalyze(tobeAnalyzed);
		for(UINT j=0; j<tobeAnalyzed.size(); j++){ if(__popMafs[j]>validMafs[i]) oneAnalyze[j] = false;}
		tdtResult res = m_rvAggreate(tdtBtable, tdtCtable, oneAnalyze, weight, false);
		if (res.mz_chi2 > maxMZ) maxMZ = res.mz_chi2;
		if (res.cmc_chi2 > maxCMC) maxCMC = res.cmc_chi2;
	}
	return tdtResult(maxMZ, maxCMC, std::numeric_limits<double>::infinity(), std::numeric_limits<double>::infinity()); // will never use MZ_P, CMC_P
}

/* rvTDT test */ 
// return pvalue for CMC-anal
vectorF rareTdt::m_noPermut() 
{
	vectorF fakeWeight(__varToBeAnalyzed.size(), 1);

	tdtResult oriBurden = m_rvAggreate(__tdtBtable, __tdtCtable, __varToBeAnalyzed, fakeWeight, true); // MZ CMC 01
        
	if(oriBurden.mz_pval == std::numeric_limits<double>::infinity()){
		vectorF res(1, std::numeric_limits<double>::infinity());
		return res;
	}

	vectorF nopermut;
	nopermut.push_back(oriBurden.cmc_pval); //CMC
	return nopermut;
}


void rareTdt::noPermut(double lowerMaf, double upperMaf, std::map<std::string, double>& pvalues, std::vector<std::string>& tests)
{

	/* pvalues.clear(); tests.clear(); */

	for(UINT i=0; i<__popMafs.size(); i++) {
		if(__popMafs[i]<=lowerMaf || __popMafs[i]>upperMaf) 
			__varToBeAnalyzed[i] = false;
	}
	// update _varToBeAnalyzed
	logRoot["tdtStatic"]["siteBeAnalyzed"] = buildJsonArray(__varToBeAnalyzed);

    vectorF np = m_noPermut(); 

	pvalues["CMC-Analytical"] = np[0]; tests.push_back("CMC-Analytical");

	return;
}



void rareTdt::updateWeights(const vectorF& weights){
	if(weights.size() == __wssWeight.size()){
		__wssWeight = weights; // use the given weights
		__external_weights = true;
		logRoot["tdtStatic"]["wssWeight"] = buildJsonArray(__wssWeight);
	}
	return;
}



/* rvTDT test */ 
// only haplotype permutation is allowed: std::string shuffleMethod = "haplo..."
// return pvalue for BRV, CMC, VT-BRV, VT-CMC, WSS
vectorF rareTdt::m_Permut(UINT adaptive, UINT PermutateTimes, double alpha, gsl_rng* gslr) 
{
	vectorF fakeWeight(__varToBeAnalyzed.size(), 1);

    tdtResult oriBurden = m_rvAggreate(__tdtBtable, __tdtCtable, __varToBeAnalyzed, fakeWeight, false); // MZ CMC 01
    tdtResult oriVT = m_VT(__tdtBtable, __tdtCtable, __varToBeAnalyzed, fakeWeight); // VT-MZ VT-CMC
    tdtResult oriWSS = m_rvAggreate(__tdtBtable, __tdtCtable, __varToBeAnalyzed, __wssWeight, false); // WSS

    //for(UINT x=0; x<__tdtBtable.size(); x++) std::cout << "B:\t" << __tdtBtable[x] << std::endl << "C:\t" << __tdtCtable[x] << std::endl;
    //std::cout << "-------------\n";

    if(oriBurden.mz_pval == std::numeric_limits<double>::infinity() || oriWSS.mz_pval == std::numeric_limits<double>::infinity() || oriVT.mz_chi2 == -99 || oriVT.cmc_chi2 == -99){
        vectorF res(5, std::numeric_limits<double>::infinity());
        return res;
    }

    vectorF oriRes;
    oriRes.push_back(oriBurden.mz_chi2); 
	oriRes.push_back(oriBurden.cmc_chi2);  // MZ CMC 
    oriRes.push_back(oriVT.mz_chi2); 
	oriRes.push_back(oriVT.cmc_chi2); 
	oriRes.push_back(oriWSS.mz_chi2); // VT-MZ; VT-CMC; WSS
    /* std::cout << oriRes << std::endl; */
	/* std::cout << "===" << std::endl; */

    vectorUI permCount(5,0);
    vectorF permPval(5,9.0);

    for(UINT i=1; i<=PermutateTimes; i++)
    {
        rvTDTshuffleOutput Obj=m_shuffle("HapoShuffle", gslr);
        //for(UINT x=0; x<Obj.tdtBtable.size(); x++) std::cout << "B:\t" << Obj.tdtBtable[x] << std::endl << "C:\t" << Obj.tdtCtable[x] << std::endl; 
       
        tdtResult permBurden = m_rvAggreate(Obj.tdtBtable, Obj.tdtCtable, __varToBeAnalyzed, fakeWeight, false); // MZ CMC
        tdtResult permVT = m_VT(Obj.tdtBtable, Obj.tdtCtable, __varToBeAnalyzed, fakeWeight); // VT-MZ, VT-CMC

		vectorF newWeights = __external_weights? __wssWeight:m_wssWeight(Obj.unTrans);
        tdtResult permWSS = m_rvAggreate(Obj.tdtBtable, Obj.tdtCtable, __varToBeAnalyzed, newWeights, false);

        vectorF permRes; 
        permRes.push_back((permBurden.mz_pval==std::numeric_limits<double>::infinity())?0:permBurden.mz_chi2); //MZ01
        permRes.push_back((permBurden.cmc_pval==std::numeric_limits<double>::infinity())?0:permBurden.cmc_chi2); //CMC01
        permRes.push_back((permVT.mz_chi2==-99)?0:permVT.mz_chi2); // VTMZ
        permRes.push_back((permVT.cmc_chi2==-99)?0:permVT.cmc_chi2); //VTCMC
        permRes.push_back((permWSS.mz_pval==std::numeric_limits<double>::infinity())?0:permWSS.mz_chi2); //WSS
        /* std::cout << i << "  " << permRes << std::endl; */
        m_updateCount(permCount, oriRes, permRes, gslr);

        if(i % adaptive != 0 || i == 0) {;}
        else{
            m_updatePvalue(permPval, permCount, i, alpha);
            bool MZok(true), CMCok(true);
            for(UINT t=0; t<permPval.size(); t++){
                if(permPval[t] > 1) {
                    if(t%2) CMCok=false; //CMC series: 1 3
                    else MZok=false;  //MZ series: 0 2 4
                }
            }
            if(__phased){
                if(MZok && CMCok) break; // CMC works for phased only
            } else {
                if(MZok) break;
            }
        }
    }
    for(UINT i=0; i<permCount.size(); i++){
        permPval[i] = (permPval[i]<=1.0)?permPval[i]:(1.0*permCount[i]+1)/(1.0*PermutateTimes+1);
    }
    return permPval;
}


void rareTdt::tdtTest(double lowerMaf, double upperMaf, UINT adaptive, UINT PermutateTimes, double alpha, std::map<std::string, double>& pvalues, std::vector<std::string>& tests, std::string test_mode="A") 
{

	__test_mode = test_mode;

	logRoot["tdtStatic"]["testMode"] = __test_mode;

	// TODO: check the parameters
	RNG rng; gsl_rng* gslr = rng.get();
	/* pvalues.clear(); tests.clear(); */

	for(UINT i=0; i<__popMafs.size(); i++) {
		if(__popMafs[i]<=lowerMaf || __popMafs[i]>upperMaf) 
			__varToBeAnalyzed[i] = false;

	}

	// update _varToBeAnalyzed
	logRoot["tdtStatic"]["siteBeAnalyzed"] = buildJsonArray(__varToBeAnalyzed);

    vectorF np = m_noPermut(); 
	vectorF hapoPval = m_Permut(adaptive, PermutateTimes, alpha, gslr);

	pvalues["CMC-Analytical"] = np[0]; tests.push_back("CMC-Analytical");
	pvalues["BRV-Haplo"] = hapoPval[0]; tests.push_back("BRV-Haplo");
	pvalues["CMC-Haplo"] = hapoPval[1]; tests.push_back("CMC-Haplo");
	pvalues["VT-BRV-Haplo"] = hapoPval[2]; tests.push_back("VT-BRV-Haplo");
	pvalues["VT-CMC-Haplo"] = hapoPval[3];  tests.push_back("VT-CMC-Haplo");
	pvalues["WSS-Haplo"] = hapoPval[4];  tests.push_back("WSS-Haplo");

	return;
}


/**
 * Check if there is a parental transmission bias
 */
void rareTdt::ParentalBiasTest(double lowerMaf, double upperMaf, std::map<std::string, double>& pvalues, std::vector<std::string>& tests)
{
	/* pvalues.clear(); tests.clear(); */

	for(UINT i=0; i<__popMafs.size(); i++) {
		if(__popMafs[i]<=lowerMaf || __popMafs[i]>upperMaf) 
			__varToBeAnalyzed[i] = false;
	}

	// update _varToBeAnalyzed
	logRoot["tdtStatic"]["siteBeAnalyzed"] = buildJsonArray(__varToBeAnalyzed);

	// return: fat_b, fat_c, mot_b, mot_c by CMC 
    	vectorF res = m_parentalTransmissionCount(__tdtBtable, __tdtCtable, __varToBeAnalyzed, true); 

	tests.push_back("pat_untrans");
	tests.push_back("pat_trans");
	tests.push_back("mat_untrans");
	tests.push_back("mat_trans");

	// updated: Mon Nov  3 17:36:26 CST 2014
	pvalues["pat_untrans"] = res[0];
	pvalues["pat_trans"] = res[1];
	pvalues["mat_untrans"] = res[2];
	pvalues["mat_trans"] = res[3];

	/* if(res[0] <= 0 || res[1] <= 0) { */
	/* 	pvalues["pat_or"] = -9;  */
	/* 	pvalues["pat_sd"] = -9;  */
	/* } else { */
	/* 	// http://onlinelibrary.wiley.com/doi/10.1046/J.1469-1809.2005.00156.x/pdf */
	/* 	pvalues["pat_or"] = res[1] * 1.0 / (res[0] * 1.0); */
	/* 	pvalues["pat_sd"] = sqrt(1.0 / res[1] + 1.0 / res[0]); */
	/* 	; */

	/* 	/\* float dir = res[3] > res[1]?1.0:-1.0; *\/ */
	/* 	/\* double pval = myfisher(res); *\/ */
	/* 	/\* pvalues["ParentalBiasTest"] = dir*pval;   *\/ */
	/* 	// "+" means more maternal transmitted events over paternal transmitted events */
	/* } */

	/* if(res[2] <= 0 || res[3] <= 0) { */
	/* 	pvalues["mat_or"] = -9;  */
	/* 	pvalues["mat_sd"] = -9;  */
	/* } else { */
	/* 	pvalues["mat_or"] = res[3] * 1.0 / (res[2] * 1.0); */
	/* 	pvalues["mat_sd"] = sqrt(1.0 / res[3] + 1.0 / res[2]); */
	/* 	// http://onlinelibrary.wiley.com/doi/10.1046/J.1469-1809.2005.00156.x/pdf */
	/* } */

	/* std::cout << res << ", pvalue: " << pvalues["Parental-Bias"] << std::endl; */

	return;
}


// parental transmission count, using CMC method 
vectorF rareTdt::m_parentalTransmissionCount(const vector2F& tdtBtable, const vector2F& tdtCtable, const vectorL& tobeAnalyzed, bool allowDenovo) const
{
	// fat_b, fat_c, mot_b, mot_c
	static const double arr[] = {0,0,0,0};
	vectorF res (arr, arr + sizeof(arr) / sizeof(arr[0]) );

	if(tdtBtable.size() != tdtCtable.size() || tdtBtable[0].size() != tobeAnalyzed.size()) 
		return res;

	for(UINT i=0; i<tdtBtable.size(); i++){

		double fb(0), fc(0), mb(0), mc(0);

		if(tdtBtable[i].size() != tdtCtable[i].size()) {
			vectorF resnull (arr, arr + sizeof(arr) / sizeof(arr[0]) );
			return resnull ;
		}

		for(UINT j=0; j<tdtBtable[i].size(); j++){
			if(tobeAnalyzed[j]) {

				// father
				if(i%2 == 0) {
					fb += tdtBtable[i][j];

					if(tdtCtable[i][j]>=100){
						double offset=(allowDenovo)?100:101;
						fc += tdtCtable[i][j] - offset;
					} else {
						fc += tdtCtable[i][j];
					}
				}
				// mother
				else {
					mb += tdtBtable[i][j];

					if(tdtCtable[i][j]>=100){
						double offset=(allowDenovo)?100:101;
						mc += tdtCtable[i][j] - offset;
					} else {	/*  */
						mc += tdtCtable[i][j];
					}
				}		
			}
		}

		if(fb + fc > 0) {
			res[0] += fb/(fb+fc);
			res[1] += fc/(fb+fc);
		}

		if(mb + mc > 0) {
			res[2] += mb/(mb+mc);
			res[3] += mc/(mb+mc);
		}

	}

	return res;
}


