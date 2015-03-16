#include "gw_utilities.h"
#include "tdtTrio.h"

#define PRINT(a) std::cout << a << std::endl;

/** Base Trio **/
shuffleOutput trio::shuffle(std::string GenoHapo, gsl_rng* gslr){
	shuffleOutput fakeOne;
	return fakeOne;
}

shuffleOutput trio::shuffleWithGenos(std::string GenoHapo, gsl_rng* gslr, vector2F& genos){
	shuffleOutput fakeOne;
	return fakeOne;
}

/** Phased Trio **/
void phasedTrio::load(const vector2F& genotype, const vectorF& popMafs, const vectorL& tobeAnalyzed, bool skipMiss) 
{
	__genotype = genotype;
	__popMafs = popMafs;
	__tobeAnalyzed = tobeAnalyzed;
	__skipMiss = skipMiss;
	__PesudoGenoReady = false;
	__isInformative=true;

	if(__genotype.size()!=6 || __genotype[0].size()!=__popMafs.size() || __tobeAnalyzed.size()!=__popMafs.size()){
		__isInformative=false; return;
	}else {
		for(UINT i=0; i<__genotype.size(); i++) {
			if(__genotype[i].size() != __popMafs.size()) {
				__isInformative=false; return;
			}
		}
	}    

	tdtTableCount();
	return;
}

void phasedTrio::tdtTableCount() 
{
	__tdtTable.assign(4, vectorF(__tobeAnalyzed.size(), 0)); // v[0-1] for father: v[0] is B count; v[1] is C count; v[2-3] for mother
	// label out uninformative sites
	for(UINT i=0; i<__genotype[0].size(); i++){
		if(__skipMiss && (__genotype[0][i]+__genotype[1][i]+__genotype[2][i]+__genotype[3][i] < 0)) {
			__tobeAnalyzed[i]=false; // skip missing
			continue;
		}
		else if((__genotype[0][i]+__genotype[1][i]<0 && __genotype[2][i]+__genotype[3][i]<0) || __genotype[4][i]+__genotype[5][i]<0) {
			__tobeAnalyzed[i]=false; // both parents missing or kid missing
			continue;
		}
		else {}

		for(UINT j =0; j<__genotype.size(); j++){
			if(__genotype[j][i] != 0 && __genotype[j][i] != 1 && __genotype[j][i] != -9) {__tobeAnalyzed[i]=false; break;}
		}
	}   
	// uninformative sites as missing: so it will not be included in unTransmitted site maf calculation
	trimNonSenseSites(__genotype, __tobeAnalyzed);

	// no missing no denovo --> missing no denovo --> missing denovo
	if (m_ChromOrigin(false, false) || m_ChromOrigin(true, false) || m_ChromOrigin(true, true)) {}
	else {
		__isInformative=false;
		return;
	} // should be medlian error

	__denovoSite.assign(__tobeAnalyzed.size(), 0); // 

	int whichFromPat = (__chrOrigin[0] <= 1)?0:1; // which child chromosome comes from father
	__unTransParentalChrs.push_back(__genotype[1-__chrOrigin[whichFromPat]]); // mark the genotype on not analyzed sites as -9
	m_phasedTdtCount(__genotype[__chrOrigin[whichFromPat]], __genotype[1-__chrOrigin[whichFromPat]], __genotype[4 + whichFromPat], 0); // last argument is used to seperate father and mother

	int whichFromMat = 1-whichFromPat;
	__unTransParentalChrs.push_back(__genotype[5-__chrOrigin[whichFromMat]]); // mark the genotype on not analyzed sites as -9
	m_phasedTdtCount(__genotype[__chrOrigin[whichFromMat]], __genotype[5-__chrOrigin[whichFromMat]], __genotype[4 + whichFromMat], 2);
	
	// in case reverse mutation, double mutations
	bool allSiteWrong(true);
	for(UINT i=0; i<__tobeAnalyzed.size(); i++) { if(__tobeAnalyzed[i]) {allSiteWrong=false; break;}}
	if(allSiteWrong) __isInformative=false;
	return;
}

void phasedTrio::m_phasedTdtCount(vectorF& transParent, vectorF& unTransParent, vectorF& Child, UINT base)
{
	for (UINT i = 0; i < Child.size(); i++)
	{
		if(!__tobeAnalyzed[i]) {continue;} //don't analyze this site; __tdtTable[base+0][i] += 0; __tdtTable[base+1][i] += 0;

		double heteroProb = __popMafs[i]*(1-__popMafs[i])*2;
		double homoProb = __popMafs[i]*__popMafs[i];
		double wildProb = (1-__popMafs[i])*(1-__popMafs[i]);
		// missing -9|-9
		if(transParent[i] + unTransParent[i] < 0){
			if(Child[i] == 0) {__tdtTable[base+0][i] += heteroProb/(wildProb+heteroProb);} //-9-9 -> 0 ==> 01 -> 0
			if(Child[i] == 1) {__tdtTable[base+1][i] += heteroProb/(homoProb+heteroProb);} //-9-9 -> 1 ==> 01 -> 1
		} else{
			if(Child[i] == 0){
				if(transParent[i] == 1) {__tobeAnalyzed[i]=false;} //reverse mutation
				else if(transParent[i] == 0 && unTransParent[i] == 1) {__tdtTable[base+0][i] += 1;} //01 -> 0
				else if(transParent[i] == 0 && unTransParent[i] == 0) {}
				else {__tobeAnalyzed[i]=false;} // error
			} else if(Child[i] == 1){
				if(transParent[i] == 1 && unTransParent[i] == 1) {}
				else if(transParent[i] == 1 && unTransParent[i] == 0) {__tdtTable[base+1][i] += 1;}
				else if(transParent[i] == 0) {__tdtTable[base+1][i] += 101; __denovoSite[i]=true;} // de novo mutation
				else {}
			} else {}
		}
	}
	return;
}

bool phasedTrio::m_ChromOrigin(bool allowMissing, bool allowDenovo)
{
	/*! * Determine the origins of Child's chromosomes
	 * Output format: <1,3>, <0,2>.  <9,9> means genotype error */
	vectorF P1(__genotype[0]), P2(__genotype[1]), M1(__genotype[2]), M2(__genotype[3]), C1(__genotype[4]), C2(__genotype[5]);
	__chrOrigin.resize(0);

	int ConeOrigin(9), CtwoOrigin(9);
	char Cone, Ctwo;
	vectorI ConeVect, CtwoVect;

	for (int i = 0; i <4; i++){
		if (ChromEqual(C1, __genotype[i], allowMissing, allowDenovo)) {ConeVect.push_back(1);}
		else {ConeVect.push_back(0);}
		if (ChromEqual(C2, __genotype[i], allowMissing, allowDenovo)) {CtwoVect.push_back(1);}
		else {CtwoVect.push_back(0);}
	}

	// Compare C1 to Parental chromosomes, determine where it comes from:
	if (accumulate(ConeVect.begin(), ConeVect.begin()+4, 0) <= 0) {Cone = 'N';} // not found in either parental 
	else{
		if (accumulate(ConeVect.begin(), ConeVect.begin()+2, 0) > 0 && accumulate(ConeVect.begin()+2, ConeVect.begin()+4, 0) <= 0) {Cone = 'P';} // only found in father
		else if (accumulate(ConeVect.begin(), ConeVect.begin()+2, 0) > 0 && accumulate(ConeVect.begin()+2, ConeVect.begin()+4, 0) > 0) {Cone = 'B';} // found in both parents
		else {Cone = 'M';} // only found in mother
	}
	//The same for C2;
	if (accumulate(CtwoVect.begin(), CtwoVect.begin()+4, 0) <= 0) {Ctwo = 'N';} // not found in either parental 
	else{
		if (accumulate(CtwoVect.begin(), CtwoVect.begin()+2, 0) > 0 && accumulate(CtwoVect.begin()+2, CtwoVect.begin()+4, 0) <= 0) {Ctwo = 'P';} // only found in father
		else if (accumulate(CtwoVect.begin(), CtwoVect.begin()+2, 0) > 0 && accumulate(CtwoVect.begin()+2, CtwoVect.begin()+4, 0) > 0) {Ctwo = 'B';} // found in both parents
		else {Ctwo = 'M';} // only found in mother
	}

	//Combine Cone and Ctwo to determine the origins of child's chromosomes;
	if (Cone == 'N' || Ctwo =='N' || (Cone == 'P' && Ctwo == 'P') || (Cone == 'M' && Ctwo == 'M')){
		return false;
	} //Genotype Errors
	else if (Cone == 'B' && Ctwo == 'B'){
		if (ConeVect[0] == 1) {ConeOrigin = 0;}
		else {ConeOrigin = 1;}
		if (CtwoVect[2] == 1) {CtwoOrigin = 2;}
		else {CtwoOrigin = 3;}
	} //child is homozygous, or parents and child are same genotype (heterozygous for sure);
	else{
		if (Cone == 'P' || Ctwo == 'M'){
			if (ConeVect[0] == 1) {ConeOrigin = 0;}
			else {ConeOrigin = 1;}
			if (CtwoVect[2] == 1) {CtwoOrigin = 2;}
			else {CtwoOrigin = 3;}
		} //Cone comes from Paternal and Ctwo comes from Maternal
		else if (Cone == 'M' || Ctwo == 'P'){
			if (ConeVect[2] == 1) {ConeOrigin = 2;}
			else {ConeOrigin = 3;}
			if (CtwoVect[0] == 1) {CtwoOrigin = 0;}
			else {CtwoOrigin = 1;}
		} //Cone comes from Maternal and Ctwo comes from Paternal
		else {;}
	}
	__chrOrigin.push_back(ConeOrigin); __chrOrigin.push_back(CtwoOrigin);
	return true;
}

bool ChromEqual(const vectorF& ChromOne, const vectorF& ChromTwo, bool allowMissing, bool allowDenovo)
{
	//We allow at most 1 de novo mutation
	if(ChromOne.size() != ChromTwo.size()) return false;
	double AllowedErrorNum = 1;
	int ErrorFound(0);
	bool StrictEqual(true), MissEqual(true), DenovoEqual(true);

	for(UINT i=0; i<ChromOne.size(); i++){
		if(ChromOne[i] != ChromTwo[i]){
			StrictEqual = false;
			if(ChromOne[i] != -9.0 && ChromTwo[i] != -9.0) { MissEqual = false; ErrorFound++;}
		}
	}
	if(ErrorFound > AllowedErrorNum) DenovoEqual = false;

	if((!allowMissing) && (!allowDenovo)) {return StrictEqual;} 
	else if(allowMissing && (!allowDenovo)) {return MissEqual;}
	else if(allowMissing && allowDenovo) {return DenovoEqual;}
	else {return false;}
}

shuffleOutput phasedTrio::shuffle(std::string GenoHapo, gsl_rng* gslr)
{
	if(!__isInformative) {std::cerr << "Error: can't shuffle on non-informative trio." << std::endl; exit(-1);}
	if(__PesudoGenoReady) {;}
	else {m_setUpPesudoGeno(gslr); __PesudoGenoReady=true;}

	vector2F shuffledGeno;
	if(GenoHapo == "GenoShuffle") shuffledGeno = m_GenoShuffle(gslr);
	else if(GenoHapo == "HapoShuffle") shuffledGeno = m_HapoShuffle(gslr);
	else {std::cerr << "Error: invalid shuffle method." << std::endl; exit(-1);}

	phasedTrio shuffleObj;
	shuffleObj.load(shuffledGeno, __popMafs, __tobeAnalyzed, __skipMiss);
	return shuffleOutput(shuffleObj.getTdtTable(), shuffleObj.getNuTransPatChrs());
}

void phasedTrio::m_setUpPesudoGeno(gsl_rng* gslr)
{
    for (UINT i=0; i<4; i++) __pesudoGeno.push_back(__genotype[i]);
    for(UINT j=0; j<__tobeAnalyzed.size(); j++){
        if(!__tobeAnalyzed[j]){
            for(UINT i=0; i<4; i++) __pesudoGeno[i][j]=-9.0; // mark all not analyzed sites as -9
        }
        else if(__genotype[0][j] + __genotype[1][j] + __genotype[2][j] + __genotype[3][j] < 0) // informative but has missing
        {
            if((__genotype[0][j] + __genotype[1][j] < 0) && (__genotype[2][j] + __genotype[3][j] < 0)) for(UINT i=0; i<4; i++) __pesudoGeno[i][j]=-9.0;
            else if(__genotype[0][j] + __genotype[1][j] < 0) { 
                int whichFromPat = (__chrOrigin[0] <= 1)?4:5;
                if(gsl_rng_uniform(gslr) < 0.5) {__pesudoGeno[0][j] = __genotype[whichFromPat][j]; __pesudoGeno[1][j] = -9.0;}
                else {__pesudoGeno[1][j] = __genotype[whichFromPat][j]; __pesudoGeno[0][j] = -9.0;}
            } // father missing
            else if(__genotype[2][j] + __genotype[3][j] < 0) {
                int whichFromMat = (__chrOrigin[0] <= 1)?5:4;
                if(gsl_rng_uniform(gslr) < 0.5) {__pesudoGeno[2][j] = __genotype[whichFromMat][j]; __pesudoGeno[3][j] = -9.0;}
                else {__pesudoGeno[3][j] = __genotype[whichFromMat][j]; __pesudoGeno[2][j] = -9.0;}
            } // mother missing
            else {}
        }
        else {}
    }
    return;
}


/** Can be used by both **/
vector2F trio::m_GenoShuffle(gsl_rng* gslr)
{
	vectorUI recordMissing(__tobeAnalyzed.size(), 0); // only record informative sites: 0 no missing, 1 father missing, 2 mother missing
	vector2F shuffledGeno(__pesudoGeno);
	if(shuffledGeno.size() != 4) {std::cerr << "Error: pesudo genotype for shuffle is not set up properly." << std::endl; exit(-1);}

	for(UINT j=0; j<shuffledGeno[0].size(); j++){
		if(!__tobeAnalyzed[j]) continue;
		vectorF fakeParents;
		for(UINT i=0; i<4; i++){
			if(shuffledGeno[i][j] < 0){
				if(gsl_rng_uniform(gslr) < __popMafs[j]) fakeParents.push_back(1);
				else fakeParents.push_back(0);
				recordMissing[j] = (i<=1)?1:2; // record missing
			}
			else fakeParents.push_back(shuffledGeno[i][j]);
		}
		random_shuffle(fakeParents.begin(), fakeParents.begin()+4);
		for(UINT i=0; i<4; i++) shuffledGeno[i][j] = fakeParents[i];
	}

	if(gsl_rng_uniform(gslr) < 0.5){
		shuffledGeno.push_back(shuffledGeno[gsl_rng_uniform_int(gslr, 2)]);
		shuffledGeno.push_back(shuffledGeno[gsl_rng_uniform_int(gslr, 2) + 2]);
	}
	else{
		shuffledGeno.push_back(shuffledGeno[gsl_rng_uniform_int(gslr, 2) + 2]);
		shuffledGeno.push_back(shuffledGeno[gsl_rng_uniform_int(gslr, 2)]);
	}

	// recode missing
	for(UINT i=0; i<recordMissing.size(); i++){
		if (recordMissing[i] == 0) {;}
		else if(recordMissing[i] == 1) {shuffledGeno[0][i] = -9; shuffledGeno[1][i] = -9;}
		else if(recordMissing[i] == 2) {shuffledGeno[2][i] = -9; shuffledGeno[3][i] = -9;}
		else {;}    
	}

	if(shuffledGeno.size() != 6) {std::cerr << "Error: genotype shuffle goes wrong." << std::endl; exit(-1);}
	else return shuffledGeno;
}

// no hapotype permutation for unphased data
vector2F trio::m_HapoShuffle(gsl_rng* gslr)
{
    vectorUI recordMissing(__tobeAnalyzed.size(), 0); // only record informative sites: 0 no missing, 1 father missing, 2 mother missing
    vector2F shuffledGeno(__pesudoGeno);
    if(shuffledGeno.size() != 4) {std::cerr << "Error: pesudo genotype for shuffle is not set up properly." << std::endl; exit(-1);}
    
    // refill missing site
    for(UINT j=0; j<shuffledGeno[0].size(); j++){
        if(!__tobeAnalyzed[j]) continue;
        
        if(shuffledGeno[0][j] < 0 || shuffledGeno[1][j] < 0 || shuffledGeno[2][j] < 0 || shuffledGeno[3][j] < 0){
            if((shuffledGeno[0][j] < 0 && shuffledGeno[1][j] < 0) || (shuffledGeno[2][j] < 0 && shuffledGeno[3][j] < 0)){std::cerr << "Error: pesudo genotype wrong, no guessed genotype in missing sites." << std::endl; exit(-1);}
            // both parent missing
            if(shuffledGeno[0][j] + shuffledGeno[1][j] + shuffledGeno[2][j] + shuffledGeno[3][j] < -9.0){std::cerr << "Error: pesudo genotype wrong, both missing are analyzed here" << std::endl; exit(-1);}
            
            int base = (shuffledGeno[0][j] + shuffledGeno[1][j] < 0)?0:2; // which parent is missing
            recordMissing[j] = (base == 0)?1:2;
            
            double infered = (shuffledGeno[base+0][j] < 0)?shuffledGeno[base+1][j]:shuffledGeno[base+0][j];
            double guessed = (gsl_rng_uniform(gslr) < __popMafs[j])?1:0;
            if(gsl_rng_uniform(gslr) < 0.5) {shuffledGeno[base+0][j] = infered; shuffledGeno[base+1][j] = guessed;}
            else {shuffledGeno[base+0][j] = guessed; shuffledGeno[base+1][j] = infered;}
        }
    }
    
    random_shuffle(shuffledGeno.begin(), shuffledGeno.begin()+4);

    if(gsl_rng_uniform(gslr) < 0.5){
        shuffledGeno.push_back(shuffledGeno[gsl_rng_uniform_int(gslr, 2)]);
        shuffledGeno.push_back(shuffledGeno[gsl_rng_uniform_int(gslr, 2) + 2]);
    }
    else{
        shuffledGeno.push_back(shuffledGeno[gsl_rng_uniform_int(gslr, 2) + 2]);
        shuffledGeno.push_back(shuffledGeno[gsl_rng_uniform_int(gslr, 2)]);
    }
    // recode missing
    for(UINT i=0; i<recordMissing.size(); i++){
        if (recordMissing[i] == 0) {;}
        else if(recordMissing[i] == 1) {shuffledGeno[0][i] = -9; shuffledGeno[1][i] = -9;}
        else if(recordMissing[i] == 2) {shuffledGeno[2][i] = -9; shuffledGeno[3][i] = -9;}
        else {;}    
    }
    
    if(shuffledGeno.size() != 6) {std::cerr << "Error: hapotype shuffle goes wrong." << std::endl; exit(-1);}
    else return shuffledGeno;
}


/** Unphased Trio **/
void unPhasedTrio::load(const vector2F& genotype, const vectorF& popMafs, const vectorL& tobeAnalyzed, bool skipMiss)
{
	__genotype = genotype;
	__popMafs = popMafs;
	__tobeAnalyzed = tobeAnalyzed;
	__skipMiss = skipMiss;
	__PesudoGenoReady = false;
	__isInformative=true;
	//data check
	if(__genotype.size() != 3 || __genotype[0].size() != __popMafs.size() || __tobeAnalyzed.size() != __popMafs.size()){
		__isInformative=false; return;
	}
	else {
		for(UINT i=0; i<__genotype.size(); i++) {
			if(__genotype[i].size() != __popMafs.size()) {
				__isInformative=false; return;
			}
		}
	}
	tdtTableCount();
	return;
}

void unPhasedTrio::tdtTableCount()
{
	__tdtTable.resize(2); 
	__tdtTable[0].assign(__tobeAnalyzed.size(), 0); __tdtTable[1].assign(__tobeAnalyzed.size(), 0); //__tdtTable[2].assign(__tobeAnalyzed.size(), 0); __tdtTable[3].assign(__tobeAnalyzed.size(), 0); 
	__unTransParentalChrs.resize(2);
	__unTransParentalChrs[0].assign(__tobeAnalyzed.size(), -9.0); __unTransParentalChrs[1].assign(__tobeAnalyzed.size(), -9.0);
	// label out uninformative sites
	for(UINT i=0; i<__genotype[0].size(); i++){
		if(__genotype[0][i] < 0 && __genotype[1][i] < 0) __tobeAnalyzed[i]=false; // both parents missing
		else if(__genotype[2][i] < 0) __tobeAnalyzed[i]=false; // kid missing
		else if(__skipMiss){
			if(__genotype[0][i]+__genotype[1][i]+__genotype[2][i] < 0) __tobeAnalyzed[i]=false; // individual missing
		}
		else {;}

		for(UINT j =0; j<__genotype.size(); j++){
			if(__genotype[j][i] != 0 && __genotype[j][i] != 1 && __genotype[j][i] != 2 && __genotype[j][i] != -9) __tobeAnalyzed[i]=false;
		}
	}
	// uninformative sites as missing
	trimNonSenseSites(__genotype, __tobeAnalyzed);

	__denovoSite.assign(__tobeAnalyzed.size(), 0);    
	for(UINT i=0; i<__genotype[0].size(); i++){
		if(!__tobeAnalyzed[i]) {continue;}
		m_unPhasedTdtCount(__genotype[0][i], __genotype[1][i], __genotype[2][i], i);
	}
	// in case reverse mutation, double mutations
	bool allSiteWrong(true);
	for(UINT i=0; i<__tobeAnalyzed.size(); i++) { if(__tobeAnalyzed[i]) {allSiteWrong=false; break;}}
	if(allSiteWrong) __isInformative=false;
	return;
}

// modify __tdtTable, __unTransParentalChrs, __denovoSite and __tobeAnalyzed; through index
// only analyze valid site, that is no double missing, no kid missing, no coding error
void unPhasedTrio::m_unPhasedTdtCount(double fat, double mot, double kid, UINT index)
{
    if(!__tobeAnalyzed[index]) {return;}
    // no missing
    if(fat >= 0 && mot >= 0){
        if(fat == 0){
            if(mot == 0){
                if(kid == 0) {__unTransParentalChrs[0][index] = 0; __unTransParentalChrs[1][index] = 0;}
                else if(kid == 1) {__unTransParentalChrs[0][index] = 0; __unTransParentalChrs[1][index] = 0; __tdtTable[1][index] += 101; __denovoSite[index]=true;} // denovo
                else {__tobeAnalyzed[index]=false;} // no double mutation
            }
            else if(mot == 1){
                if(kid == 0) {__unTransParentalChrs[0][index] = 0; __unTransParentalChrs[1][index] = 1; __tdtTable[0][index] += 1;}
                else if(kid == 1) {__unTransParentalChrs[0][index] = 0; __unTransParentalChrs[1][index] = 0; __tdtTable[1][index] += 1;}
                else if(kid == 2) {__unTransParentalChrs[0][index] = 0; __unTransParentalChrs[1][index] = 0; __tdtTable[1][index] += 102; __denovoSite[index]=true;} //transmitted and denovo
                else {;}
            }
            else if(mot == 2){
                if(kid == 0) {__tobeAnalyzed[index]=false;} // no reverse mutation
                else if(kid == 1) {__unTransParentalChrs[0][index] = 0; __unTransParentalChrs[1][index] = 1;}
                else if(kid == 2) {__unTransParentalChrs[0][index] = 0; __unTransParentalChrs[1][index] = 1; __tdtTable[1][index] += 101;__denovoSite[index]=true;} // denovo 
                else {;}
            }
            else {;}
        }
        else if(fat == 2){
            if(mot == 0){
                if(kid == 0) {__tobeAnalyzed[index]=false;} // no reverse mutation
                else if(kid == 1) {__unTransParentalChrs[0][index] = 1; __unTransParentalChrs[1][index] = 0;}
                else if(kid == 2) {__unTransParentalChrs[0][index] = 1; __unTransParentalChrs[1][index] = 0; __tdtTable[1][index] += 101;__denovoSite[index]=true;} // denovo
                else {;}
            }
            else if(mot == 1){
                if(kid == 0) {__tobeAnalyzed[index]=false;}
                else if(kid == 1) {__unTransParentalChrs[0][index] = 1; __unTransParentalChrs[1][index] = 1; __tdtTable[0][index] += 1;}
                else if(kid == 2) {__unTransParentalChrs[0][index] = 1; __unTransParentalChrs[1][index] = 0; __tdtTable[1][index] += 1;}
                else {;}
            }
            else if(mot == 2){
                if(kid == 2) {__unTransParentalChrs[0][index] = 1; __unTransParentalChrs[1][index] = 1;}
                else {__tobeAnalyzed[index]=false;}
            }
            else {;}
        }
        else if(fat == 1){
            if(mot == 0){
                if(kid == 0) {__unTransParentalChrs[0][index] = 1; __unTransParentalChrs[1][index] = 0; __tdtTable[0][index] += 1;}
                else if(kid == 1) {__unTransParentalChrs[0][index] = 0; __unTransParentalChrs[1][index] = 0; __tdtTable[1][index] += 1;}
                else if(kid == 2) {__unTransParentalChrs[0][index] = 0; __unTransParentalChrs[1][index] = 0; __tdtTable[1][index] += 102; __denovoSite[index]=true;} // denovo
                else {;}
            }
            else if(mot == 1){
                if(kid == 0) {__unTransParentalChrs[0][index] = 1; __unTransParentalChrs[1][index] = 1; __tdtTable[0][index] += 2;}
                else if(kid == 1) {__unTransParentalChrs[0][index] = 1; __unTransParentalChrs[1][index] = 0; __tdtTable[0][index] += 1; __tdtTable[1][index] += 1;}
                else if(kid == 2) {__unTransParentalChrs[0][index] = 0; __unTransParentalChrs[1][index] = 0; __tdtTable[1][index] += 2;}
                else {;}
            }
            else if(mot == 2){
                if(kid == 0) {__tobeAnalyzed[index]=false;}
                else if(kid == 1) {__unTransParentalChrs[0][index] = 1; __unTransParentalChrs[1][index] = 1; __tdtTable[0][index] += 1;}
                else if(kid == 2) {__unTransParentalChrs[0][index] = 0; __unTransParentalChrs[1][index] = 1; __tdtTable[1][index] += 1;}
                else {;}
            }
            else {;}
	}
        else {;}
    }
    // one parent missing: no denovo mutation when missing
    else if(fat == -9 || mot == -9){
        if(fat == -9 && mot == -9) {__tobeAnalyzed[index]=false; return;}
        double missing = (fat==-9.0)?fat:mot;
        mot = (fat==-9.0)?mot:fat; fat = missing; // assume missing father
        
        double heteroProb = __popMafs[index]*(1-__popMafs[index])*2;
        double homoProb = __popMafs[index]*__popMafs[index];
        double wildProb = (1-__popMafs[index])*(1-__popMafs[index]);
        
        if(mot == 0){
            if(kid == 0) {__unTransParentalChrs[1][index] = 0; __tdtTable[0][index] += heteroProb/(wildProb+heteroProb);} // -9X0 -> 0 ==> 1X0 -> 0
            else if(kid == 1) {__unTransParentalChrs[1][index] = 0; __tdtTable[1][index] += heteroProb/(homoProb+heteroProb);} // -9X0 -> 1 ==> 1X0 -> 1
            else if(kid == 2) {__tobeAnalyzed[index]=false;}
            else {;}
        }
        else if(mot == 1){
            if(kid == 0) {__unTransParentalChrs[1][index] = 1; __tdtTable[0][index] += (1+heteroProb/(wildProb+heteroProb));} // -9X1 -> 0 ==> 1X1 -> 0 or 0X1 -> 0
            else if(kid == 1) {__unTransParentalChrs[1][index] = 0; __tdtTable[0][index] += homoProb+heteroProb; __tdtTable[1][index] += wildProb+heteroProb;} // -9X1 -> 1 ==> 0X1 -> 1 or 1X1 -> 1 or 2X1 -> 1
            else if(kid == 2) {__unTransParentalChrs[1][index] = 0; __tdtTable[1][index] += (1+heteroProb/(homoProb+heteroProb));} // -9X1 -> 2 ==> 1X1 -> 2 or 2X1 -> 2
            else {;}
        }
        else if(mot == 2){
            if(kid == 0) {__tobeAnalyzed[index]=false;} 
            else if(kid == 1) {__unTransParentalChrs[1][index] = 1; __tdtTable[0][index] += heteroProb/(wildProb+heteroProb);} // -9X2 -> 1 ==> 1X2 -> 1
            else if(kid == 2) {__unTransParentalChrs[1][index] = 1; __tdtTable[1][index] += heteroProb/(homoProb+heteroProb);} // -9X2 -> 2 ==> 1X2 -> 2
            else {;}
        }
        else {;}
    }
    else{;}
}


shuffleOutput unPhasedTrio::shuffle(std::string GenoHapo, gsl_rng* gslr)
{
	if(!__isInformative) {std::cerr << "Error: can't shuffle on non-informative trio." << std::endl; exit(-1);}
	if(__PesudoGenoReady) {;}
	else { m_setUpPesudoGeno(gslr); __PesudoGenoReady=true;}
	vector2F shuffledGeno;
	if(GenoHapo == "GenoShuffle") shuffledGeno = m_GenoShuffle(gslr);
	else if(GenoHapo == "HapoShuffle") shuffledGeno = m_HapoShuffle(gslr);
	else {std::cerr << "Error: invalid shuffle method." << std::endl; exit(-1);}
	// phased to unphased conversion
	vector2F geno(shuffledGeno.begin(), shuffledGeno.begin()+3);
	for(UINT j=0; j<shuffledGeno[0].size(); j++){
		for(UINT i=0; i<3; i++){
			if(shuffledGeno[2*i][j] < 0 || shuffledGeno[2*i+1][j] < 0) geno[i][j] = -9.0;
			else geno[i][j] = shuffledGeno[2*i][j] + shuffledGeno[2*i+1][j];
		}
	}
	unPhasedTrio shuffleObj;
	shuffleObj.load(geno, __popMafs, __tobeAnalyzed, __skipMiss);
	return shuffleOutput(shuffleObj.getTdtTable(), shuffleObj.getNuTransPatChrs());
}

// for debug: shuffle and return shuffled genotypes
shuffleOutput unPhasedTrio::shuffleWithGenos(std::string GenoHapo, gsl_rng* gslr, vector2F& genos)
{
	if(!__isInformative) {std::cerr << "Error: can't shuffle on non-informative trio." << std::endl; exit(-1);}
	if(__PesudoGenoReady) {;}
	else { m_setUpPesudoGeno(gslr); __PesudoGenoReady=true;}
	vector2F shuffledGeno;
	if(GenoHapo == "GenoShuffle") shuffledGeno = m_GenoShuffle(gslr);
	else if(GenoHapo == "HapoShuffle") shuffledGeno = m_HapoShuffle(gslr);
	else {std::cerr << "Error: invalid shuffle method." << std::endl; exit(-1);}
	// phased to unphased conversion
	vector2F geno(shuffledGeno.begin(), shuffledGeno.begin()+3);
	for(UINT j=0; j<shuffledGeno[0].size(); j++){
		for(UINT i=0; i<3; i++){
			if(shuffledGeno[2*i][j] < 0 || shuffledGeno[2*i+1][j] < 0) geno[i][j] = -9.0;
			else geno[i][j] = shuffledGeno[2*i][j] + shuffledGeno[2*i+1][j];
		}
	}
	genos = geno;
	unPhasedTrio shuffleObj;
	shuffleObj.load(geno, __popMafs, __tobeAnalyzed, __skipMiss);
	return shuffleOutput(shuffleObj.getTdtTable(), shuffleObj.getNuTransPatChrs());
}


void unPhasedTrio::m_setUpPesudoGeno(gsl_rng* gslr)
{
    //for(UINT i=0; i<2; i++) std::cout << __genotype[i] << std::endl;
	__pesudoGeno.resize(4);
	for(UINT i=0; i<4; i++) __pesudoGeno[i].resize(__genotype[0].size());
	for(UINT j=0; j<__tobeAnalyzed.size(); j++){
		if(!__tobeAnalyzed[j]){
			for(UINT i=0; i<4; i++) __pesudoGeno[i][j]=-9.0; // mark all not analyzed sites as -9
		}
		else if(__genotype[0][j] < 0 && __genotype[1][j] < 0) { for(UINT i=0; i<4; i++) __pesudoGeno[i][j]=-9.0;}    
		else if(__genotype[0][j] < 0 && __genotype[1][j] >= 0)
		{
			double guessed = guessMissing(__genotype[1][j], __genotype[2][j], gslr);
			if(guessed < 0) {for(UINT i=0; i<4; i++) __pesudoGeno[i][j]=-9.0; __tobeAnalyzed[j] = false;}
			else{
				if(gsl_rng_uniform(gslr) > 0.5) { __pesudoGeno[0][j] = guessed;  __pesudoGeno[1][j] = -9.0;}
				else { __pesudoGeno[1][j] = guessed;  __pesudoGeno[0][j] = -9.0;}
			}

			if(__genotype[1][j] == 0) {__pesudoGeno[2][j] = 0; __pesudoGeno[3][j] = 0;}
			else if(__genotype[1][j] == 1) {__pesudoGeno[2][j] = 0; __pesudoGeno[3][j] = 1;}
			else if(__genotype[1][j] == 2) {__pesudoGeno[2][j] = 1; __pesudoGeno[3][j] = 1;}
			else {for(UINT i=0; i<4; i++) __pesudoGeno[i][j]=-9.0; __tobeAnalyzed[j] = false;}
		}    
		else if(__genotype[0][j] >= 0 && __genotype[1][j] < 0)
		{
			double guessed = guessMissing(__genotype[0][j], __genotype[2][j], gslr);
			if(guessed < 0) {for(UINT i=0; i<4; i++) __pesudoGeno[i][j]=-9.0; __tobeAnalyzed[j] = false;}
			else{
				if(gsl_rng_uniform(gslr) > 0.5) {__pesudoGeno[2][j] = guessed; __pesudoGeno[3][j] = -9.0;}
				else {__pesudoGeno[3][j] = guessed; __pesudoGeno[2][j] = -9.0;}
			}

			if(__genotype[0][j] == 0) {__pesudoGeno[0][j] = 0; __pesudoGeno[1][j] = 0;}
			else if(__genotype[0][j] == 1) {__pesudoGeno[0][j] = 0; __pesudoGeno[1][j] = 1;}
			else if(__genotype[0][j] == 2) {__pesudoGeno[0][j] = 1; __pesudoGeno[1][j] = 1;}
			else {for(UINT i=0; i<4; i++) __pesudoGeno[i][j]=-9.0; __tobeAnalyzed[j] = false;}
		} 
		else if(__genotype[0][j] >= 0 && __genotype[1][j] >= 0)
		{
			if(__genotype[0][j] == 0) {__pesudoGeno[0][j] = 0; __pesudoGeno[1][j] = 0;}
			else if(__genotype[0][j] == 1) {__pesudoGeno[0][j] = 0; __pesudoGeno[1][j] = 1;}
			else if(__genotype[0][j] == 2) {__pesudoGeno[0][j] = 1; __pesudoGeno[1][j] = 1;}
			else {for(UINT i=0; i<4; i++) __pesudoGeno[i][j]=-9.0; __tobeAnalyzed[j] = false;}

			if(__genotype[1][j] == 0) {__pesudoGeno[2][j] = 0; __pesudoGeno[3][j] = 0;}
			else if(__genotype[1][j] == 1) {__pesudoGeno[2][j] = 0; __pesudoGeno[3][j] = 1;}
			else if(__genotype[1][j] == 2) {__pesudoGeno[2][j] = 1; __pesudoGeno[3][j] = 1;}
			else {for(UINT i=0; i<4; i++) __pesudoGeno[i][j]=-9.0; __tobeAnalyzed[j] = false;}
		} // no missing
		else {;}
	}
    //std::cout << "-----------\n";
    //for(UINT i=0; i<4; i++) std::cout << __pesudoGeno[i] << std::endl;
    //std::cout << "===========\n";
	return;
}

double unPhasedTrio::guessMissing(double knownParent, double kid, gsl_rng* gslr)
{
    double guessedHapo(-1);
    if(knownParent == 0 && kid == 0) guessedHapo = 0;
    else if(knownParent == 0 && kid == 1) guessedHapo = 1;
    else if(knownParent == 0 && kid == 2) guessedHapo = -1; // 0X-9 -> 2
    else if(knownParent == 1 && kid == 0) guessedHapo = 0;
    else if(knownParent == 1 && kid == 1) guessedHapo = (gsl_rng_uniform(gslr) > 0.5)?1:0; // TODO: should be more explicit
    else if(knownParent == 1 && kid == 2) guessedHapo = 1;
    else if(knownParent == 2 && kid == 0) guessedHapo = -1; // 2X-9 ->0
    else if(knownParent == 2 && kid == 1) guessedHapo = 0;
    else if(knownParent == 2 && kid == 2) guessedHapo = 1;
    else {;}
    return guessedHapo;
}
