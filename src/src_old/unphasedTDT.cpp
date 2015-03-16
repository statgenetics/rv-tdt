/*!\file unphased.cpp
   \brief Deal with unphased tios

   Copyright: 2011 Zongxiao He @ Baylor College of Medicine
*/

#include "unphasedTDT.h"
#include "gsl/gsl_cdf.h"
#include "phasedTDT.h"
#include "tdt_tools.h"
#include "rareTdt.h"
bool DEBUG = false;
bool VTDEBUG = false;
bool permut = false;

unphasedTDT::unphasedTDT(const vector2F& genotypes, const vectorF& popMafs, const vectorF& samMafs, std::string& chrom, double siteMissRatio, bool SkipMiss)
{
    __twoDGenos = genotypes;
    __popMafs = popMafs;
    __samMafs = samMafs;
    __chr = chrom;
    __wrongFam = 0;
    __SkipMiss = SkipMiss;
    __Xvariants = 0;
    __varNum = __twoDGenos[0].size();
    __lowMissing.assign(__varNum, true);
    __discardFam.assign(__varNum, 0);
    __denovo.assign(__varNum, 0);
    __flip.assign(__varNum, false);
    
    bool NeedFlip(false);
    for(UINT i=0; i<__samMafs.size(); i++){
	if(__samMafs[i] > 0.5) {
	    __flip[i] = true;
	    NeedFlip = true;
	    __samMafs[i] = 1-__samMafs[i];
	}
    }
    if(NeedFlip) flipVector2F();
    
    m_SiteMissRatio();
    if(siteMissRatio > 0){
	for(UINT i=0; i<__varNum; i++) {
	    if(__MissingRatio[i] >= siteMissRatio) {__lowMissing[i] = false; __Xvariants++;}
	}
    }
    __varNumToBeAnalyzed = __varNum - __Xvariants;
    __trioGenos = m_geTrio3D(__twoDGenos, 3);

    m_Count();
    if(__SkipMiss) m_markMiss();
    // erase the sites that shouldn't be analyzed, because missing
    m_TrimTdtTable();
    return;
}

void unphasedTDT::flipVector2F()
{
  for (UINT i=0; i<__twoDGenos.size(); i++){
    for (UINT j=0; j < __twoDGenos[i].size(); j++){
	if(__flip[j]) __twoDGenos[i][j] = m_flip(__twoDGenos[i][j]);
      }
    }
  return;
}

float unphasedTDT::m_flip(float before)
{
    float after;
    if(before >= 0 && before <= 2) after = 2-before;
    else after = before;
    return after;
}

void unphasedTDT::m_markMiss()
{
    for(UINT i=0; i<__trioGenos.size(); i++){
	//walk through each site
	for(UINT j=0; j<__trioGenos[i][0].size(); j++){
	    if(__trioGenos[i][0][j]+__trioGenos[i][1][j]+__trioGenos[i][2][j] < 0) //missing
	    {
		for(UINT k=0; k<3; k++) __trioGenos[i][k][j] = 0;
	    }
	}
    }
}

void unphasedTDT::m_SiteMissRatio()
{
    //observed missing ratio
    vectorUI MissingCount(__twoDGenos[0].size(),0);
    for (UINT i = 0; i < __twoDGenos.size(); i++){
        for (UINT j = 0; j < __twoDGenos[i].size(); j++){
            if ( __twoDGenos[i][j] < 0.0) {MissingCount[j] ++;}
        }
    }
    
    for(UINT i = 0; i < MissingCount.size(); i++) __MissingRatio.push_back(float(MissingCount[i])/__twoDGenos.size());
    return;
}

void unphasedTDT::m_TrimTdtTable()
{    
    for(UINT i=0; i<__TdtCountB.size(); i++){
	for(UINT j=0; j<__TdtCountB[i].size(); j++){
	    if(!__lowMissing[j]) {__TdtCountB[i][j] = 0; __TdtCountC[i][j] = 0;}
	}
    }
    return;
}

// Remeber: the data for GenoPermut may has less variant sites compared to original one, because those sites that have too much missing or uninformative have been trimed;
unphasedTDT::unphasedTDT(const vector3F& InputData, const vectorF& popMafs, const vectorF& samMafs, std::string& chrom, bool SkipMiss)
{
    __samMafs = samMafs;
    __popMafs = popMafs;
    __chr = chrom;
    __wrongFam = 0;
    __SkipMiss = SkipMiss;
    __trioGenos = InputData;
    __varNum = __trioGenos[0][0].size();
    __lowMissing.assign(__varNum, true);
    __discardFam.assign(__varNum, 0);
    __flip.assign(__varNum, false);
    
    bool NeedFlip(false);
    for(UINT i=0; i<__samMafs.size(); i++){
	if(__samMafs[i] > 0.5) {
	    __flip[i] = true;
	    NeedFlip = true;
	    __samMafs[i] = 1-__samMafs[i];
	}
    }
    if(NeedFlip) flipVector2F();
    
    m_Count();
    return;
}

unphasedTDT::~unphasedTDT(){;}

void unphasedTDT::m_Count()
{
    // tdt count table for this gene
    for(UINT f=0; f < __trioGenos.size(); f++)
    {
	//if(permut) std::cout << "------" << std::endl;
	//if(permut) for(UINT i=0; i< __trioGenos[f].size(); i++) std::cout << __trioGenos[f][i] << std::endl;
	
	//entire parent missing should be handled over here
	bool entireMissingFather(true), entireMissingMother(true);
	//in case both parents missing
	for(UINT j=0; j<__trioGenos[f][0].size(); j++){
	    if(__lowMissing[j]) {
		if(__trioGenos[f][0][j] >= 0) entireMissingFather = false;
		if(__trioGenos[f][1][j] >= 0) entireMissingMother = false;
	    }
	}
	
	// within one family: (for unphased data, using 0,1,2,-9 coding)
	vectorF oneTrioB, oneTrioC, fatUNtrans, motUNtrans;
	for(UINT s=0; s < __trioGenos[f][0].size(); s++)
	{
	    double fat(__trioGenos[f][0][s]), mot(__trioGenos[f][1][s]), kid(__trioGenos[f][2][s]);
	    
	    vectorF untransmitted = m_untransmitted(fat, mot, kid);
	    fatUNtrans.push_back(untransmitted[0]); motUNtrans.push_back(untransmitted[1]);

	    if(__SkipMiss){
		if(fat == -9 || mot == -9 || kid == -9) {oneTrioB.push_back(0); oneTrioC.push_back(0); continue;}
	    }
	    //if(__flip[s]) {fat = m_flip(fat); mot = m_flip(mot); kid = m_flip(kid);} //flipped at begining of the class
	    
	    if((fat != 0 && fat != 1 && fat != 2 && fat != -9) || (mot != 0 && mot != 1 && mot != 2 && mot != -9) || (kid != 0 && kid != 1 && kid != 2 && kid != -9)){
		oneTrioB.push_back(0); oneTrioC.push_back(0);
		// uninformative site
		__discardFam[s]++;
	    }
	    else{
		// missing on child's genotype, or both parents'
		if(kid == -9) {oneTrioB.push_back(0); oneTrioC.push_back(0);__discardFam[s]++;}
		else if(fat == -9 && mot == -9) {oneTrioB.push_back(0); oneTrioC.push_back(0);__discardFam[s]++;}
		
		// missing on one parent's genotype
		else if(fat == -9 || mot == -9){		    
		    // don't consider de novo mutation if there is missing
		    if(fat == -9){
			vectorC Vzero = m_tdtCount(0, mot, kid, false);
			vectorC Vone = m_tdtCount(1, mot, kid, false);
			vectorC Vtwo = m_tdtCount(2, mot, kid, false);
			
			// there is de novo mutation whatever father (missing) genotype is. eg 0 X -9 --> 2
			if(IsIt(Vzero[0], 'm') + IsIt(Vone[0], 'm') + IsIt(Vtwo[0], 'm') == 3) {oneTrioB.push_back(0); oneTrioC.push_back(0); __discardFam[s]++;}
			else
			{
			    // tricky part: don't consider de novo mutation when missing
			    double Pzero = (1-__popMafs[s])*(1-__popMafs[s])*(1-IsIt(Vzero[0], 'm'));
			    double Pone = 2*(1-__popMafs[s])*(__popMafs[s])*(1-IsIt(Vone[0], 'm'));
			    double Ptwo = (__popMafs[s])*(__popMafs[s])*(1-IsIt(Vtwo[0], 'm'));
			    
			    double weight = Pzero + Pone + Ptwo;
			    double Bcount, Ccount;
			    // if want add gender information, please modify this part
			    if(__chr == "X") Bcount = Pzero*(IsIt(Vzero[1], 'b')) + Pone*(IsIt(Vone[1], 'b')) + Ptwo*(IsIt(Vtwo[1], 'b'));
			    else Bcount = Pzero*(IsIt(Vzero[0], 'b') + IsIt(Vzero[1], 'b')) + Pone*(IsIt(Vone[0], 'b') + IsIt(Vone[1], 'b')) + Ptwo*(IsIt(Vtwo[0], 'b') + IsIt(Vtwo[1], 'b'));
			    if(__chr == "X") Ccount = Pzero*(IsIt(Vzero[1], 'c')) + Pone*(IsIt(Vone[1], 'c')) + Ptwo*(IsIt(Vtwo[1], 'c'));
			    else Ccount = Pzero*(IsIt(Vzero[0], 'c') + IsIt(Vzero[1], 'c')) + Pone*(IsIt(Vone[0], 'c') + IsIt(Vone[1], 'c')) + Ptwo*(IsIt(Vtwo[0], 'c') + IsIt(Vtwo[1], 'c'));
			    
			    oneTrioB.push_back(Bcount/weight);
			    oneTrioC.push_back(Ccount/weight);
			}
		    }
		    else if(mot == -9){
			vectorC Vzero = m_tdtCount(fat, 0, kid, false);
			vectorC Vone = m_tdtCount(fat, 1, kid, false);
			vectorC Vtwo = m_tdtCount(fat, 2, kid, false);
			
			// there is de novo mutation whatever father (missing) genotype is. eg 0 X -9 --> 2
			if(IsIt(Vzero[0], 'm') + IsIt(Vone[0], 'm') + IsIt(Vtwo[0], 'm') == 3) {oneTrioB.push_back(0); oneTrioC.push_back(0); __discardFam[s]++;}
			else{
			    // tricky part: don't consider de novo mutation when missing
			    double Pzero = (1-__popMafs[s])*(1-__popMafs[s])*(1-IsIt(Vzero[0], 'm'));
			    double Pone = 2*(1-__popMafs[s])*(__popMafs[s])*(1-IsIt(Vone[0], 'm'));
			    double Ptwo = (__popMafs[s])*(__popMafs[s])*(1-IsIt(Vtwo[0], 'm'));
			    
			    double weight = Pzero + Pone + Ptwo;
			    double Bcount, Ccount;
			    if(__chr == "X") Bcount = Pzero*(IsIt(Vzero[1], 'b')) + Pone*(IsIt(Vone[1], 'b')) + Ptwo*(IsIt(Vtwo[1], 'b'));
			    else Bcount = Pzero*(IsIt(Vzero[0], 'b') + IsIt(Vzero[1], 'b')) + Pone*(IsIt(Vone[0], 'b') + IsIt(Vone[1], 'b')) + Ptwo*(IsIt(Vtwo[0], 'b') + IsIt(Vtwo[1], 'b'));
			    if(__chr == "X") Ccount = Pzero*(IsIt(Vzero[1], 'c')) + Pone*(IsIt(Vone[1], 'c')) + Ptwo*(IsIt(Vtwo[1], 'c'));
			    else Ccount = Pzero*(IsIt(Vzero[0], 'c') + IsIt(Vzero[1], 'c')) + Pone*(IsIt(Vone[0], 'c') + IsIt(Vone[1], 'c')) + Ptwo*(IsIt(Vtwo[0], 'c') + IsIt(Vtwo[1], 'c'));
			    
			    oneTrioB.push_back(Bcount/weight);
			    oneTrioC.push_back(Ccount/weight);
			}
		    }
		    // will not happen
		    else {__discardFam[s]++;}
		}
		// no missing
		else{
		    vectorC BandC = m_tdtCount(fat, mot, kid, true);
		    double Bcount, Ccount;
		    if(__chr == "X") Bcount = IsIt(BandC[1], 'b');
		    else Bcount = IsIt(BandC[0], 'b') + IsIt(BandC[1], 'b');
		    
		    if(__chr == "X") Ccount = IsIt(BandC[1], 'c');
		    else Ccount = IsIt(BandC[0], 'c') + IsIt(BandC[1], 'c');
		    oneTrioB.push_back(Bcount);
		    oneTrioC.push_back(Ccount);
		    if(IsIt(BandC[0], 'm')  + IsIt(BandC[1], 'm')  == 2) { __discardFam[s]++;}
		    //record denovo mutation, only when no missing no genotype error
		    if(IsIt(BandC[2], 'y')) __denovo[s]++;
		}
	    }
	}
	__NonTransChrom.push_back(fatUNtrans);
	__NonTransChrom.push_back(motUNtrans);
	if(entireMissingFather||entireMissingMother) continue;
	__TdtCountB.push_back(oneTrioB);
	__TdtCountC.push_back(oneTrioC);
	//if(permut) std::cout << oneTrioB << std::endl << oneTrioC << std::endl;
    }
    return;
}


// determine the tdt count within one trio for one site
// de nove == true, means consider de novo mutation here
vectorC unphasedTDT::m_tdtCount(double fat, double mot, double kid, bool allowDenovo)
{
    RNG rng;
    gsl_rng* gslr = rng.get();
    char fatBC('n'), motBC('n');
    //all mutations, may not informative, e.g. reverse mutation, double mutation
    bool mutationHappened(false);
    // mutation events that counts as C
    bool denovo(false);
    
    // the kid is 00
    if(kid == 0){
	// for father
	if(fat == 2) {mutationHappened = true;} // reverse mutation
	else if(fat == 1) {fatBC = 'b';} // father count as b
	else if(fat == 0) {;} // uninformative
	else {;}

	// for mother
	if(mot == 2) {mutationHappened = true;} // reverse mutation
	else if(mot == 1) {motBC = 'b';} // mother count as b
	else if(mot == 0) {;} // uninformative
	else {;}
    }
    // the kid is 11
    else if(kid == 2){
	// for father
	if(fat == 0){
	    if(allowDenovo) fatBC = 'c';
	    else fatBC = 'n';
	    //mutationHappened = true;
	    denovo = true;
	} // de novo mutation
	else if(fat == 1) {fatBC = 'c';} // count as c
	else if(fat == 2) {;} // uninformative
	else {;}

	// for mother
	if(mot == 0){
	    if(allowDenovo) motBC = 'c';
	    else motBC = 'n';
	    //mutationHappened = true;
	    denovo = true;
	} // de novo mutation
	else if(mot == 1) {motBC = 'c';} // count as c
	else if(mot == 2) {;} // uninformative
	else {;}
		
	// no double mutation
	if(fat == 0 && mot == 0) {motBC = 'n'; fatBC = 'n'; denovo = false;}
    }
    // the kid is 01, there ate 3*3 = 9 situations, with 6 unique
    else if(kid == 1){
	if(fat == 0 && mot == 0){
	    if(allowDenovo){
		if (gsl_rng_uniform(gslr) > 0.5) fatBC = 'c';
		else motBC = 'c';
	    }
	    //mutationHappened = true;
	    denovo = true;
        } // de novo mutation
	
	else if(fat == 0 && mot == 1) motBC = 'c';
	else if(fat == 1 && mot == 0) fatBC = 'c'; // one count as c
	
	else if((fat == 0 && mot == 2)||(fat == 2 && mot == 0)) {;} // both are uninformative
	else if(fat == 1 && mot == 1) {
	    if (gsl_rng_uniform(gslr) > 0.5) {fatBC = 'c'; motBC = 'b';}
	    else  {fatBC = 'b'; motBC = 'c';}
	    } // one parent count as c, the other count as b; because it's unphased, we don't care which parent is b
	
	else if(fat == 1 && mot == 2) fatBC = 'b';
	else if(fat == 2 && mot == 1) motBC = 'b'; // one count as b
	
	else if(fat == 2 && mot == 2) {mutationHappened = true;} // reverse mutation
	else {;}
    }
    else {;}
    
    vectorC bc(3);
    bc[0] = fatBC; bc[1] = motBC;
    if(denovo) bc[2] = 'y';
    else bc[2] = 'n';
    
    //there is mutation happened
    vectorC mutated(2, 'm');
    if(mutationHappened) return mutated;
    else if((!allowDenovo) && denovo) return mutated;
    else return bc;
}


vectorF unphasedTDT::m_untransmitted(double fat, double mot, double kid)
{
    RNG rng;
    gsl_rng* gslr = rng.get();	
    //random guess one: in case genotype error, or missing
    vectorF UNtrans(2,-9);
    if(fat == 2) UNtrans[0] = 1;
    else if(fat == 1) {if (gsl_rng_uniform(gslr) > 0.5) UNtrans[0] = 1; else UNtrans[0] = 0;}
    else if(fat == 0) UNtrans[0] = 0;
    else UNtrans[0] = -9;

    if(mot == 2) UNtrans[1] = 1;
    else if(mot == 1) {if (gsl_rng_uniform(gslr) > 0.5) UNtrans[1] = 1; else UNtrans[1] = 0;}
    else if(mot == 0) UNtrans[1] = 0;
    else UNtrans[1] = -9;
    
    if(kid == -9 || (fat == -9 && mot == -9)) return UNtrans;
    else{
	if(fat == -9)
	{
	    UNtrans[0] = -9;
	    if(mot == 0){
		if(kid == 0) {UNtrans[1] = 0;}
		else if(kid == 1) {UNtrans[1] = 0;}
		else if(kid == 2) {UNtrans[1] = -9;}
		else {;}
	    }
	    else if(mot == 1){
		if(kid == 0) {UNtrans[1] = 1;}
		else if(kid == 1) {UNtrans[1] = 0;}
		else if(kid == 2) {UNtrans[1] = 0;}
		else {;}
	    }
	    else if(mot == 2){
		if(kid == 0) {UNtrans[1] = -9;}
		else if(kid == 1) {UNtrans[1] = 1;}
		else if(kid == 2) {UNtrans[1] = 1;}
		else {;}
	    }
	    else {;}
	}
	else if(mot == -9)
	{
	    UNtrans[1] = -9;
	    if(fat == 0){
		if(kid == 0) {UNtrans[0] = 0;}
		else if(kid == 1) {UNtrans[0] = 0;}
		else if(kid == 2) {UNtrans[0] = -9;}
		else {;}
	    }
	    else if(fat == 1){
		if(kid == 0) {UNtrans[0] = 1;}
		else if(kid == 1) {UNtrans[0] = 0;}
		else if(kid == 2) {UNtrans[0] = 0;}
		else {;}
	    }
	    else if(fat == 2){
		if(kid == 0) {UNtrans[0] = -9;}
		else if(kid == 1) {UNtrans[0] = 1;}
		else if(kid == 2) {UNtrans[0] = 1;}
		else {;}
	    }
	    else {;}
	}
	// no missing
	else {
	    if(fat == 0){
		if(mot == 0){
		    if(kid == 0) {UNtrans[0] = 0; UNtrans[1] = 0;}
		    else if(kid == 1) {;}
		    else if(kid == 2) {;}
		    else {;}
		}
		else if(mot == 1){
		    if(kid == 0) {UNtrans[0] = 0; UNtrans[1] = 1;}
		    else if(kid == 1) {UNtrans[0] = 0; UNtrans[1] = 0;}
		    else if(kid == 2) {;}
		    else {;}
		}
		else if(mot == 2){
		    if(kid == 0) {;}
		    else if(kid == 1) {UNtrans[0] = 0; UNtrans[1] = 1;}
		    else if(kid == 2) {;}
		    else {;}
		}
		else {;}
	    }
	    else if(fat == 2){
		if(mot == 0){
		    if(kid == 0) {;}
		    else if(kid == 1) {UNtrans[0] = 1; UNtrans[1] = 0;}
		    else if(kid == 2) {;}
		    else {;}
		}
		else if(mot == 1){
		    if(kid == 0) {;}
		    else if(kid == 1) {UNtrans[0] = 1; UNtrans[1] = 1;}
		    else if(kid == 2) {UNtrans[0] = 1; UNtrans[1] = 0;}
		    else {;}
		}
		else if(mot == 2){
		    if(kid == 0) {;}
		    else if(kid == 1) {;}
		    else if(kid == 2) {UNtrans[0] = 1; UNtrans[1] = 1;}
		    else {;}
		}
		else {;}
	    }
	    else if(fat == 1){
		if(mot == 0){
		    if(kid == 0) {UNtrans[0] = 1; UNtrans[1] = 0;}
		    else if(kid == 1) {UNtrans[0] = 0; UNtrans[1] = 0;}
		    else if(kid == 2) {;}
		    else {;}
		}
		else if(mot == 1){
		    if(kid == 0) {UNtrans[0] = 1; UNtrans[1] = 1;}
		    else if(kid == 1) {UNtrans[0] = 1; UNtrans[1] = 0;}
		    else if(kid == 2) {UNtrans[0] = 0; UNtrans[1] = 0;}
		    else {;}
		}
		else if(mot == 2){
		    if(kid == 0) {;}
		    else if(kid == 1) {UNtrans[0] = 1; UNtrans[1] = 1;}
		    else if(kid == 2) {UNtrans[0] = 0; UNtrans[1] = 1;}
		    else {;}
		}
		else {;}
	    }
	    else {;}
	}
	return UNtrans;
    }
}


vectorF unphasedTDT::getWSSweight()
{
    return __weight;
}

// get dataset information
UINT unphasedTDT::getFamNum()
{
    return __trioGenos.size() - __wrongFam;
}

UINT unphasedTDT::getXvariNum()
{
    return __Xvariants;
}

UINT unphasedTDT::getVariNum()
{
    return __varNumToBeAnalyzed;
}

vectorL unphasedTDT::getFlip()
{
    return __flip;
}

vectorL unphasedTDT::getBeAnalyzed()
{
    return __lowMissing;
}

vectorUI unphasedTDT::getErrFam()
{
    return __discardFam;
}

vectorUI unphasedTDT::getDenovo()
{
    return __denovo;
}

vectorF unphasedTDT::getMissRatio()
{
    return __MissingRatio;
}

vectorF unphasedTDT::getTdtCountB()
{
    vectorF sum(__varNum, 0);
    for(UINT i=0; i<__TdtCountB.size(); i++)
    {
	for(UINT j=0; j<__TdtCountB[i].size(); j++) sum[j] += __TdtCountB[i][j]; 
    }
    return sum;
}

vectorF unphasedTDT::getTdtCountC()
{
    vectorF sum(__varNum, 0);
    for(UINT i=0; i<__TdtCountC.size(); i++)
    {
	for(UINT j=0; j<__TdtCountC[i].size(); j++) sum[j] += __TdtCountC[i][j];
    }
    return sum;
}

vectorF unphasedTDT::getSingleP() 
{
    vectorF Bcount = getTdtCountB();
    vectorF Ccount = getTdtCountC();
    vectorF singleP;
    
    for(UINT i=0; i<Bcount.size(); i++){
	if((Ccount[i] + Bcount[i]) == 0) singleP.push_back(-9);
	else singleP.push_back(gsl_cdf_gaussian_Q((Ccount[i]-Bcount[i])/sqrt(Ccount[i]+Bcount[i]), 1));
    }
    return singleP;
}

vector2F unphasedTDT::unphasedSST(double mafcutoff)
{
    vectorF singlePval = getSingleP();
    vector2F res;
    // SST test can't be used with VT, because VT will disturb popMafs 
    for(UINT i=0; i<singlePval.size(); i++){
	if(__lowMissing[i] && __samMafs[i] >= mafcutoff) {
	    vectorF oneSite;
	    oneSite.push_back(i+1); oneSite.push_back(singlePval[i]);
	    res.push_back(oneSite);
	}
    }
    return res;
}

// will return MZ value with upper boundary and no boundary
vectorF unphasedTDT::TdtMZ(double mafLower, double mafUpper, bool returnChi) const
{
    double tdtB(0), tdtC(0), tdtBall(0), tdtCall(0);
    
    if(DEBUG) for (UINT i = 0; i < __TdtCountB.size(); i++) std::cout << "tdtCountB: " << __TdtCountB[i] << std::endl;
    if(DEBUG) for (UINT i = 0; i < __TdtCountC.size(); i++) std::cout << "tdtCountC: " << __TdtCountC[i] << std::endl;

    for (UINT i = 0; i < __TdtCountB.size(); i++){
	for (UINT j = 0; j < __TdtCountB[i].size(); j++){
	    tdtBall += __TdtCountB[i][j]; tdtCall += __TdtCountC[i][j];
	    // ignore if out of [lower, upper]
	    if(__popMafs[j] < mafUpper && __popMafs[j] >= mafLower) {tdtB += __TdtCountB[i][j]; tdtC += __TdtCountC[i][j];}
	}
    }

    vectorF MZresult, MZchisqu;
    
    if(tdtB==0 && tdtC==0) {MZresult.push_back(std::numeric_limits<double>::infinity()); MZchisqu.push_back(std::numeric_limits<double>::infinity());}
    else{
	double MZone = gsl_cdf_gaussian_Q((tdtC-tdtB)/sqrt(tdtB+tdtC), 1);
	MZchisqu.push_back((tdtC-tdtB)/sqrt(tdtB+tdtC));
	MZresult.push_back(MZone);
    }
    
    if(tdtBall==0 && tdtCall==0) {MZresult.push_back(std::numeric_limits<double>::infinity()); MZchisqu.push_back(std::numeric_limits<double>::infinity());}
    else{
	double MZoneAll = gsl_cdf_gaussian_Q((tdtCall-tdtBall)/sqrt(tdtBall+tdtCall), 1);
	MZchisqu.push_back((tdtCall-tdtBall)/sqrt(tdtBall+tdtCall));
	MZresult.push_back(MZoneAll);
    }
    
    //std::cout << tdtB << " " << tdtC << "\t" << (tdtC-tdtB)/sqrt(tdtB+tdtC) << std::endl;
    return (returnChi)? MZchisqu:MZresult;
}


// VT: return a vector including one-sided chi-square with and without boundary
vectorF unphasedTDT::m_tdtVT(double mafLower, double mafUpper) 
{
	vectorF SortedMAF = __popMafs;
	sort(SortedMAF.begin(), SortedMAF.end());
	SortedMAF.erase(unique(SortedMAF.begin(), SortedMAF.end()), SortedMAF.end());
	__sortedpopMafs = SortedMAF;

	__sortedChiSqu.assign(SortedMAF.size(), 0);

	double VTMaxMZone(-99), VTMaxMZall(-99);

	//if(VTDEBUG) std::cout << "VT: " << SortedMAF << std::endl;

	for(UINT point = 0; point < SortedMAF.size(); point++)
	{
		// keep threshold within defined region
		//if(SortedMAF[point] < mafLower || SortedMAF[point] > mafUpper) continue;
		double TdtB(0), TdtC(0), TdtBall(0), TdtCall(0);

		// count tdtB: only record events on sites that maf less than threshold
		// same as tdtC
		for (UINT i = 0; i < __TdtCountB.size(); i++){    
			for (UINT j = 0; j < __TdtCountB[i].size(); j++){
				if(__popMafs[j] <= SortedMAF[point]){
					TdtBall += __TdtCountB[i][j]; TdtCall += __TdtCountC[i][j];
					if(__popMafs[j] <= mafUpper && __popMafs[j] >= mafLower) {TdtB += __TdtCountB[i][j]; TdtC += __TdtCountC[i][j];}
				}
				else {;}
			}
		}

		double CQone(0), CQoneAll(0);

		if(TdtB==0 && TdtC==0) {CQone = -99;}
		else {CQone = (TdtC-TdtB)/sqrt(TdtB+TdtC);}

		if(TdtBall==0 && TdtCall==0) {CQoneAll = -99;}
		else {CQoneAll = (TdtCall-TdtBall)/sqrt(TdtBall+TdtCall);}

		__sortedChiSqu[point] = CQoneAll;

		if (CQone > VTMaxMZone) VTMaxMZone = CQone;
		if (CQoneAll > VTMaxMZall) VTMaxMZall = CQoneAll;
		//std::cout <<  "old: " << SortedMAF[point] << "->" << VTMaxMZone << std::endl;
		//for (UINT j = 0; j < __TdtCountB[0].size(); j++){
			//if(__popMafs[j] <= SortedMAF[point]){
				//if(__popMafs[j] <= mafUpper && __popMafs[j] >= mafLower) {std::cout << "1\t";}
				//else {std::cout << "0\t";}
			//}
			//else {std::cout << "0\t";}
		//}
		//std::cout << std::endl;
	}

	//if(VTDEBUG) std::cout << "VTMax: " << VTMaxMZone << " " << VTMaxMZall << std::endl;
	//std::cout << __MissingRatio << std::endl;
	//std::cout << m_colTotal(__TdtCountC) << std::endl;
	////std::cout << m_colTotal(__TdtCountB) << std::endl;
	//std::cout << __popMafs << std::endl;
	//std::cout << __lowMissing << std::endl;
	//std::cout << std::endl;
	vectorF VTresult;
	VTresult.push_back(VTMaxMZone); VTresult.push_back(VTMaxMZall);
	return VTresult;
}


vectorF unphasedTDT::m_tdtWSS(double mafLower, double mafUpper)
{
    //std::cout << "WSS called: " << __NonTransChrom[0].size() << std::endl;
    //calculate the weight
    vectorI NonZero(__NonTransChrom[0].size(),0), NonMiss(__NonTransChrom[0].size(),0);
    for (UINT i=0; i<__NonTransChrom.size(); i++){
	//std::cout << __NonTransChrom[i] << std::endl;
        for (UINT j=0; j < __NonTransChrom[i].size(); j++){
            if(__NonTransChrom[i][j] > 0) NonZero[j] += __NonTransChrom[i][j];
            if(__NonTransChrom[i][j] >= 0) NonMiss[j] ++;
        }
    }
    //std::cout << "======" << std::endl;
    //
    vectorI NotMissAll(__NonTransChrom[0].size(),0);
    for (UINT i=0; i<__trioGenos.size(); i++){
        for (UINT j=0; j < __trioGenos[i].size(); j++){
	    for(UINT h=0; h < __trioGenos[i][j].size(); h++){
		if(__trioGenos[i][j][h] >= 0) NotMissAll[h] ++;
	    }
        }
    }
    
    vectorF weight;
    for (UINT i = 0; i < NonZero.size(); i++) __weight.push_back( m_WssWeightCal(NotMissAll[i], NonMiss[i], NonZero[i]));

    //std::cout << "-------" << std::endl;
    //std::cout << __trioGenos.size() << "|" << __NonTransChrom.size() << std::endl;
    //std::cout << "NoMAll: " << NotMissAll << "\nNoMiss: " << NonMiss << std::endl << "NoZero: " << NonZero << std::endl;
    //std::cout << "weight: " << __weight << std::endl;
    
    double tdtB(0), tdtC(0), tdtBall(0), tdtCall(0);
    for (UINT i = 0; i < __TdtCountB.size(); i++){
        for (UINT j = 0; j < __TdtCountB[i].size(); j++){
	    tdtBall += __weight[j]*__TdtCountB[i][j]; tdtCall += __weight[j]*__TdtCountC[i][j];
	    // ignore if out of [0, upper]
	    if(__popMafs[j] < mafUpper && __popMafs[j] >= mafLower) {tdtB += __weight[j]*__TdtCountB[i][j]; tdtC += __weight[j]*__TdtCountC[i][j];}
        }
    }
    vectorF res;
    res.push_back((tdtC+tdtB == 0)?std::numeric_limits<double>::infinity():(tdtC-tdtB)/sqrt(tdtB+tdtC));
    res.push_back((tdtCall+tdtBall == 0)?std::numeric_limits<double>::infinity():(tdtCall-tdtBall)/sqrt(tdtBall+tdtCall));

    //std::cout << tdtC << "\t" << tdtB << "\t" << (tdtC-tdtB)/sqrt(tdtB+tdtC) << std::endl;
    return res;
}


bool unphasedTDT::m_zeroFamily(double P, double M, double C)
{
    if(__SkipMiss){
	if(P == -9.0 || M == -9.0 || C == -9.0) return true;
	else if(P == M && P == C) return true;
	else return false;
    }
    else{
	if(P == -9.0 && M == -9.0) return true;
	else if(M == -9.0) return true;
	else if(P == M && P == C) return true;
	else return false;
    }
}

// trim those uninformative sites, based on missing and whether informative, to speed up permutation
void unphasedTDT::TrimUninformSites()
{
    //Before permutation, delete some variant site which are uninformative in data set
    __unInformSites.assign(__varNum, true);
    for(UINT s=0; s<__varNum; s++){
        for(UINT f=0; f<__trioGenos.size(); f++){
	    // keep this site if any trio is informative
            if(!m_zeroFamily(__trioGenos[f][0][s], __trioGenos[f][1][s], __trioGenos[f][2][s])) {__unInformSites[s]=false; break;}
        }
    }
    
    vectorL shouldTrim = __unInformSites;
    //combine with __lowMissing vector
    for(UINT i=0; i<shouldTrim.size(); i++){
	if(!__lowMissing[i]) shouldTrim[i]=true;
    }
    
    //check if there is uninformative site
    bool NeedTrim(false);
    for(UINT i=0; i<__varNum; i++) NeedTrim = (NeedTrim||shouldTrim[i]);
    if(NeedTrim){
        vector3F NewData;
        for(UINT f=0; f<__trioGenos.size(); f++){
            vectorF P, M, C;
            for(UINT s=0; s<__varNum; s++){
                if(!shouldTrim[s]) {P.push_back(__trioGenos[f][0][s]); M.push_back(__trioGenos[f][1][s]); C.push_back(__trioGenos[f][2][s]);}
            }
            vector2F NewFamily;
            NewFamily.push_back(P); NewFamily.push_back(M); NewFamily.push_back(C);
            NewData.push_back(NewFamily);
        }
        __trioGenos = NewData;
        
        vectorF NewPopMaf, NewSamMaf;
        for(UINT s=0; s<__varNum; s++){
            if(!shouldTrim[s]) {NewPopMaf.push_back(__popMafs[s]); NewSamMaf.push_back(__samMafs[s]);}
        }
        __popMafs = NewPopMaf;
	__samMafs = NewSamMaf;
    }
    return;
}

// return the p values with and without boundary
vectorF unphasedTDT::TdtPermutVT(double mafLower, double mafUpper, UINT adaptive, UINT PermutateTimes, double alpha, bool shuffleThree, std::string method)
{
    // check whether can use VT: must have same variant site number in all trios -- __varNum
    for(UINT i=0; i<__trioGenos.size(); i++){
	// after FamilyCheck, the genotypes martix within a trio should be square matrix, so just check the father's variant site number
	if(__trioGenos[i][0].size() != __varNum) {std::cout << "Error: there are different variant sites number in different trio, can't use unphased::TdtVT method, please check your data" << std::endl; exit(-1);}
    }

    // the VT static for original dataset
    vectorF oriVT = m_tdtVT(mafLower, mafUpper);
    if(oriVT[0] == -9.0 || oriVT[1] == -9.0){
	vectorF res(2, std::numeric_limits<double>::infinity());
	return res;
    }
    // trim unused sites, to speed up
    TrimUninformSites();
    
    UINT permcount(0), permcountAll(0);
    double pvalue(9.0), pvalueAll(9.0);
    __permTimes = 0;
    
    for (UINT i=1; i <= PermutateTimes; i++)
    {
        vector3F GenoShuffDate;
	if(method == "GenoShuffle") GenoShuffDate = m_GenoShuffle(shuffleThree);
	else if(method == "HapoShuffle") GenoShuffDate = m_HapoShuffle();
	else {std::cout << "Invalid shuffle method for VT" << std::endl; exit(-1);}
	
	unphasedTDT GenoPermut(GenoShuffDate, __popMafs, __samMafs, __chr, __SkipMiss);
	vectorF statistic = GenoPermut.m_tdtVT(mafLower, mafUpper);
	for(UINT t=0; t<statistic.size(); t++) {if(statistic[t]==-99) statistic[t] = 0;}
	
	RNG rng;
        gsl_rng* gslr = rng.get();
	double ok = gsl_rng_uniform(gslr); 
	
	//with boundary
	if (statistic[0] > oriVT[0]) { permcount++;}
        else if (statistic[0] == oriVT[0]) {if (ok > 0.5) permcount++;}
	else {;}
	if(pvalue > 1.0) pvalue = m_check(permcount, i, adaptive, alpha);
	
	//without boundary
	if (statistic[1] >= oriVT[1]) { permcountAll++;}
        else if (statistic[1] == oriVT[1]) {if (ok > 0.5) permcountAll++;}
	else {;}
	if(pvalueAll > 1.0) pvalueAll = m_check(permcountAll, i, adaptive, alpha);	
	
	__permTimes++;
        if (pvalue <= 1.0 && pvalueAll <= 1.0) {break;}
    }
  
    vectorF VtPermRes;
    VtPermRes.push_back((pvalue <= 1.0)? pvalue : (1.0 * permcount + 1) / (1.0 * __permTimes + 1));
    VtPermRes.push_back((pvalueAll <= 1.0)? pvalueAll : (1.0 * permcountAll + 1) / (1.0 * __permTimes + 1));
        
    return VtPermRes;
}

// return the p values with and without boundary
vectorF unphasedTDT::TdtPermutMZ(double mafLower, double mafUpper, UINT adaptive, UINT PermutateTimes, double alpha, bool shuffleThree, std::string method)
{
    // check whether can use VT: must have same variant site number in all trios -- __varNum
    for(UINT i=0; i<__trioGenos.size(); i++){
	// after FamilyCheck, the genotypes martix within a trio should be square matrix, so just check the father's variant site number
	if(__trioGenos[i][0].size() != __varNum) {std::cout << "Error: there are different variant sites number in different trio, can't use unphased::unphasedMZ method, please check your data" << std::endl; exit(-1);}
    }
    // the VT static for original dataset
    vectorF oriMZ = TdtMZ(mafLower, mafUpper, true);
    if(oriMZ[0] == std::numeric_limits<double>::infinity() || oriMZ[0] == std::numeric_limits<double>::infinity()){
	vectorF res(2, std::numeric_limits<double>::infinity());
	return res;
    }
    TrimUninformSites();
    
    UINT permcount(0), permcountAll(0);
    double pvalue(9.0), pvalueAll(9.0);
    __permTimes = 0;
    
    for (UINT i=1; i <= PermutateTimes; i++)
    {
        vector3F GenoShuffDate;
	if(method == "GenoShuffle") GenoShuffDate = m_GenoShuffle(shuffleThree);
	else if(method == "HapoShuffle") GenoShuffDate = m_HapoShuffle();
	else {std::cout << "Invalid shuffle method for VT" << std::endl; exit(-1);}
	
	unphasedTDT GenoPermut(GenoShuffDate, __popMafs, __samMafs, __chr, __SkipMiss);
	vectorF statistic = GenoPermut.TdtMZ(mafLower, mafUpper, true);
	for(UINT t=0; t<statistic.size(); t++) {if(statistic[t]==std::numeric_limits<double>::infinity()) statistic[t] = 0;}
	
	//if(VTDEBUG) std::cout << "statistic: " << statistic << std::endl;
	RNG rng;
        gsl_rng* gslr = rng.get();
	double ok = gsl_rng_uniform(gslr);
	
	//with boundary
	if (statistic[0] > oriMZ[0]) { permcount++;}
        else if (statistic[0] == oriMZ[0]) {if (ok > 0.5) permcount++;}
	else {;}
	if(pvalue > 1.0) pvalue = m_check(permcount, i, adaptive, alpha);
	
	//without boundary
	if (statistic[1] >= oriMZ[1]) { permcountAll++;}
        else if (statistic[1] == oriMZ[1]) {if (ok > 0.5) permcountAll++;}
	else {;}
	if(pvalueAll > 1.0) pvalueAll = m_check(permcountAll, i, adaptive, alpha);	
	
	__permTimes++;
        if (pvalue <= 1.0 && pvalueAll <= 1.0) {break;}
    }
  
    vectorF MzPermRes;
    MzPermRes.push_back((pvalue <= 1.0)? pvalue : (1.0 * permcount + 1) / (1.0 * __permTimes + 1));
    MzPermRes.push_back((pvalueAll <= 1.0)? pvalueAll : (1.0 * permcountAll + 1) / (1.0 * __permTimes + 1));
        
    return MzPermRes;
}


// For WSS
vectorF unphasedTDT::TdtPermutWSS(double mafLower, double mafUpper, UINT adaptive, UINT PermutateTimes, double alpha, bool shuffleThree, std::string method)
{
    for(UINT i=0; i<__trioGenos.size(); i++){
	// after FamilyCheck, the genotypes martix within a trio should be square matrix, so just check the father's variant site number
	if(__trioGenos[i][0].size() != __varNum) {std::cout << "Error: there are different variant sites number in different trio, can't use unphased::TdtWSS method, please check your data" << std::endl; exit(-1);}
    }

    // the VT static for original dataset
    vectorF oriWSS = m_tdtWSS(mafLower, mafUpper);
    if(oriWSS[0] == std::numeric_limits<double>::infinity() || oriWSS[1] == std::numeric_limits<double>::infinity()){
	return oriWSS;
    }
    
    //std::cout << oriWSS[0] << "\t" << oriWSS[1] <<  std::endl << "--------" << std::endl;
    TrimUninformSites();
    UINT permcountWSS(0), permcountWSSall(0);
    double pvalueWSS(9.0), pvalueWSSall(9.0);
    __permTimes = 0;

    for (UINT i=1; i <= PermutateTimes; i++)
    {
        vector3F GenoShuffData;
	if(method == "GenoShuffle") GenoShuffData = m_GenoShuffle(shuffleThree);
	else if (method == "HapoShuffle") GenoShuffData = m_HapoShuffle();
	else {std::cout << "ERROR: invalid shuffle method" << std::endl; exit(-1);}
	
	unphasedTDT GenoPermut(GenoShuffData, __popMafs, __samMafs, __chr, __SkipMiss);
        vectorF wss = GenoPermut.m_tdtWSS(mafLower, mafUpper);
	//std::cout << wss[0] << "\t" << wss[1] <<  std::endl;
	for(UINT t=0; t<wss.size(); t++) {if(wss[t]==std::numeric_limits<double>::infinity()) wss[t] = 0;}
	
	//if(VTDEBUG) std::cout << "statistic: " << statistic << std::endl;
	RNG rng;
        gsl_rng* gslr = rng.get();

	//with boundary
	if (wss[0] > oriWSS[0]) { permcountWSS++;}
        else if (wss[0] == oriWSS[0]) {if (gsl_rng_uniform(gslr) > 0.5) permcountWSS++;}
	else {;}
	if(pvalueWSS > 1.0) pvalueWSS = m_check(permcountWSS, i, adaptive, alpha);

	if (wss[1] > oriWSS[1]) { permcountWSSall++;}
        else if (wss[1] == oriWSS[1]) {if (gsl_rng_uniform(gslr) > 0.5) permcountWSSall++;}
	else {;}
	if(pvalueWSSall > 1.0) pvalueWSSall = m_check(permcountWSSall, i, adaptive, alpha);

	__permTimes++;
        if (pvalueWSS <= 1.0 && pvalueWSSall <= 1.0) {break;}
    }

    vectorF WssPermRes;
    WssPermRes.push_back((pvalueWSS <= 1.0)? pvalueWSS : (1.0 * permcountWSS + 1) / (1.0 * __permTimes + 1));
    WssPermRes.push_back((pvalueWSSall <= 1.0)? pvalueWSSall : (1.0 * permcountWSSall + 1) / (1.0 * __permTimes + 1));

    return WssPermRes;
}


vectorF unphasedTDT::TdtPermut(double mafLower, double mafUpper, UINT adaptive, UINT PermutateTimes, double alpha, bool shuffleThree, std::string method)
{
    permut = true;
    for(UINT i=0; i<__trioGenos.size(); i++){
	// after FamilyCheck, the genotypes martix within a trio should be square matrix, so just check the father's variant site number
	if(__trioGenos[i][0].size() != __varNum) {std::cout << "Error: there are different variant sites number in different trio, can't use unphased::TdtVT method, please check your data" << std::endl; exit(-1);}
    }

    // the VT static for original dataset
    double oriWSS = m_tdtWSS(mafLower, mafUpper)[0];
    double oriVTMZ = m_tdtVT(mafLower, mafUpper)[0];
    //std::cout << oriWSS << "\t" << oriVTMZ << "\t" <<  std::endl;
    if(oriWSS == std::numeric_limits<double>::infinity() || oriVTMZ == -99)
    {
	vectorF res(2, std::numeric_limits<double>::infinity());
	return res;
    }
    TrimUninformSites();
    
    //std::cout << "_transmitted_ev: " << getTdtCountC() << std::endl;
    //std::cout << "_untransmitted: " << getTdtCountB() << std::endl;
    //std::cout << oriWSS << "\t" << oriVTMZ << std::endl;
    UINT permcountWSS(0), permcountVTMZ(0);
    double pvalueWSS(9.0), pvalueVTMZ(9.0);
    __permTimes = 0;

    for (UINT i=1; i <= PermutateTimes; i++)
    {
        vector3F GenoShuffData;
	if(method == "GenoShuffle") GenoShuffData = m_GenoShuffle(shuffleThree);
	else if (method == "HapoShuffle") GenoShuffData = m_HapoShuffle();
	else {std::cout << "ERROR: invalid shuffle method" << std::endl; exit(-1);}
	
	unphasedTDT GenoPermut(GenoShuffData, __popMafs, __samMafs, __chr, __SkipMiss);
	double wss = GenoPermut.m_tdtWSS(mafLower, mafUpper)[0];
	if(wss==std::numeric_limits<double>::infinity()) wss = 0;
	
	double vtmz = GenoPermut.m_tdtVT(mafLower, mafUpper)[0];
	if(vtmz==-99) vtmz = 0;
	
	//std::cout << "-----" << std::endl;
	//std::cout << "_transmitted_ev: " << GenoPermut.getTdtCountC() << std::endl;
	//std::cout << "_untransmitted: " << GenoPermut.getTdtCountB() << std::endl;
	//std::cout << wss << "\t" << vtmz << std::endl;
	
	RNG rng;
        gsl_rng* gslr = rng.get();
	//std::cout << wss << "\t";
	//WSS
	if (wss> oriWSS) { permcountWSS++;}
        else if (wss == oriWSS) {if (gsl_rng_uniform(gslr) > 0.5) permcountWSS++;}
	else {;}
	if(pvalueWSS > 1.0) pvalueWSS = m_check(permcountWSS, i, adaptive, alpha);

	//VT-MZ
	if (vtmz> oriVTMZ) { permcountVTMZ++;}
        else if (vtmz == oriVTMZ) {if (gsl_rng_uniform(gslr) > 0.5) permcountVTMZ++;}
	else {;}
	if(pvalueVTMZ > 1.0) pvalueVTMZ = m_check(permcountVTMZ, i, adaptive, alpha);

	__permTimes++;
        if (pvalueWSS <= 1.0 && pvalueVTMZ <= 1.0) {break;}
    }
    //std::cout << std::endl;
    vectorF PermRes;
    PermRes.push_back((pvalueWSS <= 1.0)? pvalueWSS : (1.0 * permcountWSS + 1) / (1.0 * __permTimes + 1));
    PermRes.push_back((pvalueVTMZ <= 1.0)? pvalueVTMZ : (1.0 * permcountVTMZ + 1) / (1.0 * __permTimes + 1));

    return PermRes;
}


// T is the number of SNPs or genes
double unphasedTDT::m_check(UINT permcount, UINT iPermut, UINT adaptive, double alpha)
{
    if(iPermut % adaptive != 0 || iPermut == 0) {
        return 9.0;
    }
    double pval = permcount*1.0/iPermut*1.0;
    double sigma = sqrt(pval*(1-pval)/iPermut);
    
    double beta = 0.05;
    double gs = gsl_cdf_gaussian_Pinv(1.0-beta/2.0, sigma);
    //return (pval > alpha + 6*sigma)? pval:9.0;
    return (pval - gs > alpha)? pval:9.0;
}

vectorF unphasedTDT::getSortedMafs()
{
    return __sortedpopMafs;
}

vectorF unphasedTDT::getSortedCQ()
{
    return __sortedChiSqu;
}

UINT unphasedTDT::getPermTimes()
{
    return __permTimes;
}

//Shuffle the genotupe, that is shuffle single base and each site
vector3F unphasedTDT::m_GenoShuffle(bool shuffleThree) const
{
    vector3F NewData;
    UINT variantNum = __trioGenos[0][0].size();
    RNG rng;
    gsl_rng* gslr = rng.get();
    
    for (UINT OneFam=0; OneFam < __trioGenos.size(); OneFam++)
    {
        //for(UINT i=0; i<__trioGenos[OneFam].size(); i++) std::cout << __trioGenos[OneFam][i] << std::endl;
	//std::cout << "----" << std::endl;
	vector2F NewFamily;
        NewFamily.resize(3);
	
        for (UINT EachSite = 0; EachSite < variantNum; EachSite++)
        {
	    double fat(__trioGenos[OneFam][0][EachSite]), mot(__trioGenos[OneFam][1][EachSite]), kid(__trioGenos[OneFam][2][EachSite]);
	    
            vectorUI SiteGeno;
	    //deal with unknown genotype coding, eg 3 4...
	    if((fat != 0 && fat != 1 && fat != 2 && fat != -9) || (mot != 0 && mot != 1 && mot != 2 && mot != -9) || (kid != 0 && kid != 1 && kid != 2 && kid != -9)){
		NewFamily[0].push_back(-9); NewFamily[1].push_back(-9); NewFamily[2].push_back(-9);
		std::cout << "invalid genotype coding found in Genotype Permutation" << std::endl; exit(-1);
	    }
	    //uninformative site, keep it
	    else if(kid == -9 || (fat == -9 && mot == -9)) {
		NewFamily[0].push_back(fat); NewFamily[1].push_back(mot); NewFamily[2].push_back(kid);
	    }
	    else if(fat == -9 || mot == -9){
		// missing father's genotype
		if(fat == -9){
		    double guessedFatHapo(0);
		    if(mot == 0 && kid == 0) guessedFatHapo = 0;
		    else if(mot == 0 && kid == 1) guessedFatHapo = 1;
		    else if(mot == 0 && kid == 2) guessedFatHapo = -1; // 0X-9 -> 2
		    else if(mot == 1 && kid == 0) guessedFatHapo = 0;
		    else if(mot == 1 && kid == 1) guessedFatHapo = (gsl_rng_uniform(gslr) > 0.5)?1:0;
		    else if(mot == 1 && kid == 2) guessedFatHapo = 1;
		    else if(mot == 2 && kid == 0) guessedFatHapo = -1; // 2X-9 ->0
		    else if(mot == 2 && kid == 1) guessedFatHapo = 0;
		    else if(mot == 2 && kid == 2) guessedFatHapo = 1;
		    else {;}
			
		    // genotype error
		    if(guessedFatHapo < 0) {NewFamily[0].push_back(fat); NewFamily[1].push_back(mot); NewFamily[2].push_back(kid);}
		    else{
			if(shuffleThree){
			    // father still missing
			    NewFamily[0].push_back(-9);
			    
			    vectorF binaryGeno(3,0);
			    // how many variants they have totally
			    double tot = guessedFatHapo + mot;
			    for (UINT i = 0; i < tot; i++) binaryGeno[i] = 1;
			    //shuffle three hapotypes
			    random_shuffle (binaryGeno.begin(), binaryGeno.begin()+3);
		
			    // mother
			    NewFamily[1].push_back(binaryGeno[0] + binaryGeno[1]);
			    // kid
			    NewFamily[2].push_back(binaryGeno[rand()%2] + binaryGeno[2]);
			}
			// only shuffle mother
			else
			{
			    // father still missing
			    NewFamily[0].push_back(-9);
			    
			    vectorF binaryGeno(2,0);
			    // how many variants mother have totally
			    for (UINT i = 0; i <  mot; i++) binaryGeno[i] = 1;
			    //shuffle two hapotypes
			    random_shuffle (binaryGeno.begin(), binaryGeno.begin()+2);
		
			    // mother keep the same
			    NewFamily[1].push_back(mot);
			    // kid
			    NewFamily[2].push_back(binaryGeno[0] + guessedFatHapo);
			}
		    }
		}
		else if(mot == -9){
		    double guessedMotHapo(0);
		    if(fat == 0 && kid == 0) guessedMotHapo = 0;
		    else if(fat == 0 && kid == 1) guessedMotHapo = 1;
		    else if(fat == 0 && kid == 2) guessedMotHapo = -1; // 0X-9 -> 2
		    else if(fat == 1 && kid == 0) guessedMotHapo = 0;
		    else if(fat == 1 && kid == 1) guessedMotHapo = (gsl_rng_uniform(gslr) > 0.5)?1:0;
		    else if(fat == 1 && kid == 2) guessedMotHapo = 1;
		    else if(fat == 2 && kid == 0) guessedMotHapo = -1; // 2X-9 ->0
		    else if(fat == 2 && kid == 1) guessedMotHapo = 0;
		    else if(fat == 2 && kid == 2) guessedMotHapo = 1;
		    else {;}
			
		    // genotype error
		    if(guessedMotHapo < 0) {NewFamily[0].push_back(fat); NewFamily[1].push_back(mot); NewFamily[2].push_back(kid);}
		    else{
			// shuffle father and half mother
			if(shuffleThree)
			{
			    vectorF binaryGeno(3,0);
			    // how many variants they have totally
			    double tot = guessedMotHapo + fat;
			    for (UINT i = 0; i < tot; i++) binaryGeno[i] = 1;
			    //shuffle three hapotypes
			    random_shuffle (binaryGeno.begin(), binaryGeno.begin()+3);
		
			    // father
			    NewFamily[0].push_back(binaryGeno[0] + binaryGeno[1]);
			    // mother still missing
			    NewFamily[1].push_back(-9);
			    // kid
			    NewFamily[2].push_back(binaryGeno[rand()%2] + binaryGeno[2]);
			}
			// only shuffle father
			else
			{
			    vectorF binaryGeno(2,0);
			    // how many variants mother have totally
			    for (UINT i = 0; i <  fat; i++) binaryGeno[i] = 1;
			    //shuffle two hapotypes
			    random_shuffle (binaryGeno.begin(), binaryGeno.begin()+2);
		
			    // father keep the same 
			    NewFamily[0].push_back(fat);
			    // mother still missing
			    NewFamily[1].push_back(-9);
			    // kid
			    NewFamily[2].push_back(binaryGeno[0] + guessedMotHapo);
			}
		    }
		}
		else {std::cout << "Something goes wrong on line 813" << std::endl; exit(-1);}
	    }
	    //no missing
	    else
	    {
		//
		vectorF binaryGeno(4,0);
		// how many variants the parents have totally
		double tot = __trioGenos[OneFam][0][EachSite] + __trioGenos[OneFam][1][EachSite];
		for (UINT i = 0; i < tot; i++) binaryGeno[i] = 1;
		//shuffle parents' genotype
		random_shuffle (binaryGeno.begin(), binaryGeno.begin()+4);
		
		binaryGeno.push_back(binaryGeno[rand()%2]);
		binaryGeno.push_back(binaryGeno[rand()%2+2]);
		
		for (UINT i = 0; i < 3; i++) NewFamily[i].push_back(binaryGeno[2*i] + binaryGeno[2*i+1]);
	    }
	}
	//for(UINT i=0; i<NewFamily.size(); i++) std::cout << NewFamily[i] << std::endl;
	//std::cout << "=====" << std::endl;
        NewData.push_back(NewFamily);
    }
    return NewData;
}

// shuffle hapotype
//Shuffle the hapotype, that is shuffle the chromosome
vector3F unphasedTDT::m_HapoShuffle() const
{
    vector3F NewData;
    UINT VariantNum = __trioGenos[0][0].size();
    
    RNG rng;
    gsl_rng* gslr = rng.get();
	
    for (UINT OneFam = 0; OneFam < __trioGenos.size(); OneFam++)
    {
        //for(UINT i=0; i<__trioGenos[OneFam].size(); i++) std::cout << __trioGenos[OneFam][i] << std::endl;
	//std::cout << "#----#" << std::endl;
		
	//frist, record where missing happens, if the site is informtive
	//second, infer missing genotype in parents genotype
	//then. re-construct missing
        bool Missing(false), UninformFlag(false);
        vectorL FatMissing(VariantNum, false), MotMissing(VariantNum, false), Uninform(VariantNum, false);
	
	//reconstruct it as phased data
	vector2F NewParent; NewParent.resize(4);
	for(UINT i=0; i<VariantNum; i++){
	    //father
	    if(__trioGenos[OneFam][0][i] == 0) {NewParent[0].push_back(0); NewParent[1].push_back(0);}
	    else if(__trioGenos[OneFam][0][i] == 1) {
		if(gsl_rng_uniform(gslr) < 0.5) {NewParent[0].push_back(1); NewParent[1].push_back(0);}
		else {NewParent[0].push_back(0); NewParent[1].push_back(1);}
	    }
	    else if(__trioGenos[OneFam][0][i] == 2) {NewParent[0].push_back(1); NewParent[1].push_back(1);}
	    else {NewParent[0].push_back(-9); NewParent[1].push_back(-9);}
	    
	    //mother
	    if(__trioGenos[OneFam][1][i] == 0) {NewParent[2].push_back(0); NewParent[3].push_back(0);}
	    else if(__trioGenos[OneFam][1][i] == 1) {
		if(gsl_rng_uniform(gslr) < 0.5) {NewParent[2].push_back(1); NewParent[3].push_back(0);}
		else {NewParent[2].push_back(0); NewParent[3].push_back(1);}
	    }
	    else if(__trioGenos[OneFam][1][i] == 2) {NewParent[2].push_back(1); NewParent[3].push_back(1);}
	    else {NewParent[2].push_back(-9); NewParent[3].push_back(-9);}
	}
		
	for(UINT EachSite = 0; EachSite < VariantNum; EachSite++){
	    if(__trioGenos[OneFam][0][EachSite] == -9) {FatMissing[EachSite] = true; Missing = true;} // father missing
	    if(__trioGenos[OneFam][1][EachSite] == -9) {MotMissing[EachSite] = true; Missing = true;} // mothing missing
	    if((__trioGenos[OneFam][2][EachSite] == -9) || (FatMissing[EachSite] && MotMissing[EachSite])) {Uninform[EachSite] = true; UninformFlag = true;} // kid missing or both P missing
	    
	    double fat(__trioGenos[OneFam][0][EachSite]), mot(__trioGenos[OneFam][1][EachSite]), kid(__trioGenos[OneFam][2][EachSite]);
	    
	    //infer missing, if it's still informative
	    if(!Uninform[EachSite]){
		if(FatMissing[EachSite]){
		    double guessedFatHapo(0);
		    if(mot == 0 && kid == 0) guessedFatHapo = 0;
		    else if(mot == 0 && kid == 1) guessedFatHapo = 1;
		    else if(mot == 0 && kid == 2) guessedFatHapo = -1; // 0X-9 -> 2
		    else if(mot == 1 && kid == 0) guessedFatHapo = 0;
		    else if(mot == 1 && kid == 1) guessedFatHapo = (gsl_rng_uniform(gslr) > 0.5)?1:0;
		    else if(mot == 1 && kid == 2) guessedFatHapo = 1;
		    else if(mot == 2 && kid == 0) guessedFatHapo = -1; // 2X-9 ->0
		    else if(mot == 2 && kid == 1) guessedFatHapo = 0;
		    else if(mot == 2 && kid == 2) guessedFatHapo = 1;
		    else {;}
		    
		    if(guessedFatHapo < 0) {UninformFlag = true;}
		    else if(gsl_rng_uniform(gslr) < 0.5) {
			NewParent[0][EachSite] = guessedFatHapo; NewParent[1][EachSite] = (gsl_rng_uniform(gslr) <__popMafs[EachSite])?1:0;
			}
		    else{
			NewParent[1][EachSite] = guessedFatHapo; NewParent[0][EachSite] = (gsl_rng_uniform(gslr) <__popMafs[EachSite])?1:0;
		    }
		}
		if(MotMissing[EachSite]){
		    double guessedMotHapo(0);
		    if(fat == 0 && kid == 0) guessedMotHapo = 0;
		    else if(fat == 0 && kid == 1) guessedMotHapo = 1;
		    else if(fat == 0 && kid == 2) guessedMotHapo = -1; // 0X-9 -> 2
		    else if(fat == 1 && kid == 0) guessedMotHapo = 0;
		    else if(fat == 1 && kid == 1) guessedMotHapo = (gsl_rng_uniform(gslr) > 0.5)?1:0;
		    else if(fat == 1 && kid == 2) guessedMotHapo = 1;
		    else if(fat == 2 && kid == 0) guessedMotHapo = -1; // 2X-9 ->0
		    else if(fat == 2 && kid == 1) guessedMotHapo = 0;
		    else if(fat == 2 && kid == 2) guessedMotHapo = 1;
		    else {;}

		    if(guessedMotHapo < 0) {UninformFlag = true;}
		    else if(gsl_rng_uniform(gslr) < 0.5) {
			NewParent[2][EachSite] = guessedMotHapo; NewParent[3][EachSite] = (gsl_rng_uniform(gslr) <__popMafs[EachSite])?1:0;
			}
		    else{
			NewParent[3][EachSite] = guessedMotHapo; NewParent[2][EachSite] = (gsl_rng_uniform(gslr) <__popMafs[EachSite])?1:0;
		    }
		}
	    }
	}
	// shuffle
        random_shuffle (NewParent.begin(), NewParent.begin()+4);

        //The chromosome comes from father should be NewFamily[0 or 1]
        int ChildOneChr = rand()%2;
        int ChildTwoChr = rand()%2 +2;
        NewParent.push_back(NewParent[ChildOneChr]);
        NewParent.push_back(NewParent[ChildTwoChr]);
	
	//phased --> unphased
	vector2F NewFamily; NewFamily.resize(3);
	for(UINT i=0; i<NewParent.size()/2; i++){
	    for(UINT j=0; j<NewParent[2*i].size(); j++) NewFamily[i].push_back(NewParent[2*i][j]+NewParent[2*i+1][j]);
	}
	
        //re-construct missing
        if(Missing){
            for(UINT i=0; i<FatMissing.size(); i++){
                if(FatMissing[i]) {NewFamily[0][i] = -9.0;}
                if(MotMissing[i]) {NewFamily[1][i] = -9.0;}
            }
        }
        else {;}
	// reconstruct uninformative sites: for un-transmitted maf
	if(UninformFlag){
            for(UINT i=0; i<Uninform.size(); i++){
                if(Uninform[i]) {
		    for(UINT j = 0; j < 3; j++) NewFamily[j][i] = __trioGenos[OneFam][j][i];
		}
            }
	}
        //for(UINT i=0; i<NewFamily.size(); i++) std::cout << "100 1000 " << NewFamily[i] << std::endl;
	//std::cout << "=====" << std::endl;
	
        NewData.push_back(NewFamily);
    }

    return NewData;
}

// summarize data information
vector3F unphasedTDT::m_geTrio3D(const vector2F& genos, UINT famSize)
{
    vector3F trios;
    for(UINT i = 0; i < genos.size()/famSize; i++){
	// get one trio raw genotypes
        vector2F OneFamily;
        OneFamily.resize(famSize);
        copy(genos.begin() + i*famSize, genos.begin() + (i+1)*famSize, OneFamily.begin());
	
        if(m_FamilyCheck(OneFamily)){
            trios.push_back(OneFamily);
        }
        else{
	    __wrongFam++;
	}
    }
    return trios;
}

// return value: false means not square; true means normal
bool unphasedTDT::m_FamilyCheck(const vector2F& OneFam) const
{
    UINT VariantNum = OneFam[0].size();
    for(UINT i = 0; i < OneFam.size(); i++){
	// not a square
        if(OneFam[i].size() != VariantNum) {return false;}
    }
    return true;
}

