/*!\file unphased.cpp
   \brief Deal with phased tios

   Copyright: 2011 Zongxiao He @ Baylor College of Medicine
*/

#include "gw_utilities.h"
#include "phasedTDT.h"
#include "gsl/gsl_cdf.h"
extern bool DEBUG;
extern bool VTDEBUG;
extern bool Permut;
bool Permut = false;
bool check = true;

//phasedTDT::phasedTDT(const vector2F& genotypes, const vectorF& popMafs, const vectorF& samMafs, std::string& chrom, double siteMissRatio, bool SkipMiss)
phasedTDT::phasedTDT(const vector2F& genotypes, const vectorF& popMafs, const vectorF& samMafs, std::string& chrom, double siteMissRatio, bool SkipMiss)
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
    //__discardFam.assign(__varNum, 0);
    __denovo.assign(__varNum, 0);
    __flip.assign(__varNum, false);
    
    //std::cout << samMafs << std::endl;
    for(UINT i=0; i<__samMafs.size(); i++){
	if(__samMafs[i] > 0.5){
	    __flip[i] = true;
	    __samMafs[i] = 1-__samMafs[i];
	}
    }
    //std::cout << samMafs << std::endl;
    m_SiteMissRatio();
    if(siteMissRatio > 0){
	for(UINT i=0; i<__varNum; i++) {
	    if(__MissingRatio[i] >= siteMissRatio) {__lowMissing[i] = false; __Xvariants++;}
	}
    }
    __varNumToBeAnalyzed = __varNum - __Xvariants;
    __trioGenos = m_geTrio3D(__twoDGenos, 6);
    //if(__SkipMiss) m_markMiss();
    
    m_TdtCountForAllFamilies();
    // erase the sites that shouldn't be analyzed, because missing
    m_TrimTdtTable();
    return;
}

void phasedTDT::m_SiteMissRatio()
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

void phasedTDT::m_markMiss()
{
    for(UINT i=0; i<__trioGenos.size(); i++){
	//walk through each site
	for(UINT j=0; j<__trioGenos[i][0].size(); j++){
	    if(__trioGenos[i][0][j]+__trioGenos[i][1][j]+__trioGenos[i][2][j]+__trioGenos[i][3][j]+__trioGenos[i][4][j]+__trioGenos[i][5][j] < 0) //missing
	    {
		for(UINT k=0; k<6; k++) __trioGenos[i][k][j] = 0;
	    }
	}
    }
}

void phasedTDT::m_TrimTdtTable()
{
    for(UINT i=0; i<__TdtCountB.size(); i++){
	for(UINT j=0; j<__TdtCountB[i].size(); j++){
	    if(!__lowMissing[j]) {__TdtCountB[i][j] = 0; __TdtCountC[i][j] = 0;}
	}
    }
    return;
}

// Remeber: the data for GenoPermut may has less variant sites compared to original one, because those sites that have too much missing or uninformative have been trimed;
phasedTDT::phasedTDT(const vector3F& InputData, const vectorF& popMafs, const vectorF& samMafs, std::string& chrom, bool SkipMiss)
{
    __samMafs = samMafs;
    __popMafs = popMafs;
    __chr = chrom;
    __wrongFam = 0;
    __SkipMiss = SkipMiss;
    __trioGenos = InputData;
    __varNum = __trioGenos[0][0].size();
    __lowMissing.assign(__varNum, true);
    //__discardFam.assign(__varNum, 0);
    __flip.assign(__varNum, false);
            
    for(UINT i=0; i<__samMafs.size(); i++){
	if(__samMafs[i] > 0.5) {
	    __flip[i] = true;
	    __samMafs[i] = 1-__samMafs[i];
	}
    }

    // tdt count table for this gene
    m_TdtCountForAllFamilies();
    return;
}

phasedTDT::~phasedTDT(){;}

//Calculate TdtTable for families that passed FamilyCheck
void phasedTDT::m_TdtCountForAllFamilies()
{
    for (UINT i=0; i < __trioGenos.size(); i++)
    {
	//if(Permut) {for(UINT j=0; j< __trioGenos[i].size(); j++) std::cout << __trioGenos[i][j] << std::endl;}
	
        vectorUI ChrOrigin;
        if (m_ChromOrigin(__trioGenos[i], ChrOrigin, false, false)) {;}
        else if (m_ChromOrigin(__trioGenos[i], ChrOrigin, true, false)) {;}
        else if (m_ChromOrigin(__trioGenos[i], ChrOrigin, true, true)) {;}
        else {__wrongFam++; continue;}
        __chromOri.push_back(ChrOrigin);
	
	//set a vector for entire trio (because we split the parent to calculate B&C)
	//in case: 1) both parents missing; 2) skip the missing; 3) entire chromosome missing
	vectorL shouldBeAnalyzed(__trioGenos[i][0].size(), true);
	bool entireMissingFather(true), entireMissingMother(true);
	//in case both parents missing
	for(UINT j=0; j<__trioGenos[i][0].size(); j++){
	    if(__trioGenos[i][0][j]+__trioGenos[i][1][j]+__trioGenos[i][2][j]+__trioGenos[i][3][j] < -20){
		shouldBeAnalyzed[j] = false; // shouldn't analyze these sites
	    }
	    
	    if(__lowMissing[j]) {
		if(__trioGenos[i][0][j] >= 0 && __trioGenos[i][1][j] >= 0) entireMissingFather = false;
		if(__trioGenos[i][2][j] >= 0 && __trioGenos[i][3][j] >= 0) entireMissingMother = false;
	    }
	}
	
	//skip miss
	if(__SkipMiss){
	    for(UINT j=0; j<__trioGenos[i][0].size(); j++){
		if(__trioGenos[i][0][j]+__trioGenos[i][1][j]+__trioGenos[i][2][j]+__trioGenos[i][3][j]+__trioGenos[i][4][j]+__trioGenos[i][5][j] < 0){
		    shouldBeAnalyzed[j] = false; // shouldn't analyze these sites
		}
	    }
	}

        // No genotype error
        //For father, find out which child chromosome comes from father, that is, find out which element in ChrOrigin equals 0 or 1;
        int WhichFromPat = (ChrOrigin[0] <= 1)?0:1;
        vector2F Paternal;
        Paternal.push_back(__trioGenos[i][0]); Paternal.push_back(__trioGenos[i][1]);
        __NonTransChrom.push_back(Paternal[1-ChrOrigin[WhichFromPat]]);

        //For mother, find out which element in ChrOrigin equals 2 or 3;
        int WhichFromMat = 1-WhichFromPat;
        vector2F Maternal;
        Maternal.push_back(__trioGenos[i][2]); Maternal.push_back(__trioGenos[i][3]);
        __NonTransChrom.push_back(Maternal[3-ChrOrigin[WhichFromMat]]);
	
	//std::cout << entireMissingFather << "|" << entireMissingMother << std::endl;
	if(entireMissingFather||entireMissingMother) continue;
	m_TdtCalForOneParent(Paternal, __trioGenos[i][4 + WhichFromPat], ChrOrigin[WhichFromPat], shouldBeAnalyzed);
        m_TdtCalForOneParent(Maternal, __trioGenos[i][4 + WhichFromMat], ChrOrigin[WhichFromMat]-2, shouldBeAnalyzed);
    }
    return;
}


//Calculate transmitted or non-transmitted in one parent and the corresponding child chromosome;
//Chromosome = 0 or 1, indicates the origin of the child allele
void phasedTDT::m_TdtCalForOneParent(vector2F& Parent, vectorF& Child, int Chromosome, vectorL& shouldBeAnalyzed)
{
    vectorF Btable, Ctable;
    Btable.resize(0);
    Ctable.resize(0);

    //std::cout << "P: " << Parent[0] << std::endl << "P: " << Parent[1] << std::endl << "C: " << Child << std::endl;
    for (UINT i = 0; i < Child.size(); i++)
    {
	if(__flip[i]) {Parent[0][i] = m_flip(Parent[0][i]);  Parent[1][i] = m_flip( Parent[1][i]); Child[i] = m_flip(Child[i]);}
        //don't analyze this site
        if(!shouldBeAnalyzed[i]) {Btable.push_back(0); Ctable.push_back(0); continue;}

        double heteroProb = __popMafs[i]*(1-__popMafs[i])*2;
        double homoProb = __popMafs[i]*__popMafs[i];
        double wildProb = (1-__popMafs[i])*(1-__popMafs[i]);

        //Child missing
        if(Child[i] == -9.0) {Btable.push_back(0); Ctable.push_back(0);}
        //Non-informative at this site, missing -9|-9, this is phased data!!
        else if(Parent[Chromosome][i] == -9.0 || Parent[1-Chromosome][i] == -9.0){
            //-9-9 -> 0 ==> 01 -> 0
            if(Child[i] == 0) {Btable.push_back(heteroProb/(wildProb+heteroProb)); Ctable.push_back(0);}
            //-9-9 -> 0 ==> 01 -> 1
            if(Child[i] == 1) {Btable.push_back(0); Ctable.push_back(heteroProb/(homoProb+heteroProb));}
        }
        else{
            //no reverse mutation
            if(Child[i] == 0){
                //01 -> 0
                if(Parent[Chromosome][i] == 0 && Parent[1-Chromosome][i] == 1) {Btable.push_back(1); Ctable.push_back(0);}
                else {Btable.push_back(0); Ctable.push_back(0);}
            }
            //
            else if(Child[i] == 1){
                // 11 -> 1
                if(Parent[Chromosome][i] == 1 && Parent[1-Chromosome][i] == 1) {Btable.push_back(0); Ctable.push_back(0);}
                else if(Parent[Chromosome][i] == 0) {Btable.push_back(0); Ctable.push_back(1); __denovo[i]++;} // de novo mutation
                else if(Parent[Chromosome][i] == 1 && Parent[1-Chromosome][i] == 0) {Btable.push_back(0); Ctable.push_back(1);}
                else {Btable.push_back(0); Ctable.push_back(0);}
            }
            //Informative site: transmitted (Child[i] = 1,2...)
            else {Btable.push_back(0); Ctable.push_back(0);}
        }
    }
    //if(Permut) {std::cout << "b: " << Btable << std::endl << "c: " << Ctable << std::endl;}
    // if the parent is informative
    if(m_IsItInform(Btable) || m_IsItInform(Ctable)) {
	__TdtCountB.push_back(Btable); __TdtCountC.push_back(Ctable);
	//std::cout << "F: " << Parent[0] << std::endl<< "M: " << Parent[1] << std::endl<< "C: " << Child << std::endl;
	//std::cout <<  "B: " << Btable << std::endl <<  "C: " << Ctable << std::endl << "------" << std::endl;
    }
    return;
}

float phasedTDT::m_flip(float before)
{
    float after;
    if(before >= 0 && before <= 1) after = 1-before;
    else after = before;
    return after;
}

bool m_ChromOrigin(const vector2F& Family, vectorUI& kidAlleleOri, bool Missing, bool Denovo)
{
    /*! * Determine the origins of Child's chromosomes
     * Output: 2 ints vector, containing the origin of C1 and C2, eg <1,3>, <0,2>
     * <9,9> means genotype error
     * Implementation:
     */
    vectorF P1(Family[0]);
    vectorF P2(Family[1]);
    vectorF M1(Family[2]);
    vectorF M2(Family[3]);
    //The genotype of affected child;
    vectorF C1(Family[4]);
    vectorF C2(Family[5]);

    kidAlleleOri.resize(0);

    //The chromosome origin of C1 and C2 (0,1,2,or 3, corresponding to Family[0]~~); ConeOrigin = CtwoOrigin = 5 indicates genotype error;
    int ConeOrigin(5), CtwoOrigin(5);
    char Cone, Ctwo;
    std::vector <int> ConeVect, CtwoVect;

    for (int i = 0; i <4; i++){
        if (m_ChromEqual(C1, Family[i], Missing, Denovo)) {ConeVect.push_back(1);}
        else {ConeVect.push_back(0);}

        if (m_ChromEqual(C2, Family[i], Missing, Denovo)) {CtwoVect.push_back(1);}
        else {CtwoVect.push_back(0);}
    }

    //std::cout << "Cone: " << ConeVect << std::endl << "Ctwo: " << CtwoVect << std::endl;

    // Compare C1 to Parental chromosomes, determine where it comes from:
    // N: no parental choromosome is found to be "equal" to the chromosome of interest;
    // P: The chromosome of interest is only found in Father;
    // M: The chromosome of interest is only found in Mother;
    // B: The chromosome of interest is found in both parents;
    if (accumulate(ConeVect.begin(), ConeVect.begin()+4, 0) == 0) {Cone = 'N';}
    else if (accumulate(ConeVect.begin(), ConeVect.begin()+2, 0) > 0){
        if (accumulate(ConeVect.begin()+2, ConeVect.begin()+4, 0) == 0) {Cone = 'P';}
        else {Cone = 'B';}
    }
    else {Cone = 'M';}

    //The same for C2;
    if (accumulate(CtwoVect.begin(), CtwoVect.begin()+4, 0) == 0) {Ctwo = 'N';}
    else if (accumulate(CtwoVect.begin(), CtwoVect.begin()+2, 0) > 0){
        if (accumulate(CtwoVect.begin()+2, CtwoVect.begin()+4, 0) == 0) {Ctwo = 'P';}
        else {Ctwo = 'B';}
    }
    else {Ctwo = 'M';}

    //Combine Cone and Ctwo to determine the origins of child's chromosomes;
    //Genotype Errors;
    if (Cone == 'N' || Ctwo =='N'){
        kidAlleleOri.push_back(5); kidAlleleOri.push_back(5);
        return false;
    }
    //Another situation of genotype error: both child's chromosome come from same parent;
    else if ((Cone == 'P' && Ctwo == 'P') || (Cone == 'M' && Ctwo == 'M')){
        kidAlleleOri.push_back(5); kidAlleleOri.push_back(5);
        return false;
    }

    //No genotype error
    else if (Cone == 'B' && Ctwo == 'B'){
        //child is homozygous, or parents and child are same genotype, and heterozygous;
        if (ConeVect[0] == 1) {ConeOrigin = 0;}
        else {ConeOrigin = 1;}

        if (CtwoVect[2] == 1) {CtwoOrigin = 2;}
        else {CtwoOrigin = 3;}
    }

    //At least one chromosome has been confidentially determined;
    else{
        //Cone comes from Paternal and Ctwo comes from Maternal
        if (Cone == 'P' || Ctwo == 'M'){
            if (ConeVect[0] == 1) {ConeOrigin = 0;}
            else {ConeOrigin = 1;}

            if (CtwoVect[2] == 1) {CtwoOrigin = 2;}
            else {CtwoOrigin = 3;}
        }
        //Cone comes from Maternal and Ctwo comes from Paternal
        else if (Cone == 'M' || Ctwo == 'P'){
            if (ConeVect[2] == 1) {ConeOrigin = 2;}
            else {ConeOrigin = 3;}

            if (CtwoVect[0] == 1) {CtwoOrigin = 0;}
            else {CtwoOrigin = 1;}
        }
        else {;}
    }

    kidAlleleOri.push_back(ConeOrigin); kidAlleleOri.push_back(CtwoOrigin);
    return true;
}

//Define how many errors the "identical" chromosomes can have
bool m_ChromEqual(const vectorF& ChromOne, const vectorF& ChromTwo, bool Missing, bool Denovo)
{
    //We allow at most 1 de novo mutation,
    double AllowedErrorNum = 1;
    int ErrorFound(0);
    bool StrictEqual(true), MissEqual(true), DenovoEqual(true);

    for(UINT i = 0; i < ChromOne.size(); i++){
        if(ChromOne[i] != ChromTwo[i]){
            StrictEqual = false;
            if(ChromOne[i] != -9.0 && ChromTwo[i] != -9.0) { MissEqual = false; ErrorFound++;}
        }
    }
    if(ErrorFound > AllowedErrorNum) DenovoEqual = false;

    // strict equal
    if((!Missing) && (!Denovo)) {return StrictEqual;}
    else if(Missing && (!Denovo)) {return MissEqual;}
    else if(Missing && Denovo) {return DenovoEqual;}
    else {return false;}
}


// get dataset information
UINT phasedTDT::getFamNum()
{
    return __trioGenos.size() - __wrongFam;
}

UINT phasedTDT::getXvariNum()
{
    return __Xvariants;
}

vectorUI phasedTDT::getDenovo()
{
    return __denovo;
}

UINT phasedTDT::getVariNum()
{
    return __varNumToBeAnalyzed;
}

vectorL phasedTDT::getFlip()
{
    return __flip;
}

vectorF phasedTDT::getMissRatio()
{
    return __MissingRatio;
}

vectorL phasedTDT::getBeAnalyzed()
{
    return __lowMissing;
}

vectorF phasedTDT::getTdtCountB(){
    vectorF sum(__varNum, 0);
    for(UINT i=0; i<__TdtCountB.size(); i++){
	for(UINT j=0; j<__TdtCountB[i].size(); j++) sum[j] += __TdtCountB[i][j];
    }
    return sum;
}

vectorF phasedTDT::getTdtCountC(){
    vectorF sum(__varNum, 0);
    for(UINT i=0; i<__TdtCountC.size(); i++){
	for(UINT j=0; j<__TdtCountC[i].size(); j++) sum[j] += __TdtCountC[i][j];
    }
    return sum;
}

vectorF phasedTDT::getSingleP(){
    vectorF Bcount = getTdtCountB();
    vectorF Ccount = getTdtCountC();
    vectorF singleP;

    for(UINT i=0; i<Bcount.size(); i++){
	if((Ccount[i] + Bcount[i]) == 0) singleP.push_back(1);
	else singleP.push_back(gsl_cdf_gaussian_Q((Ccount[i]-Bcount[i])/sqrt(Ccount[i]+Bcount[i]), 1));
    }
    return singleP;
}

vector2F phasedTDT::phasedSST(double mafcutoff)
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
vectorF phasedTDT::TdtMZ(double mafLower, double mafUpper) const
{
    double tdtB(0), tdtC(0), tdtBall(0), tdtCall(0);

    for (UINT i = 0; i < __TdtCountB.size(); i++){
	//std::cout << __TdtCountB[1] << std::endl << __TdtCountC[i] << std::endl;
	for (UINT j = 0; j < __TdtCountB[i].size(); j++){
	    tdtBall += __TdtCountB[i][j]; tdtCall += __TdtCountC[i][j];
	    // ignore if out of [0, upper]
	    //if(__TdtCountB[i][j] > 0) std::cout << "B\t" << i << "\t" << __TdtCountB[i] << std::endl;
	    //if(__TdtCountC[i][j] > 0) std::cout << "C\t" << i << "\t" << __TdtCountC[i] << std::endl;
	    if(__popMafs[j] < mafUpper && __popMafs[j] >= mafLower) {tdtB += __TdtCountB[i][j]; tdtC += __TdtCountC[i][j];}
	}
    }

    vectorF MZresult;

    if(tdtB==0 && tdtC==0) {MZresult.push_back(std::numeric_limits<double>::infinity());}
    else{
	double MZone = gsl_cdf_gaussian_Q((tdtC-tdtB)/sqrt(tdtB+tdtC), 1);
	MZresult.push_back(MZone);
    }

    if(tdtBall==0 && tdtCall==0) {MZresult.push_back(std::numeric_limits<double>::infinity());}
    else{
	double MZoneAll = gsl_cdf_gaussian_Q((tdtCall-tdtBall)/sqrt(tdtBall+tdtCall), 1);
	MZresult.push_back(MZoneAll);
    }
    //std::cout << tdtB << "\t" << tdtC << "\t" << (tdtC-tdtB)/sqrt(tdtB+tdtC) << "\n";
    return MZresult;
}

vectorF phasedTDT::TdtCMC(double mafLower, double mafUpper) const
{
    if(__TdtCountB.size() != __TdtCountC.size()) {std::cerr << "\nError happened when creating 2X2 contingency table. Quit Now!\n"; exit(-1);}

    double tdtB(0), tdtC(0), tdtBall(0), tdtCall(0);
    for(UINT i=0; i<__TdtCountB.size(); i++)
    {
	//std::cout << __TdtCountB[i] << std::endl<< __TdtCountC[i] << std::endl;
	//don't analyze uninformative parent
	//bool informative(false), informativeAll(false);
        double b(0), c(0), bAll(0), cAll(0);
        for(UINT j=0; j<__TdtCountB[i].size(); j++){
	    // all sites
	    //if(__TdtCountB[i][j] == 1.0 || __TdtCountC[i][j] == 1.0) informativeAll=true;
	    bAll += floor(__TdtCountB[i][j]); cAll += floor(__TdtCountC[i][j]);
	    // ignore if out of [0, upper]
	    if(__popMafs[j] < mafUpper && __popMafs[j] >= mafLower) {
		//if(__TdtCountB[i][j] == 1 || __TdtCountC[i][j] == 1) informative=true;
		b += floor(__TdtCountB[i][j]); c += floor(__TdtCountC[i][j]);
	    }
        }
	//std::cout << tdtB << "|" << tdtC << "|" << bAll << "|" << cAll << std::endl;
        if((b+c) != 0) {tdtB += b/(b+c); tdtC += c/(b+c);} //real events exist
	else {;} //inferred event
	
        if((bAll+cAll)!=0) {tdtBall +=  bAll/(bAll+cAll); tdtCall += cAll/(bAll+cAll);}
        else {;}
    }
    
    //std::cout << tdtC << "|" << tdtB << std::endl;
    //std::cout << (tdtC-tdtB)/sqrt(tdtB+tdtC) << "|" << gsl_cdf_gaussian_Q((tdtC-tdtB)/sqrt(tdtB+tdtC), 1) << std::endl;
    vectorF CMCresult;
    if(tdtB==0 && tdtC==0) {CMCresult.push_back(std::numeric_limits<double>::infinity());}
    else{
	double CMCone = gsl_cdf_gaussian_Q((tdtC-tdtB)/sqrt(tdtB+tdtC), 1);
	CMCresult.push_back(CMCone);
    }

    if(tdtBall==0 && tdtCall==0) {CMCresult.push_back(std::numeric_limits<double>::infinity());}
    else{
	double CMConeAll = gsl_cdf_gaussian_Q((tdtCall-tdtBall)/sqrt(tdtBall+tdtCall), 1);
	CMCresult.push_back(CMConeAll);
    }
    return CMCresult;
}


vectorF phasedTDT::TdtMZCMC(double mafLower, double mafUpper) const
{
    vectorF result = TdtMZ(mafLower, mafUpper);
    result.push_back(TdtCMC(mafLower, mafUpper)[0]); result.push_back(TdtCMC(mafLower, mafUpper)[1]);
    return result;
}

bool phasedTDT::m_zeroFamily(double P1, double P2, double M1, double M2, double C1, double C2)
{
    if(__SkipMiss){
	if(P1 == -9.0 || P2 == -9.0 || M1 == -9.0 || M2 == -9.0 || C1 == -9.0 || C2 == -9.0) return true;
	else if(P1 == P2 && P1 == M1 && M1 == M2) return true;
	else return false;
    }
    else{
	if((P1 == -9.0 || P2 == -9.0) && (M1 == -9.0 || M2 == -9.0)) return true;
	else if(M1 == -9.0 || M2 == -9.0) return true;
	else if(P1 == P2 && P1 == M1 && M1 == M2) return true;
	else return false;
    }
}

vectorF phasedTDT::m_tdtWSS(double mafLower, double mafUpper)
{
    //calculate the weight
    vectorI NonZero(__NonTransChrom[0].size(),0), NonMiss(__NonTransChrom[0].size(),0);
    for (UINT i=0; i<__NonTransChrom.size(); i++){
        for (UINT j=0; j < __NonTransChrom[i].size(); j++){
            if(__NonTransChrom[i][j] > 0) NonZero[j] ++;
            if(__NonTransChrom[i][j] >= 0) NonMiss[j] ++;
	    //std::cout << __NonTransChrom[i] << std::endl;
        }
    }
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
    for (UINT i = 0; i < NonZero.size(); i++) __weight.push_back( m_WssWeightCal(int(NotMissAll[i]/2), int(NonMiss[i]/2), NonZero[i]));
    
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
    //std::cout << "-------" << std::endl;

    return res;
}

// VT: return a vector including one-sided chi-square with and without boundary
vectorF phasedTDT::m_tdtVT(double mafLower, double mafUpper, std::string aggregate)
{
    vectorF SortedMAF = __popMafs;

    sort(SortedMAF.begin(), SortedMAF.end());
    SortedMAF.erase(unique(SortedMAF.begin(), SortedMAF.end()), SortedMAF.end());
    __sortedpopMafs = SortedMAF;

    __sortedChiSqu.assign(SortedMAF.size(), 0);

    double VTMaxMZone(-99), VTMaxMZall(-99);
    //if(VTDEBUG) std::cout << "VT: " << SortedMAF << std::endl;

    for(UINT point = 0; point < SortedMAF.size(); point++){
	// keep threshold within defined region
	//if(SortedMAF[point] < mafLower || SortedMAF[point] > mafUpper) continue;
        double TdtB(0), TdtC(0), TdtBall(0), TdtCall(0);
	// count tdtB: only record events on sites that maf less than threshold
	// same as tdtC
        for (UINT i = 0; i < __TdtCountB.size(); i++){
            double b(0), c(0), bAll(0), cAll(0);
            for (UINT j = 0; j < __TdtCountB[i].size(); j++){
                if(__popMafs[j] <= SortedMAF[point]){
		    if(aggregate == "MZ") {bAll += __TdtCountB[i][j]; cAll += __TdtCountC[i][j];}
		    else if(aggregate == "CMC") {bAll += floor(__TdtCountB[i][j]); cAll += floor(__TdtCountC[i][j]);}
		    else {std::cerr << "invalid aggregate argument for VT method" << std::endl; exit(-1);}
		    
		    if(__popMafs[j] <= mafUpper && __popMafs[j] >= mafLower) {
			if(aggregate == "MZ") {b += __TdtCountB[i][j]; c += __TdtCountC[i][j];}
			else if(aggregate == "CMC") {b += floor(__TdtCountB[i][j]); c += floor(__TdtCountC[i][j]);}
			else {std::cerr << "invalid aggregate argument for VT method" << std::endl; exit(-1);}
		    }
		}
            }

            if(aggregate == "MZ") {TdtB += b; TdtC += c; TdtBall += bAll; TdtCall += cAll;}
            else if(aggregate == "CMC") {
                if(b==0&&c==0) continue;
                TdtB += b/(b+c); TdtC += c/(b+c);
		TdtBall += bAll/(bAll+cAll); TdtCall += cAll/(bAll+cAll);
                }
            else {std::cerr << "invalid aggregate argument for VT method" << std::endl; exit(-1);}
        }

	double CQone(0), CQoneAll(0);
	if(TdtB==0 && TdtC==0) {CQone = -99;}
	else {CQone = (TdtC-TdtB)/sqrt(TdtB+TdtC);}
	__sortedChiSqu[point] = CQone;

	if(TdtBall==0 && TdtCall==0) {CQoneAll = -99;}
	else {CQoneAll = (TdtCall-TdtBall)/sqrt(TdtBall+TdtCall);}

	if (CQone > VTMaxMZone) VTMaxMZone = CQone;
	if (CQoneAll > VTMaxMZall) VTMaxMZall = CQoneAll;
	
	//std::cout << SortedMAF[point] << ": " << TdtC << "|" << TdtB << "|" << CQone << "|" << VTMaxMZone << std::endl ;
    }

    //if(VTDEBUG) std::cout << "VTMax: " << VTMaxMZone << " " << VTMaxMZall << std::endl;
    vectorF VTresult;
    VTresult.push_back(VTMaxMZone); VTresult.push_back(VTMaxMZall);
    return VTresult;
}

// trim those uninformative sites, based on missing and whether informative, to speed up permutation
void phasedTDT::TrimUninformSites()
{
    //Before permutation, delete some variant site which are uninformative in data set
    __unInformSites.assign(__varNum, true);
    for(UINT s=0; s<__varNum; s++){
        for(UINT f=0; f<__trioGenos.size(); f++){
	    // keep this site if any trio is informative
            if(!m_zeroFamily(__trioGenos[f][0][s], __trioGenos[f][1][s], __trioGenos[f][2][s], __trioGenos[f][3][s], __trioGenos[f][4][s], __trioGenos[f][5][s])) {__unInformSites[s]=false; break;}
        }
    }
    //std::cout << __unInformSites << std::endl;
    vectorL shouldTrim = __unInformSites;
    //combine with __lowMissing vector
    for(UINT i=0; i<shouldTrim.size(); i++){
	if(!__lowMissing[i]) shouldTrim[i]=true;
    }
    //std::cout << "before trimming: " << std::endl << __popMafs << std::endl << __samMafs << std::endl;
    //std::cout << shouldTrim << std::endl;
    //std::cout << __unInformSites << std::endl;
    //check if there is uninformative site
    bool NeedTrim(false);
    for(UINT i=0; i<__varNum; i++) NeedTrim = (NeedTrim||shouldTrim[i]);
    if(NeedTrim){
        vector3F NewData;
        for(UINT f=0; f<__trioGenos.size(); f++){
            vectorF P1, P2, M1, M2, C1, C2;
            for(UINT s=0; s<__varNum; s++){
                if(!shouldTrim[s])
		{P1.push_back(__trioGenos[f][0][s]); P2.push_back(__trioGenos[f][1][s]); M1.push_back(__trioGenos[f][2][s]);
		 M2.push_back(__trioGenos[f][3][s]); C1.push_back(__trioGenos[f][4][s]); C2.push_back(__trioGenos[f][5][s]);}
            }
            vector2F NewFamily;
            NewFamily.push_back(P1); NewFamily.push_back(P2); NewFamily.push_back(M1); NewFamily.push_back(M2); NewFamily.push_back(C1); NewFamily.push_back(C2);
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
    //std::cout << "after  trimming: " << std::endl << __popMafs << std::endl << __samMafs << std::endl;
    return;
}

// For WSS
vectorF phasedTDT::TdtPermutWSS(double mafLower, double mafUpper, UINT adaptive, UINT PermutateTimes, double alpha, std::string method)
{
    for(UINT i=0; i<__trioGenos.size(); i++){
	// after FamilyCheck, the genotypes martix within a trio should be square matrix, so just check the father's variant site number
	if(__trioGenos[i][0].size() != __varNum) {std::cout << "Error: there are different variant sites number in different trio, can't use phased::TdtWSS method, please check your data" << std::endl; exit(-1);}
    }

    // the VT static for original dataset
    vectorF oriWSS = m_tdtWSS(mafLower, mafUpper);
    if(oriWSS[0] == std::numeric_limits<double>::infinity() || oriWSS[1] == std::numeric_limits<double>::infinity()){
	return oriWSS;
    }
    
    TrimUninformSites();

    UINT permcountWSS(0), permcountWSSall(0);
    double pvalueWSS(9.0), pvalueWSSall(9.0);
    __permTimes = 0;

    for (UINT i=1; i <= PermutateTimes; i++)
    {
        vector3F GenoShuffData;
	if(method == "GenoShuffle") GenoShuffData = m_GenoShuffle();
	else if (method == "HapoShuffle") GenoShuffData = m_HapoShuffle();
	else {std::cout << "ERROR: invalid shuffle method" << std::endl; exit(-1);}
	
	phasedTDT GenoPermut(GenoShuffData, __popMafs, __samMafs, __chr, __SkipMiss);
        vectorF wss = GenoPermut.m_tdtWSS(mafLower, mafUpper);
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

// For VT
vectorF phasedTDT::TdtPermutVT(double mafLower, double mafUpper, UINT adaptive, UINT PermutateTimes, double alpha, std::string VTaggregate, std::string method)
{
    for(UINT i=0; i<__trioGenos.size(); i++){
	// after FamilyCheck, the genotypes martix within a trio should be square matrix, so just check the father's variant site number
	if(__trioGenos[i][0].size() != __varNum) {std::cout << "Error: there are different variant sites number in different trio, can't use phased::TdtVT method, please check your data" << std::endl; exit(-1);}
    }

    // uninformative gene
    vectorF oriVT = m_tdtVT(mafLower, mafUpper, VTaggregate);
    if(oriVT[0] == -99){
	vectorF res(2, std::numeric_limits<double>::infinity());
	return res;
    }
    
    TrimUninformSites();

    UINT permcountVT(0), permcountVTall(0);
    double pvalueVT(9.0), pvalueVTall(9.0);
    __permTimes = 0;

    for (UINT i=1; i <= PermutateTimes; i++)
    {
	//std::cout<< "-------------" << std::endl;
        vector3F GenoShuffData;
	if(method == "GenoShuffle") GenoShuffData = m_GenoShuffle();
	else if (method == "HapoShuffle") GenoShuffData = m_HapoShuffle();
	else {std::cout << "ERROR: invalid shuffle method" << std::endl; exit(-1);}
	
	phasedTDT GenoPermut(GenoShuffData, __popMafs, __samMafs, __chr, __SkipMiss);
	vectorF vt = GenoPermut.m_tdtVT(mafLower, mafUpper, VTaggregate);
	for(UINT t=0; t<vt.size(); t++) {if(vt[t]==-99) vt[t] = 0;}

	//if(VTDEBUG) std::cout << "statistic: " << statistic << std::endl;
	RNG rng;
        gsl_rng* gslr = rng.get();

	//without boundary
	if (vt[0] > oriVT[0]) { permcountVT++;}
        else if (vt[0] == oriVT[0]) {if (gsl_rng_uniform(gslr) > 0.5) permcountVT++;}
	else {;}
	if(pvalueVT > 1.0) pvalueVT = m_check(permcountVT, i, adaptive, alpha);

	if (vt[1] > oriVT[1]) { permcountVTall++;}
        else if (vt[1] == oriVT[1]) {if (gsl_rng_uniform(gslr) > 0.5) permcountVTall++;}
	else {;}
	if(pvalueVTall > 1.0) pvalueVTall = m_check(pvalueVTall, i, adaptive, alpha);

	__permTimes++;
        if (pvalueVT <= 1.0 && pvalueVTall <= 1.0) {break;}
    }

    vectorF VtPermRes;
    VtPermRes.push_back((pvalueVT <= 1.0)? pvalueVT : (1.0 * permcountVT + 1) / (1.0 * __permTimes + 1));
    VtPermRes.push_back((pvalueVTall <= 1.0)? pvalueVTall : (1.0 * pvalueVTall + 1) / (1.0 * __permTimes + 1));

    return VtPermRes;
}

// For MZ
vectorF phasedTDT::TdtPermutMZ(double mafLower, double mafUpper, UINT adaptive, UINT PermutateTimes, double alpha, std::string method)
{
    Permut = true;
    for(UINT i=0; i<__trioGenos.size(); i++){
	// after FamilyCheck, the genotypes martix within a trio should be square matrix, so just check the father's variant site number
	if(__trioGenos[i][0].size() != __varNum) {std::cout << "Error: there are different variant sites number in different trio, can't use phased::TdtMZ method, please check your data" << std::endl; exit(-1);}
    }

    // the VT static for original dataset
    vectorF oriMZ = TdtMZ(mafLower, mafUpper);
    TrimUninformSites();
    
    //std::cout << oriMZ << std::endl;
    UINT permcountMZ(0), permcountMZall(0);
    double pvalueMZ(9.0), pvalueMZall(9.0);
    __permTimes = 0;
    
    for (UINT i=1; i <= PermutateTimes; i++)
    {
        vector3F GenoShuffData;
	if(method == "GenoShuffle") GenoShuffData = m_GenoShuffle();
	else if (method == "HapoShuffle") GenoShuffData = m_HapoShuffle();
	else {std::cout << "ERROR: invalid shuffle method" << std::endl; exit(-1);}
	
	phasedTDT GenoPermut(GenoShuffData, __popMafs, __samMafs, __chr, __SkipMiss);
        vectorF MZ = GenoPermut.TdtMZ(mafLower, mafUpper);
		
	RNG rng;
        gsl_rng* gslr = rng.get();
	
	//with boundary
	if (MZ[0] > oriMZ[0]) { permcountMZ++;}
        else if (MZ[0] == oriMZ[0]) {if (gsl_rng_uniform(gslr) > 0.5) permcountMZ++;}
	else {;}
	if(pvalueMZ > 1.0) pvalueMZ = m_check(permcountMZ, i, adaptive, alpha);
	
	if (MZ[1] > oriMZ[1]) { permcountMZall++;}
        else if (MZ[1] == oriMZ[1]) {if (gsl_rng_uniform(gslr) > 0.5) permcountMZall++;}
	else {;}
	if(pvalueMZall > 1.0) pvalueMZall = m_check(permcountMZall, i, adaptive, alpha);
	
	//std::cout << MZ << std::endl;
	__permTimes++;
        if (pvalueMZ <= 1.0 && pvalueMZall <= 1.0) {break;}
    }
  
    vectorF MZPermRes;
    MZPermRes.push_back((pvalueMZ <= 1.0)? pvalueMZ : (1.0 * permcountMZ + 1) / (1.0 * __permTimes + 1));
    MZPermRes.push_back((pvalueMZall <= 1.0)? pvalueMZall : (1.0 * permcountMZall + 1) / (1.0 * __permTimes + 1));
    
    return MZPermRes;
}

// For WSS VT-MZ VT-CMC
vectorF phasedTDT::TdtPermut(double mafLower, double mafUpper, UINT adaptive, UINT PermutateTimes, double alpha, std::string method)
{
    for(UINT i=0; i<__trioGenos.size(); i++){
	// after FamilyCheck, the genotypes martix within a trio should be square matrix, so just check the father's variant site number
	if(__trioGenos[i][0].size() != __varNum) {std::cout << "Error: there are different variant sites number in different trio, can't use phased::TdtPermut method, please check your data" << std::endl; exit(-1);}
    }

    // the VT static for original dataset
    double oriWSS = m_tdtWSS(mafLower, mafUpper)[0];
    double oriVTMZ = m_tdtVT(mafLower, mafUpper, "MZ")[0];
    double oriVTCMC = m_tdtVT(mafLower, mafUpper, "CMC")[0];
    std::cout << oriWSS << "\t" << oriVTMZ << "\t" << oriVTCMC <<  std::endl;
    if(oriWSS == std::numeric_limits<double>::infinity() || oriVTMZ == -99 || oriVTCMC == -99)
    {
	vectorF res(3, std::numeric_limits<double>::infinity());
	return res;
    }
    TrimUninformSites();

    UINT permcountWSS(0), permcountVTMZ(0), permcountVTCMC(0);
    double pvalueWSS(9.0), pvalueVTMZ(9.0), pvalueVTCMC(9.0);
    __permTimes = 0;

    for (UINT i=1; i <= PermutateTimes; i++)
    {
        vector3F GenoShuffData;
	if(method == "GenoShuffle") GenoShuffData = m_GenoShuffle();
	else if (method == "HapoShuffle") GenoShuffData = m_HapoShuffle();
	else {std::cout << "ERROR: invalid shuffle method" << std::endl; exit(-1);}
	
	phasedTDT GenoPermut(GenoShuffData, __popMafs, __samMafs, __chr, __SkipMiss);

	double wss = GenoPermut.m_tdtWSS(mafLower, mafUpper)[0];
	//std::cout << wss << "\t";
	if(wss==std::numeric_limits<double>::infinity()) wss = 0;
	
	double vtmz = GenoPermut.m_tdtVT(mafLower, mafUpper, "MZ")[0];
	//std::cout << vtmz << "\t";
	if(vtmz==-99) vtmz = 0;
	
	double vtcmc = GenoPermut.m_tdtVT(mafLower, mafUpper, "CMC")[0];
	//std::cout << vtcmc << "\n";
	if(vtcmc==-99) vtcmc = 0;

	//std::cout << wss << "\t" << vtmz << "\t" << vtcmc << std::endl;
	
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

	//VT-CMC
	if (vtcmc> oriVTCMC) { permcountVTCMC++;}
        else if (vtcmc == oriVTCMC) {if (gsl_rng_uniform(gslr) > 0.5) permcountVTCMC++;}
	else {;}
	if(pvalueVTCMC > 1.0) pvalueVTCMC = m_check(permcountVTCMC, i, adaptive, alpha);

	__permTimes++;
        if (pvalueWSS <= 1.0 && pvalueVTMZ <= 1.0 && pvalueVTCMC<= 1.0) {break;}
    }
    //std::cout << std::endl;
    vectorF PermRes;
    PermRes.push_back((pvalueWSS <= 1.0)? pvalueWSS : (1.0 * permcountWSS + 1) / (1.0 * __permTimes + 1));
    PermRes.push_back((pvalueVTMZ <= 1.0)? pvalueVTMZ : (1.0 * permcountVTMZ + 1) / (1.0 * __permTimes + 1));
    PermRes.push_back((pvalueVTCMC <= 1.0)? pvalueVTCMC : (1.0 * permcountVTCMC + 1) / (1.0 * __permTimes + 1));

    return PermRes;
}

// T is the number of SNPs or genes
double phasedTDT::m_check(UINT permcount, UINT iPermut, UINT adaptive, double alpha)
{
    if(iPermut % adaptive != 0 || iPermut == 0) {
        return 9.0;
    }

    //double alpha = 0.05;
    //double pval = permcount*1.0/iPermut*1.0;
    //double z = gsl_cdf_gaussian_Pinv(1.0-alpha/2.0, 1.0);
    //double zsq = z*z;
    //double plw = (pval + zsq / (2.0*iPermut) - z * sqrt((pval*(1.0-pval)+zsq/(4.0*iPermut))/(1.0*iPermut))) / (1.0+zsq/(1.0*iPermut));

    double pval = permcount*1.0/iPermut*1.0;
    double sigma = sqrt(pval*(1-pval)/iPermut);

    double beta = 0.05;
    double gs = gsl_cdf_gaussian_Pinv(1.0-beta/2.0, sigma);
    //return (pval > alpha + 6*sigma)? pval:9.0;
    return (pval - gs > alpha)? pval:9.0;
}

vectorF phasedTDT::getWSSweight()
{
    return __weight;
}

vectorF phasedTDT::getSortedMafs()
{
    return __sortedpopMafs;
}

vectorF phasedTDT::getSortedCQ()
{
    return __sortedChiSqu;
}

UINT phasedTDT::getPermTimes()
{
    return __permTimes;
}

//Shuffle the genotupe, that is shuffle single base and each site
vector3F phasedTDT::m_GenoShuffle() const
{
    if(__trioGenos.size() != __chromOri.size()) {std::cerr << "unmatched family size and __chromOri vector size, unable to genotype shuffle. Quit Now!" << std::endl; exit(-1);}
    vector3F NewData;
    UINT VariantNum = __trioGenos[0][0].size();

    for (UINT OneFam = 0; OneFam < __trioGenos.size(); OneFam++)
    {
        //std::cout << "Family " << OneFam+1 << "; before shuffle: " << std::endl;
        //for(UINT i=0; i<__trioGenos[OneFam].size(); i++) std::cout << __trioGenos[OneFam][i] << std::endl;
	//std::cout << "----" << std::endl;
        //genotype error in this family, cannot determine chromosome origin
        if(__chromOri[OneFam][0] == 5 || __chromOri[OneFam][1] == 5) {continue;}

        vector2F NewParents;
        NewParents.resize(4);
        for (UINT i = 0; i < 4; i++) NewParents[i].resize(VariantNum);

        //whether missing happens, and whether uninformative because missing
        bool Missing(false), UninformFlag(false);
        vectorL FatMissing(VariantNum, false), MotMissing(VariantNum, false), Uninform(VariantNum, false);

        for (UINT EachSite = 0; EachSite < VariantNum; EachSite++)
        {
	    //kid missing or parents both missing
	    if((__trioGenos[OneFam][4][EachSite] == -9 || __trioGenos[OneFam][5][EachSite] == -9) || ((__trioGenos[OneFam][0][EachSite] == -9 || __trioGenos[OneFam][1][EachSite] == -9) && (__trioGenos[OneFam][2][EachSite] == -9 || __trioGenos[OneFam][3][EachSite] == -9)))
	    {
		Uninform[EachSite] = true;
		for(UINT i = 0; i < 4; i++) NewParents[i][EachSite] = -9;
		UninformFlag = true;
	    }
            else{
                //Don't need to permut if father = mother!, even father = mother = -9|-9
                double one = __trioGenos[OneFam][0][EachSite];
                bool allthesame(true);
                for (UINT i = 1; i < 4; i++){
                    if(__trioGenos[OneFam][i][EachSite] != one) allthesame = false;
                }

                if(allthesame) {for(UINT i = 0; i < 4; i++) NewParents[i][EachSite] = one;}
                // no missing
                else if(__trioGenos[OneFam][0][EachSite] != -9 && __trioGenos[OneFam][1][EachSite] != -9 && __trioGenos[OneFam][2][EachSite] != -9 && __trioGenos[OneFam][3][EachSite] != -9)
                {
                    vectorF SiteGeno;
                    for (UINT i = 0; i < 4; i++) SiteGeno.push_back(__trioGenos[OneFam][i][EachSite]);
                    random_shuffle (SiteGeno.begin(), SiteGeno.begin()+4);
                    for (UINT i = 0; i < 4; i++) NewParents[i][EachSite] = SiteGeno[i];
                }
                // missing father
                else if(__trioGenos[OneFam][0][EachSite] == -9 || __trioGenos[OneFam][1][EachSite] == -9)
                {
                    vectorF SiteGeno;
                    UINT whichFromFather = __chromOri[OneFam][0] <= 1?0:1;
		    // doesn't matter, we will go back this site again if there is missing
                    SiteGeno.push_back(__trioGenos[OneFam][4+whichFromFather][EachSite]); SiteGeno.push_back(__trioGenos[OneFam][4+whichFromFather][EachSite]);
                    //mother's genotype
                    SiteGeno.push_back(__trioGenos[OneFam][2][EachSite]); SiteGeno.push_back(__trioGenos[OneFam][3][EachSite]);
                    //just shuffle 3 hapotypes
                    random_shuffle (SiteGeno.begin()+1, SiteGeno.begin()+4);
                    for (UINT i = 0; i < 4; i++) NewParents[i][EachSite] = SiteGeno[i];
                    //record missing events
                    FatMissing[EachSite] = true;
                    Missing = true;
                }
                // missing mother
                else if(__trioGenos[OneFam][2][EachSite] == -9 || __trioGenos[OneFam][3][EachSite] == -9)
                {
                    vectorF SiteGeno;
                    //father's genotype
                    SiteGeno.push_back(__trioGenos[OneFam][0][EachSite]); SiteGeno.push_back(__trioGenos[OneFam][1][EachSite]);
                    UINT whichFromMother = __chromOri[OneFam][0] >= 2?0:1;
                    SiteGeno.push_back(__trioGenos[OneFam][4+whichFromMother][EachSite]); SiteGeno.push_back(__trioGenos[OneFam][4+whichFromMother][EachSite]);
                    //just shuffle 3 hapotypes
                    random_shuffle (SiteGeno.begin(), SiteGeno.begin()+3);
                    for (UINT i = 0; i < 4; i++) NewParents[i][EachSite] = SiteGeno[i];
                    //record missing events
                    MotMissing[EachSite] = true;
                    Missing = true;
                }
                else {for(UINT i = 0; i < 4; i++) NewParents[i][EachSite] = 0;}
            }
        }
        //The chromosome comes from father should be NewFamily[0 or 1]
        int ChildOneChr = rand()%2;
        int ChildTwoChr = rand()%2 +2;
        NewParents.push_back(NewParents[ChildOneChr]);
        NewParents.push_back(NewParents[ChildTwoChr]);

        //remark missing
        if(Missing){
            for(UINT i=0; i<FatMissing.size(); i++){
                if(FatMissing[i]) {NewParents[0][i] = -9.0; NewParents[1][i] = -9.0;}
                if(MotMissing[i]) {NewParents[2][i] = -9.0; NewParents[3][i] = -9.0;}
            }
        }
        else {;}
	// reconstruct uninformative sites: for un-transmitted maf
	//std::cout << Uninform << std::endl;
	if(UninformFlag){
            for(UINT i=0; i<Uninform.size(); i++){
                if(Uninform[i]) {
		    for(UINT j = 0; j < 6; j++) NewParents[j][i] = __trioGenos[OneFam][j][i];
		}
            }
	}
        //std::cout << "Family " << OneFam+1 << "; after shuffle: " << std::endl;
        //for(UINT i=0; i<NewParents.size(); i++) std::cout << NewParents[i] << std::endl;
	//std::cout << "=====" << std::endl;
        NewData.push_back(NewParents);
    }
    return NewData;
}

// shuffle hapotype
//Shuffle the hapotype, that is shuffle the chromosome
vector3F phasedTDT::m_HapoShuffle() const
{
    if(__trioGenos.size() != __chromOri.size()) {std::cerr << "unmatched family size and __chromOri vector size, unable to hapotype shuffle. Quit Now!" << std::endl; exit(-1);}
    vector3F NewData;
    UINT VariantNum = __trioGenos[0][0].size();
    
    RNG rng;
    gsl_rng* gslr = rng.get();
	
    for (UINT OneFam = 0; OneFam < __trioGenos.size(); OneFam++)
    {
        //for(UINT i=0; i<__trioGenos[OneFam].size(); i++) std::cout << __trioGenos[OneFam][i] << std::endl;
	//std::cout << "----" << std::endl;
	
	if(__chromOri[OneFam][0] == 5 || __chromOri[OneFam][1] == 5) {continue;}
	
	//frist, record where missing happens, if the site is informtive
	//second, infer missing genotype in parents genotype
	//then. re-construct missing
        bool Missing(false), UninformFlag(false);
        vectorL FatMissing(VariantNum, false), MotMissing(VariantNum, false), Uninform(VariantNum, false);
	
	vector2F NewParent; NewParent.resize(4);
        copy( __trioGenos[OneFam].begin(), __trioGenos[OneFam].begin()+4, NewParent.begin() );
	
	for(UINT EachSite = 0; EachSite < VariantNum; EachSite++){
	    if(__trioGenos[OneFam][0][EachSite] == -9 || __trioGenos[OneFam][1][EachSite] == -9) {FatMissing[EachSite] = true; Missing = true;} // father missing
	    if(__trioGenos[OneFam][2][EachSite] == -9 || __trioGenos[OneFam][3][EachSite] == -9) {MotMissing[EachSite] = true; Missing = true;} // mothing missing
	    if((__trioGenos[OneFam][4][EachSite] == -9 || __trioGenos[OneFam][5][EachSite] == -9) || (FatMissing[EachSite] && MotMissing[EachSite])) {Uninform[EachSite] = true; UninformFlag = true;} // kid missing or both P missing
	    
	    //infer missing, if it's still informative
	    if(!Uninform[EachSite]){
		if(FatMissing[EachSite]){
		    UINT whichFromFather = __chromOri[OneFam][0] <= 1?0:1;
		    NewParent[__chromOri[OneFam][whichFromFather]][EachSite] = __trioGenos[OneFam][4+whichFromFather][EachSite];
		    if(gsl_rng_uniform(gslr) < __popMafs[EachSite]) NewParent[1-__chromOri[OneFam][whichFromFather]][EachSite] = 1;
		    else NewParent[1-__chromOri[OneFam][whichFromFather]][EachSite] = 0;
		}
		if(MotMissing[EachSite]){
		    UINT whichFromMother = __chromOri[OneFam][0] >= 2?0:1;
		    NewParent[__chromOri[OneFam][whichFromMother]][EachSite] = __trioGenos[OneFam][4+whichFromMother][EachSite];
		    if(gsl_rng_uniform(gslr) < __popMafs[EachSite]) NewParent[5-__chromOri[OneFam][whichFromMother]][EachSite] = 1;
		    else NewParent[5-__chromOri[OneFam][whichFromMother]][EachSite] = 0;
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

        //re-construct missing
        if(Missing){
            for(UINT i=0; i<FatMissing.size(); i++){
                if(FatMissing[i]) {NewParent[0][i] = -9.0; NewParent[1][i] = -9.0;}
                if(MotMissing[i]) {NewParent[2][i] = -9.0; NewParent[3][i] = -9.0;}
            }
        }
        else {;}
	// reconstruct uninformative sites: for un-transmitted maf
	if(UninformFlag){
            for(UINT i=0; i<Uninform.size(); i++){
                if(Uninform[i]) {
		    for(UINT j = 0; j < 6; j++) NewParent[j][i] = __trioGenos[OneFam][j][i];
		}
            }
	}
        //for(UINT i=0; i<NewParent.size(); i++) std::cout << NewParent[i] << std::endl;
	//std::cout << "=====" << std::endl;
	
        NewData.push_back(NewParent);
    }

    return NewData;
}

vector3F phasedTDT::m_geTrio3D(const vector2F& genos, UINT famSize)
{
    vector3F trios;
    for(UINT i = 0; i < genos.size()/famSize; i++)
    {
	// get one trio raw genotypes
        vector2F OneFamily;
        OneFamily.resize(famSize);
        copy(genos.begin() + i*famSize, genos.begin() + (i+1)*famSize, OneFamily.begin());

        if(m_FamilyCheck(OneFamily)){
            // cut off this site if one of trio has missing in this site
            trios.push_back(OneFamily);
        }
        else{
	    __wrongFam++;
	}
    }
    return trios;
}

// return value: false means not square; true means normal
bool phasedTDT::m_FamilyCheck(const vector2F& OneFam) const
{
    UINT VariantNum = OneFam[0].size();

    for(UINT i = 0; i < OneFam.size(); i++)
    {
	// not a square
        if(OneFam[i].size() != VariantNum) {return false;}
    }
    return true;
}
