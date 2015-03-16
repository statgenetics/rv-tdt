
#include <sstream>
template <typename T>
std::string n2s ( T Number )
{
  std::stringstream ss;
  ss << Number;
  return ss.str();
}

void zx_createMissing(vector2F& genotypes, double indMissRatio, double siteMissRatio, gsl_rng* gslr);
vector2F zx_genGenos(const vector2F& affectedOffsprings, const vectorF& keepedMafs, gsl_rng* gslr);
vector2F zx_splitOffspringGeno(const vectorF& offspringHap, gsl_rng* gslr);
vectorF zx_generateHapo(const vectorF& mafs, const vectorL& Missing, gsl_rng* gslr);
vector2F zx_phased2unphased(const vector2F& phased);

vector2F zx_phased2unphased(const vector2F& phased)
{
    vector2F unphased;
    for(UINT i=0; i<phased.size()/2; i++){
        vectorF oneInd;
        for(UINT j=0; j<phased[2*i].size(); j++){
            oneInd.push_back(phased[2*i][j] + phased[2*i+1][j]);
        }
        unphased.push_back(oneInd);
    }
    return unphased;
}


vector2F zx_genGenos(const vector2F& affectedOffsprings, const vectorF& keepedMafs, gsl_rng* gslr)
{
    /*! * Generate 2d genotypes from simulated case data
     * If the simulated genotype marks missing site with -9: keep to be -9 in all member in that trio
     * Implementation:
     */
    vector2F genotypes;
    
    for(UINT ind=0; ind<affectedOffsprings.size(); ind++){
        vectorL MissingList(affectedOffsprings[ind].size(), false);
        for(UINT j=0; j<affectedOffsprings[ind].size(); j++){
            if(affectedOffsprings[ind][j] == -9.0) {MissingList[j] = true;}
        }
        // simulated data code as 0/1/2, split it into two vector with 0/1 code
        vector2F offsprings = zx_splitOffspringGeno(affectedOffsprings[ind], gslr);
        
        vectorF randomOne = zx_generateHapo(keepedMafs, MissingList, gslr);
        vectorF randomTwo = zx_generateHapo(keepedMafs, MissingList, gslr);
 
        double whoGoesToFather = gsl_rng_uniform(gslr);
        //offsprings[0] goes into father;
        if(whoGoesToFather < 0.5){
            double whichFatAllele = gsl_rng_uniform(gslr);
            //offspring[0] goes into father[0];
            if(whichFatAllele < 0.5){
                genotypes.push_back(offsprings[0]);
                genotypes.push_back(randomOne);
            }
            else{
                genotypes.push_back(randomOne);
                genotypes.push_back(offsprings[0]);
            }
    
            double whichMotAllele = gsl_rng_uniform(gslr);
            //offspring[1] goes into mother[0]
            if(whichMotAllele < 0.5){
                genotypes.push_back(offsprings[1]);
                genotypes.push_back(randomTwo);
            }
            else{
                genotypes.push_back(randomTwo);
                genotypes.push_back(offsprings[1]);
            }
        }
        //offspring[0] goes into mother;
        else{
            double whichFatAllele = gsl_rng_uniform(gslr);
            //offspring[1] goes into father[0];
            if(whichFatAllele < 0.5){
                genotypes.push_back(offsprings[1]);
                genotypes.push_back(randomOne);
            }
            else{
                genotypes.push_back(randomOne);
                genotypes.push_back(offsprings[1]);
            }
    
            double whichMotAllele = gsl_rng_uniform(gslr);
            if(whichMotAllele < 0.5){
                genotypes.push_back(offsprings[0]);
                genotypes.push_back(randomTwo);
            }
            else{
                genotypes.push_back(randomTwo);
                genotypes.push_back(offsprings[0]);
            }
        }
        // the kid
        genotypes.push_back(offsprings[0]); genotypes.push_back(offsprings[1]);
    }
    return genotypes;
}

vector2F zx_splitOffspringGeno(const vectorF& offspringHap, gsl_rng* gslr)
{
    vectorF hapOne, hapTwo;
    for(UINT i = 0; i < offspringHap.size(); i++){
        if(offspringHap[i] == 0.0){
            hapOne.push_back(0.0);
            hapTwo.push_back(0.0);
        }
        else if(offspringHap[i] == 1.0){
            double runif = gsl_rng_uniform(gslr);
            if(runif < 0.5){
                hapOne.push_back(0.0);
                hapTwo.push_back(1.0);
            }
            else{
                hapOne.push_back(1.0);
                hapTwo.push_back(0.0);
            }
        }
        else if(offspringHap[i] == 2.0){
            hapOne.push_back(1.0);
            hapTwo.push_back(1.0);
        }
        else if(offspringHap[i] == -9.0){
            hapOne.push_back(-9.0);
            hapTwo.push_back(-9.0);
        }
        else { std::cerr << "Only (0,1,2 or -9) genotype code is allowed, please check the input file" << std::endl; exit(-1);}
    }

    vector2F splited;
    splited.push_back(hapOne); splited.push_back(hapTwo);
    return splited;
}

vectorF zx_generateHapo(const vectorF& mafs, const vectorL& Missing, gsl_rng* gslr)
{
    vectorF hapo;
    for(UINT i = 0; i < Missing.size(); i++){
        if(Missing[i] == true) {hapo.push_back(-9.0);}
        else{
            double randomNum = gsl_rng_uniform(gslr);
            if (randomNum <= mafs[i]) {hapo.push_back(1.0);}
            else {hapo.push_back(0.0);}
        }
    }
    return hapo;
}

void zx_createMissing(vector2F& genotypes, double indMissRatio, double siteMissRatio, gsl_rng* gslr)
{
    if(indMissRatio > 1 || indMissRatio < 0 || siteMissRatio > 1 || siteMissRatio < 0) {std::cerr << "\nInvalid missing arguments. Quit Now!\n"; exit(-1);}
    
    if(indMissRatio > 0){
        //total family number genotypes.size()/6
        for(UINT p=0; p<genotypes.size()/6; p++){
            double prob = gsl_rng_uniform(gslr);
            if(prob <= indMissRatio){
                double runif = gsl_rng_uniform(gslr);
                //very tricky: have one missing parent per trio
                if(runif <= 0.5){
                    for(UINT j=0; j < genotypes[p*6].size(); j++) {genotypes[p*6][j] = -9; genotypes[p*6+1][j] = -9;}
                } //missing father
                else {for(UINT j=0; j < genotypes[p*6+2].size(); j++) {genotypes[p*6+2][j] = -9; genotypes[p*6+3][j] = -9;}}
            }
        }
    }
    
    if(siteMissRatio > 0){
        for(UINT s=0; s<genotypes[0].size(); s++){
            /*double ForThisSite = gsl_rng_uniform(gslr)*siteMissRatio*2; //mean is siteMissRatio, uniformly distributed in [0, siteMissRatio*2]*/
            //total parent number
            for(UINT p=0; p<genotypes.size()/3; p++){
                double prob = gsl_rng_uniform(gslr);
                //only parent missing
                if(prob <= siteMissRatio) {genotypes[(p/2)*6+(p%2)*2][s] = -9; genotypes[(p/2)*6+(p%2)*2+1][s] = -9;}
            }
        }
    }
    return;
}

// variants statics summary: by order: dataIdx originalNum  AnalyzedSite VariantNotInTrio, VariantInParent, VariantInKid, CommonVariantInTrio, #ofTrioHasVariantInParent, #ofTrioHasVariantInKid, MissInParent, MissInSite
vectorF zx_vStatSum(UINT dataIdx, double boundary, vectorF& mafs, vectorF& keepedMafs, vector2F& genotypesTDT)
{
    UINT siteNum = genotypesTDT[0].size();
    vectorL VariantInParents(siteNum, false), VariantInKid(siteNum, false);
    vectorL trioHasVarInParent(genotypesTDT.size()/6, false), trioHasVarInKid(genotypesTDT.size()/6, false);
    vectorF MissInSite(siteNum, 0);
    double MissInParent(0);

    for(UINT i=0; i<genotypesTDT.size()/6; i++){
        for(UINT j=0; j<genotypesTDT[6*i].size(); j++)
        {
            if(genotypesTDT[6*i][j]>0 || genotypesTDT[6*i+1][j]>0 || genotypesTDT[6*i+2][j]>0 || genotypesTDT[6*i+3][j]>0) {VariantInParents[j] = true; trioHasVarInParent[i] = true;} //variant found in parents
            if(genotypesTDT[6*i+4][j]>0 || genotypesTDT[6*i+5][j]>0) {VariantInKid[j] = true; trioHasVarInKid[i] = true;} //variant found in kid
            for(UINT p=0; p<4; p++) {if(genotypesTDT[6*i+p][j] < 0) MissInSite[j] += 1;} //record missing in variant site
        }
        for(UINT p=0; p<4; p++) {if(accumulate(genotypesTDT[6*i+p].begin(), genotypesTDT[6*i+p].end(), 0) == int((-9)*siteNum)) MissInParent += 1;} //record missing in entire chromosome
    }

    // rare variant static
    vectorL CommonVariantInTrio;
    vectorF staticSum(11,0); //
    for(UINT i=0; i<VariantInParents.size(); i++){
        if((!VariantInParents[i]) && (!VariantInKid[i])) staticSum[3] += 1;
        if(VariantInParents[i]) staticSum[4] += 1;
        if(VariantInKid[i]) staticSum[5] += 1;
        if((VariantInParents[i]||VariantInKid[i]) && keepedMafs[i] >= boundary) staticSum[6] += 1;
    }
    for(UINT i=0; i<genotypesTDT.size()/6; i++){
        if(trioHasVarInParent[i]) staticSum[7]++;
        if(trioHasVarInKid[i]) staticSum[8]++;
    }
    staticSum[0] = dataIdx+1;
    staticSum[1] = mafs.size();
    staticSum[2] = keepedMafs.size();
    staticSum[9] = MissInParent/((genotypesTDT.size()/6)*4); //family number: genotypesTDT.size()/6
    staticSum[10] = accumulate(MissInSite.begin(), MissInSite.end(), 0)/float((genotypesTDT.size()/6)*4*siteNum); //total parental haplotypes: (genotypesTDT.size()/6)*4
    return staticSum;
}

UINT ldGenHaplo(vectorF& prob , gsl_rng* gslr);
vectorF breakLD(vectorF& prev, double prob, gsl_rng* gslr);

UINT ldGenHaplo(vectorF& prob , gsl_rng* gslr)
{
    double runif = gsl_rng_uniform(gslr);
    double low(0), high(prob[0]);
    UINT index(0);
    for(UINT i=0; i<prob.size(); i++){
        if(low <= runif && runif <= high) { index = i; break; }
        else {
            low=high;
            high+=prob[i+1];
        }
    }
    /*std::cout << runif << "--> " << index << std::endl;*/
    return index;
}


vectorF breakLD(vectorF& prev, double prob, gsl_rng* gslr)
{
    vectorF after;
    for(UINT i=0; i<prev.size(); i++){
        if(prev[i] > 0) {
            double runif = gsl_rng_uniform(gslr);
            if(runif < prob) after.push_back(1-prev[i]);
            else after.push_back(prev[i]);
        } else after.push_back(prev[i]);
    }
    return after;
}

