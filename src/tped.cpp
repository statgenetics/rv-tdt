#include "tped.h"

pedFileReader::~pedFileReader() {;}

pedFileReader::pedFileReader(std::string genoFile, bool phased)
{
	__phased = phased;
	if(!m_scanTped(genoFile, __rawGenos)){ 
		std::cerr << "\tUnable to source [ " << genoFile << " ]. Quit.\n"; exit(-1);
	}

	// by default, no mapping file; consider as single gene to set up __genes;
	mapSF noMafInfo;
	BOOST_FOREACH(pairSvF pa, __rawGenos){
		noMafInfo.insert(std::make_pair(pa.first, std::numeric_limits<double>::quiet_NaN()));
	}
	ngsGene singleGene("YourGene", __rawGenos, noMafInfo); //TODO: should use regular expression to figure out gene name from input genoFile name
	__genes.push_back(singleGene);
	return;
}

bool pedFileReader::m_scanTped(std::string filename, mapSvF& geno)
{
	vectorSS data;
	if(!loadFile(filename, data)){std::cerr << "\tUnable to source [ " << filename << " ]. Quit." << std::endl; return false;}
	__indNum = (data[0].size() - 1)/(__phased?2:1);
	// data check: same length
	UINT size = data[0].size();  
	for(UINT i=0; i<data.size(); i++){
		if(data[i].size() != size) {
			std::cerr << "Invalid genotype data, must be same length for each variant" << std::endl; return false;
		}

		std::string key = data[i][0];
		vectorF genotypes;
		for(UINT j=1; j<data[i].size(); j++) genotypes.push_back(atof(data[i][j].c_str()));
		if(geno.find(key) == geno.end()){
			geno[key] = genotypes;
		} else {__errorMessage += "duplicate genotype in variant " + key + "\n";}
	}
	return true;
}

void pedFileReader::loadMapFile(std::string mapFile)
{
  vectorSS data;
  if(!loadFile(mapFile, data)){std::cerr << "\tUnable to source [ " << mapFile << " ]. Quit." << std::endl; exit(-1);}
  
  // first setup geneTosnp map and snpTomaf map
  mapSF mafMap;
  mapSvS geneMap;
  for(UINT i=0; i<data.size(); i++){
	  if(data[i].size() < 2) continue; // if column number is less than 2, skip
	  std::string gene = data[i][0];
	  std::string snp = data[i][1];
	  double maf;
	  // if maf == NA, replace with singleton maf 
	  if(data[i].size() == 2) maf = 1.0/(__indNum*2.0);
	  else if(atof(data[i][2].c_str()) != 0) maf = atof(data[i][2].c_str());
	  else maf = 1.0/(__indNum*2.0); // invalid number, including 'NA'
	  if(mafMap.find(snp) == mafMap.end()) mafMap[snp] = maf;
	  if(geneMap.find(gene) == geneMap.end()){
		  vectorS snps; snps.push_back(snp);
		  geneMap[gene] = snps;
	  } else {
		  if(std::find(geneMap[gene].begin(), geneMap[gene].end(), snp) == geneMap[gene].end()) geneMap[gene].push_back(snp);
		  else {__errorMessage += "duplicate variant " + snp + "in gene " + gene + " in mapping file\n";}
	  }
  }
  
  // then split __rawGenos to setup __genes; consider as multiple genes
  __genes.clear(); // clear content: "YourGene" by default
  std::map<std::string, std::vector< std::string > >::const_iterator iter = geneMap.begin();
  for(; iter != geneMap.end(); ++iter){
	  std::string geneName = iter->first;
	  vectorS snpList = iter->second;
	  //std::cout << snpList << std::endl;
	  mapSvF genos;
	  mapSF mafs;
	  BOOST_FOREACH(std::string snp, snpList){
		  if(__rawGenos.find(snp) != __rawGenos.end()) genos.insert(std::make_pair(snp, __rawGenos[snp]));
		  else {__errorMessage += "SNP " + snp + " not found for " + geneName + " in genotype file\n";}

		  if(mafMap.find(snp) != mafMap.end()) mafs.insert(std::make_pair(snp, mafMap[snp]));
		  else {__errorMessage += "SNP " + snp + " not found (to have maf) for " + geneName + "i n genotype file\n";}
	  }

	  if(genos.size() != mafs.size()) {__errorMessage += "gene " + geneName + " has incompatible genotype data and maf information\n"; continue;}
	  else {
		  ngsGene oneGene(geneName, genos, mafs); 
		  __genes.push_back(oneGene);
	  }
  } //ngsGene(std::string _name, mapSvF _genoP, mapSF _maf){
  return;
}

// update: keep phenotype list to be string; 5/17/2012
void pedFileReader::loadPhenoFile(std::string phenoFile)
{
	vectorSS data;
	if(!loadFile(phenoFile, data)){std::cerr << "\tUnable to source [ " << phenoFile << " ]. Quit." << std::endl; exit(-1);}

	for(UINT i=0; i<data.size(); i++){
		std::string ind = data[i][0];
		vectorS phenos;
		for(UINT j=1; j<data[i].size(); j++) phenos.push_back(data[i][j]);
		if(__phenotypes.find(ind) == __phenotypes.end()) {
			__phenotypes[ind] = phenos; 
			//std::cout << __phenotypes[ind] << std::endl;
			__indOrder.push_back(ind);
		}
		else __errorMessage += "duplicate individual " + ind + " in phenotype file\n";
	}

	// check if the phenotype valid
	std::map<std::string, std::vector< double > >::const_iterator iter = __rawGenos.begin();
	int ph = __phased?2:1;
	if(__phenotypes.size()*ph != (iter->second).size()) {
		__errorMessage += "unable to load phenotype data" + phenoFile + ", unmatched size with genotype\n";
		__phenotypes.clear();
	}
	return;
}


void pedFileReader::loadMapPheno(std::string mapFile, std::string phenoFile)
{
	loadMapFile(mapFile);
	loadPhenoFile(phenoFile);
	return;
}

inline vector2F transpose(vector2F& before)
{
	vector2F after(before[0].size(), vectorF(before.size(), 0));
	for(UINT i=0; i<before.size(); i++){
		for(UINT j=0; j<before[i].size(); j++) after[j][i] = before[i][j];
	}
	return after;
}

void pedFileReader::getGenotypes(vector3F& genos, vectorSS& snps, vectorS& geneNames, vector2F& mafs)
{
	genos.resize(0); geneNames.resize(0); mafs.resize(0);
	BOOST_FOREACH(ngsGene gene, __genes){
		vector2F oneGene;
                vectorS snpList;
		vectorF oneMaf;
		BOOST_FOREACH(pairSvF snp, gene.genoP){
			if(gene.maf.find(snp.first) == gene.maf.end()) {continue;}
			else{
				oneGene.push_back(snp.second);
                                snpList.push_back(snp.first);
				oneMaf.push_back(gene.maf[snp.first]);
			}
		}
		mafs.push_back(oneMaf);
                snps.push_back(snpList);
		geneNames.push_back(gene.geneName);
		genos.push_back(transpose(oneGene)); // need to transpose matrix
	}
	return;
}

void pedFileReader::getPhenotype(vectorS& inds, mapSvS& phenos)
{
	inds.resize(0); 
	phenos.clear();
	inds = __indOrder; 
	phenos = __phenotypes;
	//std::cout << "Dump out" << std::endl;
	//mapSvS::const_iterator it = phenos.begin();
	//while(it != phenos.end()){
		//std::cout << it->first << ": " << it->second << std::endl;
		//it++;
	//}

	return;
}

void pedFileReader::dumpOutGenos(std::string folder)
{
	std::string folderPrefix;
	if(makedir(folder)){
		folderPrefix = "./" +folder + "/";
	} else { std::cerr << "Error: can't create folder " << folder << std::endl; exit(-1);}

	BOOST_FOREACH(ngsGene gene, __genes){
		std::ofstream gout((folderPrefix + gene.geneName + ".tped").c_str()), mout((folderPrefix + gene.geneName + ".map").c_str());
		BOOST_FOREACH(pairSvF snp, gene.genoP){
			if(gene.maf.find(snp.first) == gene.maf.end()) {continue;}
			else{
				gout << snp.first << " " << snp.second << std::endl;
				mout << gene.geneName << " " << snp.first;
				//if(gene.maf[snp.first] != gene.maf[snp.first]) {;}
				if(gene.maf[snp.first] == -99) {mout << " NA";}
				else mout << " " << gene.maf[snp.first];
				mout << std::endl;
			}
		}
		gout.close(); mout.close();
	}
	return;
}
