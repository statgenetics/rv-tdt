#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <stdio.h>
#include <string>
#include <cstring>
#include <string>
#include <unistd.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_statistics.h>

#include <fstream>
#include <iostream>

#include "fb_skat.h"
#define SHOW(a) std::cout << #a << ": " << (a) << std::endl


using namespace std;

/* float *pass; */

int* vFtoHervI(vectorF& vf){
	UINT size = vf.size();
	int* res = new int[size + 10];

	for(UINT i=0; i<size; i++){
		res[i+1] = (int) vf[i];
	}
	return res;
}

float* vFtoHervF(vectorF& vf){
	UINT size = vf.size();
	float* res = new float[size + 10];

	for(UINT i=0; i< size; i++){
		res[i+1] = vf[i];
	}
	return res;
}


int number_of_lines (const char* filename){
    int nol = 0;
    ifstream myfile (filename);
    string line;
    if (filename){
		while(myfile.peek() != EOF){
			getline(myfile, line);
			nol++;
		}

		return nol;
    }
    else{ return -1;}

}


float fbSkat::get_skat_stat_fast(float **geno, float *Y, float *weight, int nbrmarkers, float offset)
{
	int i, j; //, k;
	float Qstat=0, s;
	if (ro==1) {
		for (i=1; i<=n; i++)
		{
			s=0;
			for (j=1; j<=nbrmarkers; j++)
				s=s+weight[j-1]*geno[i][j];
			Qstat=Qstat+(Y[i]-offset)*s;
		}
		return (Qstat*Qstat);
	}
	if (ro==0)
    {
		for (j=1; j<=nbrmarkers; j++)
        {
			s=0;
			for (i=1; i<=n; i++) {
				s=s+(Y[i]-offset)*geno[i][j];
				/* std::cout << geno[i][j] << " " << Y[i] << " | "; */
			}
			/* std::cout << std::endl; */
			/* SHOW(weight[j-i]); */
			/* SHOW(s); */
			Qstat=Qstat+weight[j-1]*weight[j-1]*s*s;
			/* SHOW(Qstat);  */
		}
		return (Qstat);
    }
	return(0);
}

int mendel_error(int p1, int p2, int ch)
{
	if ((p1==0) && (p2==0) && (ch==1)) return 1;
	if ((p1==2) && (p2==2) && (ch==1)) return 1;
	if ((p1==2) && (ch==0)) return 1;
	if ((p2==2) && (ch==0)) return 1;
	if ((p1==2) && (p2==2) && (ch==0)) return 1;
	if ((p1==0) && (ch==2)) return 1;
	if ((p2==0) && (ch==2)) return 1;
	if ((p1==0) && (p2==0) && (ch==2)) return 1;
	return 0;
}


inline float   abs1(float x)    {return (x>0.0 ? x :-x);}
inline int      abs1(int x)       {return (x>0 ? x :-x);}

float expected_geno_offspring(int ch, int p1, int p2)
{
	if ((p1 == p2) && (abs1(ch - p1) == 2)) return(0);
	else {
		if((p1 == 1) && (p2 == 1))
			return(pow(0.5,(abs1(ch - 1) + 1)));
		if(abs1(p1 - p2) == 2) {
			if(ch == 1)
				return(1);
			else return(0);
		}
		if((p1 == p2) && (p1 != 1)) {
			if(ch == p1)
				return(1);
			else return(0);
		}
		if(abs1(p1 - p2) == 1) {
			if(ch == p1)
				return(0.5);
			if(ch == p2)
				return(0.5);
			return(0);
		}
	}
	return(0);
}


fbSkat::fbSkat(std::string gene_name_, vector2F genotype, vectorF Y_, vectorF pass_, vectorF weight_ns_)
{
	// insert empty
	n = genotype.size()/3;
	x_genotype  =new int*[n+10];
	p1_genotype =new int*[n+10];
	p2_genotype =new int*[n+10];

	nall = genotype[0].size();

	for(int i=0; i < n; i++){

		x_genotype[i+1] = new int[nall+10];
		p1_genotype[i+1] = new int[nall+10];
		p2_genotype[i+1] = new int[nall+10];
		for (int j=0; j<nall; j++){
			p1_genotype[i+1][j+1] = genotype[3*i][j];
			p2_genotype[i+1][j+1] = genotype[3*i + 1][j];
			x_genotype[i+1][j+1] = genotype[3*i + 2][j];
		}

	    /* p1_genotype[i+1] = vFtoHervI(genotype[3*i]); */
	    /* p2_genotype[i+1] = vFtoHervI(genotype[3*i+1]); */
		/* x_genotype[i+1] = vFtoHervI(genotype[3*i+2]); */
	}

	/* for(int i=1; i<= n; i++){ */
	/* 	for(int j=1; j<=nall; j++){ */
	/* 		std::cout << x_genotype[i][j] << " "; */
	/* 	} */
	/* 	std::cout << std::endl; */
	/* } */

	Y = vFtoHervF(Y_);
	pass = vFtoHervF(pass_);
	weight_ns = vFtoHervF(weight_ns_);
	gene_name = gene_name_;
	return;
}

/**
 * Don't need geneFile, since method runs for each gene
 */

vectorF fbSkat::SKAT()
{
	/* int ig, b, p1, p2, k, m, i, j, nbrmarkers, nbrmarkers_temp, h, curr_mark,  nfam, min_mark, b1; */
	int ig, b, p1, p2, k, m, i, j, nbrmarkers, nbrmarkers_temp, curr_mark, b1, nfam;
	float Ymean;
	/* char a[100]; */
	int ngenes;
	/* float MAF, mend_th; */
	/* FILE *fout, *fin, *fmend; */
	/* float  *Y;              */
	/* Y =new float[10000]; */
	int *permut=new int[10000];
	const gsl_rng_type * T;
	gsl_rng * r;
	gsl_rng_env_setup();

	T = gsl_rng_default;
	r = gsl_rng_alloc (T);


	/* const char* pedigreeFile = argv[1]; // read pedigree file; 9 is missing data; 1 is affected; 0 unaffected; data.ped */
	/* const char* variantPassFile = argv[2]; // pass information for each variant; if this file is not available, set all to 1.; variant_pass.txt */
	/* const char* genesFile = argv[3]; // read gene file : gene names & start/end positions; genes.txt */
	/* const char* weightsFile = argv[4]; //read weights for variants; e.g. 1 for nonsynonymous and 0 otherwise; if this information is missing, set all weights to 1. weights.txt */
	/* const char* resultFile = argv[5];   //results_ */
	/* const char* mendelianErrorsFile = argv[6]; //mendelianErrorsFile; set to mendelian_errors_ */
	/* npermut=atoi(argv[7]); */
	/* MAF=atof(argv[8]); */
	/* mend_th=atof(argv[9]); */
	/* ro=atoi(argv[10]); */
	/* min_mark=atoi(argv[11]); */


	//./FB-SKAT.arg 144 411078 17609
	/* n=number_of_lines(pedigreeFile)/3; */
	/* nall=number_of_lines(weightsFile); */
	/* ngenes=number_of_lines(genesFile); */

	/* n = x_genotype. - 1; */
	/* nall = weight_ns.size() - 1; */
	ngenes = 1;

	/* cout << "FB-SKAT Version 2 " << endl; */
	/* cout << "Last updated Mar 18, 2013 " << endl; */
	/* cout << "Number of Variants = " << nall << endl; */
	/* cout << "Number of Regions = " << ngenes << endl; */
	/* cout << "Number of Trios = " << n << endl; */

	/* cout << "Pedigree file = " << pedigreeFile << endl; */
	/* cout << "Variant Pass file = " << variantPassFile << endl; */
	/* cout << "Gene file = " << genesFile << endl; */
	/* cout << "Weight file = " << weightsFile << endl; */
	/* cout << "Result File = " << resultFile << endl; */
	/* cout << "Mendelian Error File = " << mendelianErrorsFile << endl; */

	/////


	// alloc_data
	mend=new float[nall+10];
	weight=new float[nall+10];
	/* weight_ns=new float[nall+10]; */
	/* pass=new float[nall+10]; */
	gene=new char*[ngenes+10];
	for (i=1; i<=ngenes; i++)
		gene[i]=new char[100];
	start_gene=new int[ngenes+10];
	end_gene=new int[ngenes+10];

	Qstat_p=new float[npermut+10];


	// read pedigree file; 9 is missing data; 1 is affected; 0 unaffected
	/* fin=fopen(pedigreeFile,"r"); */
	/* for (i=1; i<=n; i++) */
	/* { */
	/* 	char b[100]; */
	/* 	fscanf(fin,"%s %s %s %s %s %f",a,b,a,a,a,&Y[i]); */
	/* 	for (j=1; j<=nall; j++) */
	/* 		fscanf(fin,"%d",&x_genotype[i][j]); */
	/* 	fscanf(fin,"%s %s %s %s %s %s",a,a,a,a,a,a); */
	/* 	for (j=1; j<=nall; j++) */
	/* 		fscanf(fin,"%d",&p1_genotype[i][j]); */
	/* 	fscanf(fin,"%s %s %s %s %s %s",a,a,a,a,a,a); */
	/* 	for (j=1; j<=nall; j++) */
	/* 		fscanf(fin,"%d",&p2_genotype[i][j]); */
	/* } */
	/* fclose(fin); */


	// pass information for each variant; if this file is not available, set all to 1.
	/* fin=fopen(variantPassFile,"r"); */
	/* for (i=1; i<=nall; i++) */
	/* 	fscanf(fin,"%f",&pass[i]); */
	/* fclose(fin); */

	/* std::cout << "pass:\n"; */
	/* for (i=1; i<=nall; i++) */
	/* 	std::cout << pass[i] << " "; */
	/* std::cout << std::endl; */

	//Mendelian errors
	//remove marker j if too many mendelian errors
	/* char numefile[200]; */
	/* char buffer[100]; */

	/* strcpy(numefile,mendelianErrorsFile); */
	/* sprintf(buffer,"%f",MAF); */
	/* strcat(numefile,buffer); */
	/* sprintf(buffer,"%s","_"); */
	/* strcat(numefile,buffer); */
	/* sprintf(buffer,"%f",mend_th); */
	/* strcat(numefile,buffer); */
	/* sprintf(buffer,"%s","_"); */
	/* strcat(numefile,buffer); */
	/* sprintf(buffer,"%d",ro); */
	/* strcat(numefile,buffer); */
	/* sprintf(buffer,"%s",".txt"); */
	/* strcat(numefile,buffer); */

	/* fmend=fopen(numefile,"w"); */

	for (j=1; j<=nall; j++)
	{
		mend[j]=0; 
		nfam=0;
		for (i=1; i<=n; i++)
	    {
			if (x_genotype[i][j]+p1_genotype[i][j]+p2_genotype[i][j]<9) {
				nfam=nfam+1;
				if (mendel_error(p1_genotype[i][j],p2_genotype[i][j],x_genotype[i][j])>0) {
					p1_genotype[i][j]=0; 
					p2_genotype[i][j]=0; 
					x_genotype[i][j]=0; 
					mend[j]=mend[j]+1;
				}
			}
	    }
		/* int ili=(int)mend[j]; */
		if (nfam>0) mend[j]=mend[j]/(1.0*nfam);
		/* fprintf(fmend,"For marker %d there are %d (%f) mendel errors\n",j,ili,mend[j]); */
	}
	/* fclose(fmend); */

	for (i=1; i<=n; i++)
		for (j=1; j<=nall; j++)
			if (x_genotype[i][j]+p1_genotype[i][j]+p2_genotype[i][j]>=9) {
				p1_genotype[i][j]=0; p2_genotype[i][j]=0; x_genotype[i][j]=0;
			}

	/* std::cout << "geno:\n"; */
	/* for (i=1; i<=nall; i++) */
	/* 	std::cout << x_genotype[1][i] << " "; */
	/* std::cout << std::endl; */
	/* for (i=1; i<=nall; i++) */
	/* 	std::cout << p1_genotype[1][i] << " "; */
	/* std::cout << std::endl; */
	/* for (i=1; i<=nall; i++) */
	/* 	std::cout << p2_genotype[1][i] << " "; */
	/* std::cout << std::endl; */

	// read gene file : gene names & start/end positions

	/* fin=fopen(genesFile,"r"); */
	/* for (i=1; i<=ngenes; i++) */
	/* 	fscanf(fin,"%s %d %d",gene[i],&start_gene[i],&end_gene[i]); */
	/* fclose(fin); */
	char * gname = new char[gene_name.size() + 1];
	strcpy(gname, gene_name.c_str());
	gene[1] = gname;
	start_gene[1] = 1;
	end_gene[1] = nall;

	// read weights for variants; e.g. 1 for nonsynonymous and 0 otherwise; if this information is missing, set all weights to 1.

	/* fin=fopen(weightsFile,"r"); */
	/* for (i=1; i<=nall; i++) */
	/* 	fscanf(fin,"%f",&weight_ns[i]); */
	/* fclose(fin); */

	/* std::cout << "weight:\n"; */
	/* for (i=1; i<=nall; i++) */
	/* 	std::cout << weight_ns[i] << " "; */
	/* std::cout << std::endl; */

	float   Echp1p2[10][10];
	for(p1=0;p1<=2;p1++)
		for(p2=0;p2<=2;p2++)
			Echp1p2[p1][p2]=0.0;

	for(p1=0;p1<=2;p1++)
		for(p2=0;p2<=2;p2++)
			for(m=0;m<=2;m++)
				Echp1p2[p1][p2]+=m*expected_geno_offspring(m,p1,p2);



	//write result
	/* char numfile2[200]; */
	/* char buffer2[100]; */
	/* strcpy(numfile2,resultFile); */
	/* sprintf(buffer2,"%f",MAF); */
	/* strcat(numfile2,buffer2); */
	/* sprintf(buffer2,"%s","_"); */
	/* strcat(numfile2,buffer2); */
	/* sprintf(buffer2,"%f",mend_th); */
	/* strcat(numfile2,buffer2); */
	/* sprintf(buffer2,"%s","_"); */
	/* strcat(numfile2,buffer2); */
	/* sprintf(buffer2,"%d",ro); */
	/* strcat(numfile2,buffer2); */
	/* sprintf(buffer2,"%s",".txt"); */
	/* strcat(numfile2,buffer2); */

	/* fout=fopen(numfile2,"w"); */
	/* for (ig=1; ig<=ngenes; ig++) */
	/* { */
	ig = 1;
	start=start_gene[ig];
	end=end_gene[ig];
	nbrmarkers_temp=end-start+1;
	zz=new int[nbrmarkers_temp+10];
	pf=new float[nbrmarkers_temp+10];
	pf_temp=new float[nbrmarkers_temp+10];
	geno_temp=new float*[n+10];
	for (i=1; i<=n; i++)
		geno_temp[i]=new float[nbrmarkers_temp+10];

	for (j=1; j<=nbrmarkers_temp; j++)
	{ zz[j]=0; pf_temp[j]=0;}

	for (i=1; i<=n; i++)
	{
		for (j=start; j<=(start+nbrmarkers_temp-1); j++)
		{
			geno_temp[i][j-start+1]=x_genotype[i][j]- Echp1p2[p1_genotype[i][j]][p2_genotype[i][j]];
			pf_temp[j-start+1]=pf_temp[j-start+1]+p1_genotype[i][j]+p2_genotype[i][j];
			zz[j-start+1]=zz[j-start+1]+x_genotype[i][j];
		}
	}

	for (j=start; j<=(start+nbrmarkers_temp-1); j++)
		pf_temp[j-start+1]=(pf_temp[j-start+1]+1.0)/(4.0*n+1.0);

	geno=new float*[n+10];
	for (i=1; i<=n; i++)
		geno[i]=new float[nbrmarkers_temp+10];
	geno_p=new float*[n+10];
	for (i=1; i<=n; i++)
		geno_p[i]=new float[nbrmarkers_temp+10];

	curr_mark=0;
	/* float stat2=0; */
	/* std::cout << start << " " << nbrmarkers_temp << std::endl; */
	for (j=start; j<=(start+nbrmarkers_temp-1); j++)
	{
		// weight_ns[j]>0 OR weight_ns[j]>=0 ???
		if ((zz[j-start+1]+pf_temp[j-start+1]>0) && (pf_temp[j-start+1]<=MAF) && (mend[j]<=mend_th) && (pass[j]>0) && (weight_ns[j]>0)) {
			//if ((zz[j-start+1]+pf_temp[j-start+1]>0) && (pf_temp[j-start+1]<=MAF) && (mend[j]<=mend_th) && (pass[j]>0) && (weight_ns[j]>=0)) {
			/* std::cout << "\tanalyzed"; */
			curr_mark=curr_mark+1;
			pf[curr_mark]=pf_temp[j-start+1];
			for (i=1; i<=n; i++)
				geno[i][curr_mark]=geno_temp[i][j-start+1];
		}
		else {
			// invalid site
			/* std::cout << j << " " << zz[j-start+1]+pf_temp[j-start+1] << " " << pf_temp[j-start+1] << " " << mend[j] << " " << pass[j] << " " << weight_ns[j] << std::endl; */
		}
	}

	nbrmarkers=curr_mark;
	//std::cout << "nbrmarkers: " << nbrmarkers << std::endl;
	vectorF res_wrong;
	for(UINT x = 0; x < 4; x++) 
		res_wrong.push_back(NAN);
	if (nbrmarkers<min_mark) {
		return res_wrong;	
	}

	//can also have a weight that is data-dependent
	for (k=0; k<nbrmarkers; k++)
		weight[k]=gsl_ran_beta_pdf(pf[k+1],1,25);

	/* for (i=0; i<nbrmarkers; i++) */
	/* 	std::cout << weight[i] << " "; */
	/* std::cout << std::endl; */

	Ymean=0;
		
	for (i=0; i<n; i++){
		/* std::cout << Y[i+1] << " "; */
		Ymean=Ymean+Y[i+1];
	}
	//if all affected
	if ((Ymean>=n) && (Ymean<=n))  Ymean=0.05;
	else Ymean=Ymean/n;

	/* std::cout << geno[1][1] << " " << Ymean << std::endl; */

	float Qstat=get_skat_stat_fast(geno,Y,weight,nbrmarkers,Ymean);

	/* std::cout << Qstat << std::endl; */

	for (b1=1; b1<=npermut; b1++)
	{
		for (i=1; i<=n; i++)
			permut[i]=2*gsl_ran_binomial(r,0.5,1)-1;
		
		for (i=0; i<n; i++)
			for (k=0; k<nbrmarkers; k++)
				geno_p[i+1][k+1]=permut[i+1]*geno[i+1][k+1];

		Qstat_p[b1]=get_skat_stat_fast(geno_p,Y,weight,nbrmarkers,Ymean);
	}

	float mu_Q=0;
	for (b=1; b<=npermut; b++)
		mu_Q=mu_Q+Qstat_p[b];
	mu_Q=mu_Q/npermut;
	float var_Q=0;
	for (b=1; b<=npermut; b++)
		var_Q=var_Q+(Qstat_p[b]-mu_Q)*(Qstat_p[b]-mu_Q);
	var_Q=var_Q/npermut;
	float mu4_Q=0;
	for (b=1; b<=npermut; b++)
		mu4_Q=mu4_Q+(Qstat_p[b]-mu_Q)*(Qstat_p[b]-mu_Q)*(Qstat_p[b]-mu_Q)*(Qstat_p[b]-mu_Q);
	mu4_Q=mu4_Q/npermut;
	float gamma=mu4_Q/(var_Q*var_Q)-3;
	float newdf=12.0/gamma;

	/* SHOW(npermut); */
	/* SHOW(gamma); */
	/* SHOW(Qstat); */
	/* SHOW(var_Q); */
	/* SHOW(mu_Q); */
	/* SHOW(newdf); */
	/* std::cout << std::endl; */
	/* std::cout << "\n" << gene[ig] << " " << start_gene[ig] << " " << end_gene[ig] << " " << Qstat*2*mu_Q/var_Q << " " << 2*mu_Q*mu_Q/var_Q << " " << (Qstat-mu_Q)*sqrt(2*newdf)/sqrt(var_Q)+newdf << " " << newdf << std::endl; */

	// TODO: why this would happen
	if (abs(var_Q) < 1.0e-6) {
		return res_wrong;
	}

	vectorF res;
	res.push_back( Qstat*2*mu_Q/var_Q );
	res.push_back( 2*mu_Q*mu_Q/var_Q );
	res.push_back( (Qstat-mu_Q)*sqrt(2*newdf)/sqrt(var_Q)+newdf );
	res.push_back( newdf);
	return res;
/*	}  */
}
