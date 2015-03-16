/**
 * \file        fb_skat.h
 * \author      Zongxiao He (zongxiah@bcm.edu)
 * \copyright   Copyright 2013, Leal Group
 * \date        2013-07-10
 *
 * \brief
 *
 *
 */


#ifndef FB_SKAT_H
#define FB_SKAT_H

#include <vector>
#include <string>


typedef std::vector<double> vectorF;
typedef std::vector<int> vectorI;
typedef std::vector< std::vector<double> > vector2F;
typedef std::vector< std::vector <int> > vector2I;
typedef unsigned int UINT;

class fbSkat
{
public:
	fbSkat(std::string gene_name, vector2F genotype, vectorF Y, vectorF pass, vectorF weight_ns);
	void TestLoad(int nper_, float MAF_, float mend_th_, int ro_, int min_mark_){
		npermut = nper_;
		MAF = MAF_;
		mend_th = mend_th_;
		ro = ro_;
		min_mark = min_mark_;
		return;
	}
	vectorF SKAT();

private:
		int npermut;
		float MAF;
		float mend_th;
		int ro;
		int min_mark;

		// insert a empty element at first place
		/* vector2I x_genotype, p1_genotype, p2_genotype; */
		/* vectorF Y; */
		/* vectorF pass, weight_ns; */
		std::string gene_name;

		float get_skat_stat_fast(float **geno, float *Y, float *weight, int nbrmarkers, float offset);

		// variable as global
		float *pass;
		char **gene;
		float *mend;
		float *Qstat_p;
		float **geno, **geno_temp, **geno_p;
		int **x_genotype, **p1_genotype, **p2_genotype;
		float *weight, *weight_ns, *Y;
		float *pf, *pf_temp;
		int start, end;
		int *start_gene, *end_gene;
		int *zz;
		int  n, nall;
};

int number_of_lines (const char* filename);
int mendel_error(int p1, int p2, int ch);
float expected_geno_offspring(int ch, int p1, int p2);

int* vFtoHervI(vectorF& vf);
float* vFtoHervF(vectorF& vf);

#endif /* FB_SKAT_H */
