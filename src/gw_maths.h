// =====================================================================================
// 
//       Filename:  gw_maths.h
// 
//    Description:  useful maths functions
// 
//        Version:  1.0
//        Created:  01/04/2011 10:49:48 PM
//       Revision:  none
//       Compiler:  g++
// 
//         Author:  Gao Wang (gw), wangow@gmail.com
//                  Baylor College of Medicine, Texas, USA
//        License:  GNU General Public License <http://www.gnu.org/licenses/>
//                  Copyright (c) 2011, Gao Wang
// 
// =====================================================================================


#ifndef GWMATHS_H
#define GWMATHS_H

#include<vector>

double gw_dmax(double , double ); 
int gw_imax(int , int );
double gw_dmin(double , double );
int gw_imin(int , int );
void gw_round(double&myValue, double PRECISION);
bool fEqual(double a, double b);

//!- Statistics

double gw_sum(const std::vector<double>&); 
double gw_prod(const std::vector<double>&); 
double gw_mean(const std::vector<double>&); 
double gw_var(const std::vector<double>&); 
double gw_sd(const std::vector<double>&);

//extern double lnfact_table[]; 
//double gw_ln_factorial(double );
//double gw_lnchoose(double, double );
//double gw_hypergeometric_pmf(unsigned int , unsigned int , unsigned int , unsigned int );
//double gw_hypergeometric_cmf(unsigned int , unsigned int , unsigned int , unsigned int );

double fexact_two_sided_pvalue(const std::vector<int>&);
double Mann_Whitneyu(double*, int , double*, int );
#endif
///:~
