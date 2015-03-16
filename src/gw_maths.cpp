// =====================================================================================
// 
//       Filename:  gw_maths.cpp
// 
//    Description:  useful maths functions
// 
//        Version:  1.0
//        Created:  01/04/2011 10:49:19 PM
//       Revision:  none
//       Compiler:  g++
// 
//         Author:  Gao Wang (gw), wangow@gmail.com
//                  Baylor College of Medicine, Texas, USA
//        License:  GNU General Public License <http://www.gnu.org/licenses/>
//                  Copyright (c) 2011, Gao Wang
// 
// =====================================================================================


#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <functional>
#include <cmath>
#include <cstdlib>
#include <limits>
#include "gw_maths.h"
#include "fisher2.h"

//!- Basic

double gw_dmax(double a, double b) 
{
  if (a < b) return b;
  else return a;
}


int gw_imax(int a, int b) 
{
  if (a < b) return b;
  else return a;
}


double gw_dmin(double a, double b) 
{
  if (a < b) return a;
  else return b;
}


int gw_imin(int a, int b) 
{
  if (a < b) return a;
  else return b;
}


void gw_round(double& myValue, double PRECISION) 
{
  double myRemainder = fmod(myValue, PRECISION);
  if (myRemainder > ( 0.5 * PRECISION ))
  {
    myValue = (myValue - myRemainder + PRECISION);
  }
  else
  {
    myValue = (myValue  - myRemainder);
  }
  return;
}

bool fEqual(double a, double b) {
  return std::fabs(a - b) < std::numeric_limits<double>::epsilon();
}

//!- Statistics

double gw_sum(const std::vector<double>& mydata) 
{
  double dataSum(0.0);
  for (unsigned int i = 0; i != mydata.size(); ++i) 
    dataSum = dataSum + mydata[i];
  return dataSum;
}

double gw_prod(const std::vector<double>& mydata) 
{
  double dataProd(1.0);
  for (unsigned int i = 0; i != mydata.size(); ++i) 
    dataProd = dataProd * mydata[i];
  return dataProd;
}


double gw_mean(const std::vector<double>& mydata) 
{
  double dataSum(0.0);
  for (unsigned int i = 0; i != mydata.size(); ++i) 
    dataSum = dataSum + mydata[i];
  return dataSum / double(mydata.size());
}


double gw_var(const std::vector<double>& mydata) 
{
  double mMean = gw_mean(mydata);
  double dataSum(0.0);
  for (unsigned int i = 0; i != mydata.size(); ++i) 
    dataSum = dataSum + pow(mydata[i] - mMean, 2.0);
  return dataSum / double( mydata.size() - 1 );
}


double gw_sd(const std::vector<double>& mydata) 
{
  double mVar = gw_var(mydata);
  return sqrt(mVar);
}


//!- Factorials, combinatory and hypergeometric pmf/cmf
/*
#include "lnfact_table"
double gw_ln_factorial(double n) 
{
	
//    double rlnFact = 0.0;
//	if (n == 0 || n == 1) return rlnFact;
//    for (int i = 2; i <= n; ++i) rlnFact += log(double(i));
//	return rlnFact;


//  return lnfact_table[int(n)];
}


double gw_lnchoose(double n, double k) 
{
    
	double lnChoose = 0.0;
	if (n == 0 || n == k || k == 0) 
    return lnChoose;
	lnChoose = gw_ln_factorial(n) - gw_ln_factorial(n-k) - gw_ln_factorial(k);
	return lnChoose;
}

double gw_hypergeometric_pmf(unsigned int k, unsigned int n1, unsigned int n2, unsigned int t) 
{
  // prob(k) =  choose(n1, k) choose(n2, t-k) / choose(n1+n2,t)

  
     if (k > n1 || k > t || t - k > n2 || t > n1 + n2) {
     cout << "error!" << endl; 
     exit(-1);
     }
     
  double dn1 = double(n1), dk = double(k), dn2 = double(n2), dt = double(t);

  double c1 = gw_lnchoose(dn1,dk);
  double c2 = gw_lnchoose(dn2,dt-dk);
  double c3 = gw_lnchoose(dn1+dn2,dt);
  double hPmf = exp(c1 + c2 - c3);

  return hPmf;
}


double gw_hypergeometric_cmf(unsigned int k, unsigned int n1, unsigned int n2, unsigned int t) 
{
  // prob(<= k) = \sum_0^k choose(n1, k) choose(n2, t-k) / choose(n1+n2,t)

    
//        if (k > n1 || k > t || t - k > n2 || t > n1 + n2) {
//        cout << "error!" << endl; 
//        exit(-1);
//        }
  double hCmf = gw_hypergeometric_pmf(k, n1, n2, t);
  double hPmf = hCmf;
  for (double i = k; i >= 0; --i) {
    double diffactor = i / (n1 - i + 1.0) * (n2 - t + i) / (t - i + 1.0);
    hPmf = hPmf * diffactor;
    hCmf += hPmf;
  }

  return hCmf;
}

*/


//!- Fisher's test for 2 by 2 tables

double fexact_two_sided_pvalue(const std::vector<int> &twotwoTable) 
{
  // This is specific for 2x2 tables. contingency_table = matrix(twotwoTable, 2, 2, byrow = T)

  double contingency_table[4] = {0, 0, 0, 0};
  for (int i = 0; i != 4; ++i) contingency_table[i] = twotwoTable[i];
  //stuff for Fisher's test
  int nrow = 2;
  int ncol = 2;
  double expected = -1.0;
  double percnt = 100.0;
  double emin = 0.0;
  double prt = 0.0;
  double pvalue = 0.0;
  int workspace = 300000;

  fexact(&nrow, &ncol, contingency_table, &nrow, &expected, &percnt, &emin, &prt, &pvalue, &workspace);
  return pvalue;
}



//!- Mann-Whitney rank test statistic "U"
// www.alglib.net
// http://en.wikipedia.org/wiki/Mann-Whitney_U
double Mann_Whitneyu(double x[], int n, double y[], int m)
{
  int i;
  int k;
  int t;
  double tmp;
  int tmpi;
  int ns;
  double u;
  int w;

  // Prepare
  
  if( n <= 5 || m <= 5 )
  {
    std::cout << "Sample size too small" << std::endl;
    exit(1);
  }
  ns = n+m;
  double r[ns-1];
  int c[ns-1];
  for(i = 0; i <= n-1; i++)
  {
    r[i] = x[i];
    c[i] = 0;
  }
  for(i = 0; i <= m-1; i++)
  {
    r[n+i] = y[i];
    c[n+i] = 1;
  }

  // sort {R, C}, QS: smaller scores ranks higher
  
  if( ns!=1 )
  {
    i = 2;
    do
    {
      t = i;
      while(t!=1)
      {
        k = t/2;
        if( r[k-1]>=r[t-1] )
        {
          t = 1;
        }
        else
        {
          tmp = r[k-1];
          r[k-1] = r[t-1];
          r[t-1] = tmp;
          tmpi = c[k-1];
          c[k-1] = c[t-1];
          c[t-1] = tmpi;
          t = k;
        }
      }
      i = i+1;
    }
    while(i<=ns);
    i = ns-1;
    do
    {
      tmp = r[i];
      r[i] = r[0];
      r[0] = tmp;
      tmpi = c[i];
      c[i] = c[0];
      c[0] = tmpi;
      t = 1;
      while(t!=0)
      {
        k = 2*t;
        if( k>i )
        {
          t = 0;
        }
        else
        {
          if( k<i )
          {
            if( r[k] > r[k-1] )
            {
              k = k+1;
            }
          }
          if( r[t-1] >= r[k-1] )
          {
            t = 0;
          }
          else
          {
            tmp = r[k-1];
            r[k-1] = r[t-1];
            r[t-1] = tmp;
            tmpi = c[k-1];
            c[k-1] = c[t-1];
            c[t-1] = tmpi;
            t = k;
          }
        }
      }
      i = i-1;
    }
    while(i>=1);
  }

  // Compute U
  
  u = 0;
  w = 1;
  for(i = 0; i <= ns-1; i++)
  {
    if (i == 0)
      w = 1;
    else {
      if (r[i] > r[i-1])
        w = i + 1;
    }

    if(c[i] == 0)
    {
      //ranks (sum of) for cases
      //std::cout << w << " ";
      //std::cout << r[i] << " ";
      u = u+w;
    }
  }

  //u = n*m+m*(m+1)/2-u;
  //std::cout << u << " ";
  return u;
}
///:~
