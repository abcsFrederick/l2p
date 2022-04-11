
#include <ctype.h>
#include <errno.h>
#include <limits.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "pathworks.h"



static inline double lngamm(const double z)
// Reference: "Lanczos, C. 'A precision approximation
// of the gamma double ', J. SIAM Numer. Anal., B, 1, 86-96, 1964."
// Translation of  Alan Miller's FORTRAN-implementation
// See http://lib.stat.cmu.edu/apstat/245
{
  double x = 0;
  x += 0.1659470187408462e-06/(z+7);
  x += 0.9934937113930748e-05/(z+6);
  x -= 0.1385710331296526    /(z+5);
  x += 12.50734324009056     /(z+4);
  x -= 176.6150291498386     /(z+3);
  x += 771.3234287757674     /(z+2);
  x -= 1259.139216722289     /(z+1);
  x += 676.5203681218835     /(z);
  x += 0.9999999999995183;
//   return(Math.log(x)-5.58106146679532777-z+(z-0.5)*Math.log(z+6.5));
   return(log(x)-5.58106146679532777-z+(z-0.5)*log(z+6.5));
}

static inline double lnfact(const int n)
{
  if(n<=1) return(0);
  return(lngamm(n+1));
}

static inline double lnbico(const int n,const int k)
{
  double ret;
  ret =lnfact(n)-lnfact(k)-lnfact(n-k);

  return(ret);
}

static inline double hyper_323(const int n11,const int n1_,const int n_1,int n)
{
    double d;
    d = lnbico(n1_,n11) + lnbico(n-n1_,n_1-n11) - lnbico(n,n_1);
    d = exp(d);
    return d;
}

static int sn11,sn1_,sn_1,sn;
static double sprob;

static double hyper0(int n11i,int n1_i,int n_1i,int ni)
{

// printf("in hyper0 %d %d %d %d \n",n11i,n1_i,n_1i,ni);
  if(!(n1_i|n_1i|ni))
  {
//printf("in hyper0 NOT %d %d %d %d \n",n11i,n1_i,n_1i,ni);
    if(!(n11i % 10 == 0))
    {
      if(n11i==sn11+1)
      {
        sprob *= ((double)(sn1_-sn11)/(double)(n11i))*((double)(sn_1-sn11)/(double)(n11i+sn-sn1_-sn_1));
        sn11 = n11i;
        return sprob;
      }
      if(n11i==sn11-1)
      {
        sprob *= ((double)(sn11)/(double)(sn1_-n11i))*((double)(sn11+sn-sn1_-sn_1)/(double)(sn_1-n11i));
        sn11 = n11i;
        return sprob;
      }
    }
    sn11 = n11i;
  }
  else
  {
//printf("in hyper0 else %d %d %d %d \n",n11i,n1_i,n_1i,ni);
    sn11 = n11i;
    sn1_=n1_i;
    sn_1=n_1i;
    sn=ni;
  }
// printf("in hyper0 before hyper_323 %d %d %d %d\n",sn11,sn1_,sn_1,sn);
  sprob = hyper_323(sn11,sn1_,sn_1,sn);
// printf("hyper returns sprob = %10.15f after hyper_323\n",sprob);
  return sprob;
}

static double  hyper(const int n11)
{
  return(hyper0(n11,0,0,0));
}

static double sleft,sright,sless,slarg;

static double exact(const int n11,const int n1_,const int n_1,const int n)
{
  int i,j;
  double p;
  double prob;
  int max=n1_;

// printf("in exact %d %d %d %d \n",n11,n1_, n_1,n);
  if(n_1<max) max=n_1;
  int min = n1_+n_1-n;
  if(min<0) min=0;
  if(min==max)
  {
    sless = sright= sleft = slarg = 1;
    return 1;
  }

  prob=hyper0(n11,n1_,n_1,n);
// printf("in exact prob=%20.17f \n",prob);
  sleft=0;
  p=hyper(min);
  for(i=min+1; p<0.99999999*prob; i++)
  {
    sleft += p;
    p=hyper(i);
  }
  i--;
  if(p<1.00000001*prob) sleft += p;
  else i--;
  sright=0;
  p=hyper(max);
  for(j=max-1; p<0.99999999*prob; j--)
  {
    sright += p;
    p=hyper(j);
  }
  j++;
  if(p<1.00000001*prob) sright += p;
  else j++;
  if(abs(i-n11)<abs(j-n11))
  {
    sless = sleft;
    slarg = 1 - sleft + prob;
  }
  else
  {
    sless = 1 - sright + prob;
    slarg = sright;
  }
  return prob;
}


static double left,right,twotail;
//static int n11old=-1;
//static int n12old=-1;
//static int n21old=-1;
//static int n22old=-1;

double exact22(int n11_,int n12_,int n21_,int n22_)
{
#if 0
  double prob = 0.0;
  double n11_ = parseInt("0"+n11,10);
  double n12_ = parseInt("0"+n12,10);
  double n21_ = parseInt("0"+n21,10);
  double n22_ = parseInt("0"+n22,10);
   if((n11old==n11_i) && (n12old==n12_) && (n21old==n21_) && (n22old==n22_)) return;
  n11old=n11_;
  n12old=n12_;
  n21old=n21_;
  n22old=n22_;
#endif
  if(n11_<0) n11_ *= -1;
  if(n12_<0) n12_ *= -1;
  if(n21_<0) n21_ *= -1;
  if(n22_<0) n22_ *= -1;

  int n1_ = n11_+n12_;
  int n_1 = n11_+n21_;
  int n   = n11_ +n12_ +n21_ +n22_;

  // prob = exact(n11_,n1_,n_1,n); don't need return value
  ( void ) exact(n11_,n1_,n_1,n);
// printf("prob after exact is  %30.25f\n",prob);

  left    = sless;
  right   = slarg;
  twotail = sleft+sright;
// printf("left=%20.15f right=%20.15f \n",sleft,sright);
  if(twotail>1) twotail=1;

  return twotail;
// printf("%d %d %d %d %12.8f prob=%20.15f twotail=%20.15f\n",(int)n11_, (int)n12_, (int)n21_, (int)n22_, twotail,prob,twotail);
/*
  document.form1.output.value +=
  newline+
  " TABLE = [ " +
  n11_+" , "+
  n12_+" , "+
  n21_+" , "+
  n22_+" ]" + newline +
  "Left   : p-value = "+ left + newline +
  "Right  : p-value = "+ right + newline +
  "2-Tail : p-value = "+ twotail +
  newline +   "------------------------------------------";
*/
}

// test oneside
static double exact_oneside(int n11,int n1_,int n_1,int n)
{
  int i,j;
  double p;
  double prob;
  int max=n1_;

// printf("in exact %d %d %d %d \n",n11,n1_, n_1,n);
  if(n_1<max) max=n_1;
  int min = n1_+n_1-n;
  if(min<0) min=0;
  if(min==max)
  {
    sless = sright= sleft = slarg = 1;
    return 1;
  }

  prob=hyper0(n11,n1_,n_1,n);
// printf("in exact prob=%20.17f \n",prob);
  sleft=0;
  p=hyper(min);
  for(i=min+1; p<0.99999999*prob; i++)
  {
    sleft += p;
    p=hyper(i);
  }
  i--;
  if(p<1.00000001*prob) sleft += p;
  else i--;
  sright=0;
  p=hyper(max);
  for(j=max-1; p<0.99999999*prob; j--)
  {
    sright += p;
    p=hyper(j);
  }
  j++;
  if(p<1.00000001*prob) sright += p;
  else j++;
  if(abs(i-n11)<abs(j-n11))
  {
    sless = sleft;
    slarg = 1 - sleft + prob;
  }
  else
  {
    sless = 1 - sright + prob;
    slarg = sright;
  }
  return prob;
}

double exact22_oneside(int n11_,int n12_,int n21_,int n22_, int dbg)
{
#if 0
  double prob = 0.0;
  double n11_ = parseInt("0"+n11,10);
  double n12_ = parseInt("0"+n12,10);
  double n21_ = parseInt("0"+n21,10);
  double n22_ = parseInt("0"+n22,10);
   if((n11old==n11_i) && (n12old==n12_) && (n21old==n21_) && (n22old==n22_)) return;
  n11old=n11_;
  n12old=n12_;
  n21old=n21_;
  n22old=n22_;
#endif
  if(n11_<0) n11_ *= -1;
  if(n12_<0) n12_ *= -1;
  if(n21_<0) n21_ *= -1;
  if(n22_<0) n22_ *= -1;

  int n1_ = n11_+n12_;
  int n_1 = n11_+n21_;
  int n   = n11_ +n12_ +n21_ +n22_;

  // prob = exact(n11_,n1_,n_1,n); don't need return value
  ( void ) exact_oneside(n11_,n1_,n_1,n);
// printf("prob after exact is  %30.25f\n",prob);

  left    = sless;
  right   = slarg;
  twotail = sleft+sright;
#if 0
  if(twotail>1) twotail=1;
  return twotail;
#else
if (dbg)
{
fprintf(stderr,"sleft=%20.15f sright=%20.15f \n",sleft,sright);
}
  if(sright>1) sright=1;
  return sright;
#endif
// printf("%d %d %d %d %12.8f prob=%20.15f twotail=%20.15f\n",(int)n11_, (int)n12_, (int)n21_, (int)n22_, twotail,prob,twotail);
/*
  document.form1.output.value +=
  newline+
  " TABLE = [ " +
  n11_+" , "+
  n12_+" , "+
  n21_+" , "+
  n22_+" ]" + newline +
  "Left   : p-value = "+ left + newline +
  "Right  : p-value = "+ right + newline +
  "2-Tail : p-value = "+ twotail +
  newline +   "------------------------------------------";
*/
}

