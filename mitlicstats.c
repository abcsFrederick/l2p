
/*
vi mitlicstats.c ; gcc -Werror -o mitlicstats mitlicstats.c -lm
fisher exact code from : https://github.com/gatoravi/fisher-exact
benjamini-hochberg code from : https://rosettagit.org/drafts/p-value-correction/
*/

#include "stdio.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
// #include "kfunc.h"


/* The MIT License

   Copyright (C) 2010, 2013 Genome Research Ltd.
   Copyright (C) 2011 Attractive Chaos <attractor@live.co.uk>

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/


#if 0
#include "kfunc.h"
#else
/* The MIT License

   Copyright (C) 2010, 2013 Genome Research Ltd.
   Copyright (C) 2011 Attractive Chaos <attractor@live.co.uk>

   Permission is hereby granted, free of charge, to any person obtaining
   a copy of this software and associated documentation files (the
   "Software"), to deal in the Software without restriction, including
   without limitation the rights to use, copy, modify, merge, publish,
   distribute, sublicense, and/or sell copies of the Software, and to
   permit persons to whom the Software is furnished to do so, subject to
   the following conditions:

   The above copyright notice and this permission notice shall be
   included in all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
   EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
   MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
   NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
   BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
   CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
*/

#ifndef HTSLIB_KFUNC_H
#define HTSLIB_KFUNC_H

#ifdef __cplusplus
extern "C" {
#endif

/* Log gamma function
 * \log{\Gamma(z)}
 * AS245, 2nd algorithm, http://lib.stat.cmu.edu/apstat/245
 */
static double kf_lgamma(double z);

/* complementary error function
 * \frac{2}{\sqrt{\pi}} \int_x^{\infty} e^{-t^2} dt
 * AS66, 2nd algorithm, http://lib.stat.cmu.edu/apstat/66
 */
double kf_erfc(double x);

/* The following computes regularized incomplete gamma functions.
 * Formulas are taken from Wiki, with additional input from Numerical
 * Recipes in C (for modified Lentz's algorithm) and AS245
 * (http://lib.stat.cmu.edu/apstat/245).
 *
 * A good online calculator is available at:
 *
 *   http://www.danielsoper.com/statcalc/calc23.aspx
 *
 * It calculates upper incomplete gamma function, which equals
 * kf_gammaq(s,z)*tgamma(s).
 */

double kf_gammap(double s, double z);
double kf_gammaq(double s, double z);

/* Regularized incomplete beta function. The method is taken from
 * Numerical Recipe in C, 2nd edition, section 6.4. The following web
 * page calculates the incomplete beta function, which equals
 * kf_betai(a,b,x) * gamma(a) * gamma(b) / gamma(a+b):
 *
 *   http://www.danielsoper.com/statcalc/calc36.aspx
 */
// static double kf_betai(double a, double b, double x);

/*
 *    n11  n12  | n1_
 *    n21  n22  | n2_
 *   -----------+----
 *    n_1  n_2  | n
 */

#ifdef __cplusplus
}
#endif

#endif
#endif

/* Log gamma function
 * \log{\Gamma(z)}
 * AS245, 2nd algorithm, http://lib.stat.cmu.edu/apstat/245
 */

static double kf_lgamma(double z)
{
	double x = 0;
	x += 0.1659470187408462e-06 / (z+7);
	x += 0.9934937113930748e-05 / (z+6);
	x -= 0.1385710331296526     / (z+5);
	x += 12.50734324009056      / (z+4);
	x -= 176.6150291498386      / (z+3);
	x += 771.3234287757674      / (z+2);
	x -= 1259.139216722289      / (z+1);
	x += 676.5203681218835      / z;
	x += 0.9999999999995183;
	return log(x) - 5.58106146679532777 - z + (z-0.5) * log(z+6.5);
}

/* complementary error function
 * \frac{2}{\sqrt{\pi}} \int_x^{\infty} e^{-t^2} dt
 * AS66, 2nd algorithm, http://lib.stat.cmu.edu/apstat/66
 */
double kf_erfc(double x)
{
	const double p0 = 220.2068679123761;
	const double p1 = 221.2135961699311;
	const double p2 = 112.0792914978709;
	const double p3 = 33.912866078383;
	const double p4 = 6.37396220353165;
	const double p5 = .7003830644436881;
	const double p6 = .03526249659989109;
	const double q0 = 440.4137358247522;
	const double q1 = 793.8265125199484;
	const double q2 = 637.3336333788311;
	const double q3 = 296.5642487796737;
	const double q4 = 86.78073220294608;
	const double q5 = 16.06417757920695;
	const double q6 = 1.755667163182642;
	const double q7 = .08838834764831844;
	double expntl, z, p;
	z = fabs(x) * M_SQRT2;
	if (z > 37.) return x > 0.? 0. : 2.;
	expntl = exp(z * z * - .5);
	if (z < 10. / M_SQRT2) // for small z
	    p = expntl * ((((((p6 * z + p5) * z + p4) * z + p3) * z + p2) * z + p1) * z + p0)
			/ (((((((q7 * z + q6) * z + q5) * z + q4) * z + q3) * z + q2) * z + q1) * z + q0);
	else p = expntl / 2.506628274631001 / (z + 1. / (z + 2. / (z + 3. / (z + 4. / (z + .65)))));
	return x > 0.? 2. * p : 2. * (1. - p);
}

/* The following computes regularized incomplete gamma functions.
 * Formulas are taken from Wiki, with additional input from Numerical
 * Recipes in C (for modified Lentz's algorithm) and AS245
 * (http://lib.stat.cmu.edu/apstat/245).
 *
 * A good online calculator is available at:
 *
 *   http://www.danielsoper.com/statcalc/calc23.aspx
 *
 * It calculates upper incomplete gamma function, which equals
 * kf_gammaq(s,z)*tgamma(s).
 */

#define KF_GAMMA_EPS 1e-14
#define KF_TINY 1e-290

// regularized lower incomplete gamma function, by series expansion
static double _kf_gammap(double s, double z)
{
	double sum, x;
	int k;
	for (k = 1, sum = x = 1.; k < 100; ++k) {
		sum += (x *= z / (s + k));
		if (x / sum < KF_GAMMA_EPS) break;
	}
	return exp(s * log(z) - z - kf_lgamma(s + 1.) + log(sum));
}
// regularized upper incomplete gamma function, by continued fraction
static double _kf_gammaq(double s, double z)
{
	int j;
	double C, D, f;
	f = 1. + z - s; C = f; D = 0.;
	// Modified Lentz's algorithm for computing continued fraction
	// See Numerical Recipes in C, 2nd edition, section 5.2
	for (j = 1; j < 100; ++j) {
		double a = j * (s - j), b = (j<<1) + 1 + z - s, d;
		D = b + a * D;
		if (D < KF_TINY) D = KF_TINY;
		C = b + a / C;
		if (C < KF_TINY) C = KF_TINY;
		D = 1. / D;
		d = C * D;
		f *= d;
		if (fabs(d - 1.) < KF_GAMMA_EPS) break;
	}
	return exp(s * log(z) - z - kf_lgamma(s) - log(f));
}

double kf_gammap(double s, double z)
{
	return z <= 1. || z < s? _kf_gammap(s, z) : 1. - _kf_gammaq(s, z);
}

double kf_gammaq(double s, double z)
{
	return z <= 1. || z < s? 1. - _kf_gammap(s, z) : _kf_gammaq(s, z);
}

/* Regularized incomplete beta function. The method is taken from
 * Numerical Recipe in C, 2nd edition, section 6.4. The following web
 * page calculates the incomplete beta function, which equals
 * kf_betai(a,b,x) * gamma(a) * gamma(b) / gamma(a+b):
 *
 *   http://www.danielsoper.com/statcalc/calc36.aspx
 */
static double kf_betai_aux(double a, double b, double x)
{
	double C, D, f;
	int j;
	if (x == 0.) return 0.;
	if (x == 1.) return 1.;
	f = 1.; C = f; D = 0.;
	// Modified Lentz's algorithm for computing continued fraction
	for (j = 1; j < 200; ++j) {
		double aa, d;
		int m = j>>1;
		aa = (j&1)? -(a + m) * (a + b + m) * x / ((a + 2*m) * (a + 2*m + 1))
			: m * (b - m) * x / ((a + 2*m - 1) * (a + 2*m));
		D = 1. + aa * D;
		if (D < KF_TINY) D = KF_TINY;
		C = 1. + aa / C;
		if (C < KF_TINY) C = KF_TINY;
		D = 1. / D;
		d = C * D;
		f *= d;
		if (fabs(d - 1.) < KF_GAMMA_EPS) break;
	}
	return exp(kf_lgamma(a+b) - kf_lgamma(a) - kf_lgamma(b) + a * log(x) + b * log(1.-x)) / a / f;
}

double kf_betai(double a, double b, double x)
{
	return x < (a + 1.) / (a + b + 2.)? kf_betai_aux(a, b, x) : 1. - kf_betai_aux(b, a, 1. - x);
}

#ifdef KF_MAIN
#include <stdio.h>
int main(int argc, char *argv[])
{
	double x = 5.5, y = 3;
	double a, b;
	printf("erfc(%lg): %lg, %lg\n", x, erfc(x), kf_erfc(x));
	printf("upper-gamma(%lg,%lg): %lg\n", x, y, kf_gammaq(y, x)*tgamma(y));
	a = 2; b = 2; x = 0.5;
	printf("incomplete-beta(%lg,%lg,%lg): %lg\n", a, b, x, kf_betai(a, b, x) / exp(kf_lgamma(a+b) - kf_lgamma(a) - kf_lgamma(b)));
	return 0;
}
#endif


// log\binom{n}{k}
static double lbinom(int n, int k)
{
    if (k == 0 || n == k) return 0;
    return lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1);
}

// n11  n12  | n1_
// n21  n22  | n2_
//-----------+----
// n_1  n_2  | n

// hypergeometric distribution
static double hypergeo(int n11, int n1_, int n_1, int n)
{
    return exp(lbinom(n1_, n11) + lbinom(n-n1_, n_1-n11) - lbinom(n, n_1));
}

typedef struct {
    int n11, n1_, n_1, n;
    double p;
} hgacc_t;

// incremental version of hypergenometric distribution
static double hypergeo_acc(int n11, int n1_, int n_1, int n, hgacc_t *aux)
{
    if (n1_ || n_1 || n) {
        aux->n11 = n11; aux->n1_ = n1_; aux->n_1 = n_1; aux->n = n;
    } else { // then only n11 changed; the rest fixed
        if (n11%11 && n11 + aux->n - aux->n1_ - aux->n_1) {
            if (n11 == aux->n11 + 1) { // incremental
                aux->p *= (double)(aux->n1_ - aux->n11) / n11
                    * (aux->n_1 - aux->n11) / (n11 + aux->n - aux->n1_ - aux->n_1);
                aux->n11 = n11;
                return aux->p;
            }
            if (n11 == aux->n11 - 1) { // incremental
                aux->p *= (double)aux->n11 / (aux->n1_ - n11)
                    * (aux->n11 + aux->n - aux->n1_ - aux->n_1) / (aux->n_1 - n11);
                aux->n11 = n11;
                return aux->p;
            }
        }
        aux->n11 = n11;
    }
    aux->p = hypergeo(aux->n11, aux->n1_, aux->n_1, aux->n);
    return aux->p;
}

double kt_fisher_exact(int n11, int n12, int n21, int n22, double *_left, double *_right, double *two)
{
    int i, j, max, min;
    double p, q, left, right;
    hgacc_t aux;
    int n1_, n_1, n;

    n1_ = n11 + n12; n_1 = n11 + n21; n = n11 + n12 + n21 + n22; // calculate n1_, n_1 and n
    max = (n_1 < n1_) ? n_1 : n1_; // max n11, for right tail
    min = n1_ + n_1 - n;    // not sure why n11-n22 is used instead of min(n_1,n1_)
    if (min < 0) min = 0; // min n11, for left tail
    *two = *_left = *_right = 1.;
    if (min == max) return 1.; // no need to do test
    q = hypergeo_acc(n11, n1_, n_1, n, &aux); // the probability of the current table
    // left tail
    p = hypergeo_acc(min, 0, 0, 0, &aux);
    for (left = 0., i = min + 1; p < 0.99999999 * q && i<=max; ++i) // loop until underflow
        left += p, p = hypergeo_acc(i, 0, 0, 0, &aux);
    --i;
    if (p < 1.00000001 * q) left += p;
    else --i;
    // right tail
    p = hypergeo_acc(max, 0, 0, 0, &aux);
    for (right = 0., j = max - 1; p < 0.99999999 * q && j>=0; --j) // loop until underflow
        right += p, p = hypergeo_acc(j, 0, 0, 0, &aux);
    ++j;
    if (p < 1.00000001 * q) right += p;
    else ++j;
    // two-tail
    *two = left + right;
    if (*two > 1.) *two = 1.;
    // adjust left and right
    if (abs(i - n11) < abs(j - n11)) right = 1. - left + q;
    else left = 1.0 - right + q;
    if (left > 1.0) left = 1.0;
    if (right > 1.0) right = 1.0;
    *_left = left; *_right = right;
    return q;
}


#define SIZE 50
#define each_i(start, end) for (i = start; i < end; ++i)

typedef enum { UP, DOWN } direction;

typedef struct { int index; double value; } iv1;

typedef struct { int index; int value; } iv2;


int compare_iv1(const void *a, const void *b) {
    double aa = ((iv1 *)a) -> value;
    double bb = ((iv1 *)b) -> value;
    if (aa > bb) return 1;
    if (aa < bb) return -1;
    return 0;
}

int compare_iv1_desc(const void *a, const void *b) {
    return -compare_iv1(a, b);
}

int compare_iv2(const void *a, const void *b) {
    return ((iv2 *)a) -> value - ((iv2 *)b) -> value;
}

void ratchet(double *pa, int size, direction dir) {
    int i;
    double m = pa[0];
    if (dir == UP) {
        each_i(1, size) {
            if (pa[i] > m) pa[i] = m;
            m = pa[i];
        }
    }
    else {
        each_i(1, size) {
            if (pa[i] < m) pa[i] = m;
            m = pa[i];
        }
    }
    each_i(0, size) if (pa[i] > 1.0) pa[i] = 1.0;
}


static void schwartzian(const double *p, double *pa, int size, direction dir) 
{
    int i;
#if 0
    int order[SIZE];
    int order2[SIZE];
    iv1 iv1s[SIZE];
    iv2 iv2s[SIZE];
    double pa2[SIZE];
#endif
    int *order;
    int *order2;
    iv1 *iv1s;
    iv2 *iv2s;
    double *pa2;

    order = malloc(sizeof(int) * size); 
    order2 = malloc(sizeof(int) * size); 
    iv1s = malloc(sizeof(iv1) * size); 
    iv2s = malloc(sizeof(iv2) * size); 
    pa2 = malloc(sizeof(double) * size); 

    each_i(0, SIZE) { iv1s[i].index = i; iv1s[i].value = p[i]; }
    if (dir == UP)
        qsort(iv1s, SIZE, sizeof(iv1s[0]), compare_iv1_desc);
    else
        qsort(iv1s, SIZE, sizeof(iv1s[0]), compare_iv1);
    each_i(0, SIZE) order[i] = iv1s[i].index;
    each_i(0, SIZE) pa[i] *= p[order[i]];
    ratchet(pa, size, dir);
    each_i(0, SIZE) { iv2s[i].index = i; iv2s[i].value = order[i]; }
    qsort(iv2s, SIZE, sizeof(iv2s[0]), compare_iv2);
    each_i(0, SIZE) order2[i] = iv2s[i].index;
    each_i(0, SIZE) pa2[i] = pa[order2[i]];
    each_i(0, SIZE) pa[i] = pa2[i];

    if (order) free ( order ) ;
    if (order2) free( order2 ); 
    if (iv1s) free( iv1s ); 
    if (iv2s) free( iv2s ); 
    if (pa2) free( pa2 ); 
}

void adjust(const double *p, double *pa, int size) 
{
    int i;
// always BH
     each_i(0, size) pa[i] = (double)size / (size - i);
     schwartzian(p, pa, size, UP);

#if 0
    if (!strcmp(type, "Benjamini-Hochberg")) 
    {
     each_i(0, size) pa[i] = (double)size / (size - i);
     schwartzian(p, pa, UP);
    }
    else if (!strcmp(type, "Benjamini-Yekutieli")) {
        double q = 0.0;
        each_i(1, SIZE + 1) q += 1.0 / i;
        each_i(0, SIZE) pa[i] = q * SIZE / (SIZE - i);
        schwartzian(p, pa, UP);
    }
    else if (!strcmp(type, "Bonferroni")) {
        each_i(0, SIZE) pa[i] = (p[i] * SIZE > 1.0) ? 1.0 : p[i] * SIZE;
    }
    else if (!strcmp(type, "Hochberg")) {
        each_i(0, SIZE) pa[i]  = i + 1.0;
        schwartzian(p, pa, UP);
    }
    else if (!strcmp(type, "Holm")) {
        each_i(0, SIZE) pa[i] = SIZE - i;
        schwartzian(p, pa, DOWN);
    }
    else if (!strcmp(type, "Hommel")) {
        int i, j;
        int order[SIZE];
        int order2[SIZE];
        iv1 iv1s[SIZE];
        iv2 iv2s[SIZE];
        double s[SIZE];
        double q[SIZE];
        double pa2[SIZE];
        int indices[SIZE];
        each_i(0, SIZE) { iv1s[i].index = i; iv1s[i].value = p[i]; }
        qsort(iv1s, SIZE, sizeof(iv1s[0]), compare_iv1);
        each_i(0, SIZE) order[i] = iv1s[i].index;
        each_i(0, SIZE) s[i] = p[order[i]];
        double min = s[0] * SIZE;
        each_i(1, SIZE) {
            double temp = s[i] / (i + 1.0);
            if (temp < min) min = temp;
        }
        each_i(0, SIZE) q[i] = min;
        each_i(0, SIZE) pa2[i] = min;
        for (j = SIZE - 1; j >= 2; --j) {
            each_i(0, SIZE) indices[i] = i;
            int upper_start = SIZE - j + 1;      /* upper indices start index */
            int upper_size = j - 1;              /* size of upper indices */
            int lower_size = SIZE - upper_size;  /* size of lower indices */
            double qmin = j * s[indices[upper_start]] / 2.0;
            each_i(1, upper_size) {
                double temp = s[indices[upper_start + i]] * j / (2.0 + i);
                if (temp < qmin) qmin = temp;
            }
            each_i(0, lower_size) {
                double temp = s[indices[i]] * j;
                q[indices[i]] = (temp < qmin) ? temp : qmin;
            }
            each_i(0, upper_size) q[indices[upper_start + i]] = q[SIZE - j];
            each_i(0, SIZE) if (pa2[i] < q[i]) pa2[i] = q[i];
        }
        each_i(0, SIZE) { iv2s[i].index = i; iv2s[i].value = order[i]; }
        qsort(iv2s, SIZE, sizeof(iv2s[0]), compare_iv2);
        each_i(0, SIZE) order2[i] = iv2s[i].index;
        each_i(0, SIZE) pa[i] = pa2[order2[i]];
    }
    else if (!strcmp(type, "Šidák")) {
        each_i(0, SIZE) pa[i] = 1.0 - pow(1.0 - p[i], SIZE);
    }
    else {
        printf("\nSorry, do not know how to do '%s' correction.\n", type);
        printf("Perhaps you want one of these?:\n");
        each_i(0, 7) printf("  %s\n", types[i]);
        exit(1);
    }
#endif
}

void bh_adjusted(const double *p, double *pa, int size) 
{
    adjust(p, pa,size);
#if 0
    if (check(p)) 
    {
    adjust(p, pa,size);
    }
    else 
    {
        printf("p-values must be in range 0.0 to 1.0\n");
        exit(1);
    }
#endif
}

#if 0
int check(const double* p,int size) 
{
    int i;

    for (i=0;i<size;i++)
    {
        if (p[i] < 0.0 || p[i] > 1.0) return 0;
    }
    return 1;
}

double p_values[50] = {
        4.533744e-01, 7.296024e-01, 9.936026e-02, 9.079658e-02, 1.801962e-01,
        8.752257e-01, 2.922222e-01, 9.115421e-01, 4.355806e-01, 5.324867e-01,
        4.926798e-01, 5.802978e-01, 3.485442e-01, 7.883130e-01, 2.729308e-01,
        8.502518e-01, 4.268138e-01, 6.442008e-01, 3.030266e-01, 5.001555e-02,
        3.194810e-01, 7.892933e-01, 9.991834e-01, 1.745691e-01, 9.037516e-01,
        1.198578e-01, 3.966083e-01, 1.403837e-02, 7.328671e-01, 6.793476e-02,
        4.040730e-03, 3.033349e-04, 1.125147e-02, 2.375072e-02, 5.818542e-04,
        3.075482e-04, 8.251272e-03, 1.356534e-03, 1.360696e-02, 3.764588e-04,
        1.801145e-05, 2.504456e-07, 3.310253e-02, 9.427839e-03, 8.791153e-04,
        2.177831e-04, 9.693054e-04, 6.610250e-05, 2.900813e-02, 5.735490e-03
};

   ** TESTING CODE
int main_bh() 
{
    int i;
    double *pa;
    int size;


    size = 50;
    pa = malloc(sizeof(double)*size); //  double pa[SIZE] = { 0.0 };
    for (i=0;i<size;i++)
    {
        pa[i] = (double)size / (size - i);
    }
    schwartzian(p_values, pa, size, UP);
    for (i=0;i<size;i++)
    {
        if (!(i % 5)) printf("\n[%2d]  ", i);
            printf("%1.10f ", pa[i]);
    }
    printf("\n");
    if (pa) free(pa);
    return 0;
}


    // testing routine
int main_fe() 
{
    char s[1002];
    double fisher_left_p, fisher_right_p, fisher_twosided_p;
    FILE *fp;
    FILE *fpo;
    int a,b,c,d;

    fpo = fopen("fetest.R","w");
    fp = fopen("fe_testcases","r");
    fprintf(fpo,"sink(\"outr\")\n");
    while (fgets(s,1000,fp))
    {
        sscanf(s,"%d %d %d %d",&a,&b,&c,&d);     
        kt_fisher_exact(a, b, c, d,  &fisher_left_p, &fisher_right_p,&fisher_twosided_p);
        printf("%18.16f %18.16f %18.16f %d %d %d %d\n",fisher_twosided_p,fisher_right_p,fisher_left_p,a,b,c,d);
fprintf(fpo,"a=%d\n",a);
fprintf(fpo,"b=%d\n",b);
fprintf(fpo,"c=%d\n",c);
fprintf(fpo,"d=%d\n",d);
fprintf(fpo,"x=c(a,b,c,d)\n");
fprintf(fpo,"data  = matrix(x,nrow=2)\n");
fprintf(fpo,"y1=fisher.test(data)\n");
fprintf(fpo,"p1=y1$p.value\n");
fprintf(fpo,"y2=fisher.test(data,alternative=\"greater\")\n");
fprintf(fpo,"p2=y2$p.value\n");
fprintf(fpo,"y3=fisher.test(data,alternative=\"less\")\n");
fprintf(fpo,"p3=y3$p.value\n");
fprintf(fpo,"cat(sprintf(\"%%18.16f %%18.16f %%18.16f %%d %%d %%d %%d\\n\",p1,p2,p3,a,b,c,d))\n");
    }
    fclose(fp);
    fprintf(fpo,"sink()\n");
    fclose(fpo);
    return 0;
}

int main()
{
    main_fe() ;
    main_bh() ;
}
//x=c( 9,35,1274,36863); data  = matrix(x,nrow=2); fisher.test(data)
#endif
