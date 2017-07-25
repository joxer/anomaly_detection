package anomaly
/*
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_sort.h>
#include <math.h>
#include <stdio.h>

#ifndef min
# define min(a, b)      ((a) > (b) ? (b) : (a))
#endif

#ifndef sign
# define sign(a)        ((a) >= 0 ? 1.0 : -1.0)
#endif


double poly(const double *cc, int nord, double x)
{

int j;
double p, ret_val;

ret_val = cc[0];
if (nord > 1) {
p = x * cc[nord-1];
for (j = nord - 2; j > 0; j--)
p = (p + cc[j]) * x;

ret_val += p;
}
return ret_val;
}

void swilk(int *init,
double *x, int *n, int *n1, int *n2,
double *a,
double *w, double *pw, int *ifault)
{

const static double zero = 0.f;
const static double one = 1.f;
const static double two = 2.f;

const static double small = 1e-19f;

const static double g[2] = { -2.273f,.459f };
const static double
c1[6] = { 0.f,.221157f,-.147981f,-2.07119f, 4.434685f, -2.706056f },
c2[6] = { 0.f,.042981f,-.293762f,-1.752461f,5.682633f, -3.582633f };
const static double c3[4] = { .544f,-.39978f,.025054f,-6.714e-4f };
const static double c4[4] = { 1.3822f,-.77857f,.062767f,-.0020322f };
const static double c5[4] = { -1.5861f,-.31082f,-.083751f,.0038915f };
const static double c6[3] = { -.4803f,-.082676f,.0030302f };
const static double c7[2] = { .164f,.533f };
const static double c8[2] = { .1736f,.315f };
const static double c9[2] = { .256f,-.00635f };

double r__1;

int i, j, ncens, i1, nn2;

double zbar, ssassx, summ2, ssumm2, gamma, delta, range;
double a1, a2, an, bf, ld, m, s, sa, xi, sx, xx, y, w1;
double fac, asa, an25, ssa, z90f, sax, zfm, z95f, zsd, z99f, rsn, ssx, xsx;

*pw = 1.;
if (*w >= 0.) {
*w = 1.;
}
if (*n < 3) {   *ifault = 1; return;
}

an = (double) (*n);
nn2 = *n / 2;
if (*n2 < nn2) {    *ifault = 3; return;
}
if (*n1 < 3) {  *ifault = 1; return;
}
ncens = *n - *n1;
if (ncens < 0 || (ncens > 0 && *n < 20)) {  *ifault = 4; return;
}
if (ncens > 0) {
delta = (double) ncens / an;
if (delta > .8f) {  *ifault = 5; return;
}
} else { delta = 0.f; }

--a;
if (! (*init)) {
if (*n == 3) {
const static double sqrth = .70710678f;
a[1] = sqrth;
} else {
an25 = an + .25;
summ2 = zero;
for (i = 1; i <= *n2; ++i) {
a[i] = gsl_cdf_ugaussian_Pinv ((i - .375f) / an25);
r__1 = a[i];
summ2 += r__1 * r__1;
}
summ2 *= two;
ssumm2 = sqrt(summ2);
rsn = one / sqrt(an);
a1 = poly(c1, 6, rsn) - a[1] / ssumm2;

if (*n > 5) {
i1 = 3;
a2 = -a[2] / ssumm2 + poly(c2, 6, rsn);
fac = sqrt((summ2 - two * (a[1] * a[1]) - two * (a[2] * a[2]))
/ (one - two * (a1 * a1) - two * (a2 * a2)));
a[2] = a2;
} else {
i1 = 2;
fac = sqrt((summ2 - two * (a[1] * a[1])) /
( one  - two * (a1 * a1)));
}
a[1] = a1;
for (i = i1; i <= nn2; ++i)
a[i] /= - fac;
}
*init = (1);
}

if (*w < zero) {
w1 = 1. + *w;
*ifault = 0;
goto L70;
}

range = x[*n1 - 1] - x[0];
if (range < small) {
*ifault = 6;    return;
}

*ifault = 0;
xx = x[0] / range;
sx = xx;
sa = -a[1];
j = *n - 1;
for (i = 1; i < *n1; --j) {
xi = x[i] / range;
if (xx - xi > small) {
*ifault = 7;
}
sx += xi;
++i;
if (i != j)
sa += sign(i - j) * a[min(i,j)];
xx = xi;
}
if (*n > 5000) {
*ifault = 2;
}


sa /= *n1;
sx /= *n1;
ssa = ssx = sax = zero;
j = *n - 1;
for (i = 0; i < *n1; ++i, --j) {
if (i != j)
asa = sign(i - j) * a[1+min(i,j)] - sa;
else
asa = -sa;
xsx = x[i] / range - sx;
ssa += asa * asa;
ssx += xsx * xsx;
sax += asa * xsx;
}

ssassx = sqrt(ssa * ssx);
w1 = (ssassx - sax) * (ssassx + sax) / (ssa * ssx);
L70:
*w = 1. - w1;

if (*n == 3) {
const static double pi6 = 1.90985931710274;
const static double stqr= 1.04719755119660;
*pw = pi6 * (asin(sqrt(*w)) - stqr);
if(*pw < 0.) *pw = 0.;
return;
}
y = log(w1);
xx = log(an);
if (*n <= 11) {
gamma = poly(g, 2, an);
if (y >= gamma) {
*pw = 1e-99;
return;
}
y = -log(gamma - y);
m = poly(c3, 4, an);
s = exp(poly(c4, 4, an));
} else {
m = poly(c5, 4, xx);
s = exp(poly(c6, 3, xx));
}

if (ncens > 0) {
const static double three = 3.f;

const static double z90 = 1.2816f;
const static double z95 = 1.6449f;
const static double z99 = 2.3263f;
const static double zm = 1.7509f;
const static double zss = .56268f;
const static double bf1 = .8378f;

const static double xx90 = .556;
const static double xx95 = .622;

ld = -log(delta);
bf = one + xx * bf1;
r__1 = pow(xx90, (double) xx);
z90f = z90 + bf * pow(poly(c7, 2, r__1), (double) ld);
r__1 = pow(xx95, (double) xx);
z95f = z95 + bf * pow(poly(c8, 2, r__1), (double) ld);
z99f = z99 + bf * pow(poly(c9, 2, xx), (double)ld);


zfm = (z90f + z95f + z99f) / three;
zsd = (z90 * (z90f - zfm) +
z95 * (z95f - zfm) + z99 * (z99f - zfm)) / zss;
zbar = zfm - zsd * zm;
m += zbar * s;
s *= zsd;
}

*pw = gsl_cdf_gaussian_Q(y - m, s);


return;
}


shapiro_wilk(double d_data[], int d_n,double arg[]){
int init = 0;
int n = d_n;
int n1 = d_n;
int n2 = d_n/2;
int error = 0;
double d_w = 0.0;
double d_pValue = 0.0;
double a[n2];


swilk(&init, d_data, &n, &n1, &n2, a, &d_w, &d_pValue, &error);
arg[0] = d_w;
arg[1] = d_pValue;
arg[2] = error;
}
 */
// #cgo LDFLAGS: -lgsl  -lgslcblas -lm
import "C"
import (
	"github.com/montanaflynn/stats"
	"math"
	"unsafe"

	"fmt"
	"sort"
)



type Anomaly struct {
	series            map[string][]float64
	features          []string
	internal_mean_std [][]float64
}

type anomalyBuilder struct{
	series map[string][]float64
	features []string
}

func (a *anomalyBuilder)setSeries(series map[string][] float64) *anomalyBuilder{
	a.series = series
	return a
}

func (a *anomalyBuilder)setFeatures(features []string) *anomalyBuilder {
	a.features = features
	return a
}

func (a *anomalyBuilder)build() *Anomaly {
	ret := Anomaly{a.series,a.features,nil}
	ret.internal_mean_std = make([][]float64,len(ret.features))
	return &ret
}


func New() *anomalyBuilder {
	return &anomalyBuilder{}
}


func (a *Anomaly) train(){

	for i, idx := range(a.features) {
		mean,_ := stats.Mean(a.series[idx])
		std,_ := stats.StandardDeviation(a.series[idx])
		a.internal_mean_std[i] = []float64{mean, std}
	}


}

func (a *Anomaly) fit(fit_value []float64) float64{
	values := 1.0

	for i := range(a.internal_mean_std){
		value := -math.Pow((fit_value[i]-a.internal_mean_std[i][0]),2)/(2*math.Pow(a.internal_mean_std[i][1],2))
		iret := (1/(math.Sqrt(2*math.Pi*math.Pow(a.internal_mean_std[i][1],2))))*math.Exp(value)
		values = values*iret
	}
	return values
}

func (anomaly *Anomaly) ShapiroWilk(series []float64)  {
	sort.Float64s(series)
	var tmp = make([]C.double,len(series))
	for i := 0 ; i < len(series);i++{
		tmp[i] = C.double(series[i])
	}
	tmp2 := make([]C.double,3)
	ptr := unsafe.Pointer(&tmp[0])
	ptr2 := unsafe.Pointer(&tmp2[0])
	C.shapiro_wilk((*C.double)(ptr),C.int(len(series)),(*C.double)(ptr2))

	fmt.Print(tmp2)


}