package anomaly

import (
  "github.com/montanaflynn/stats"
	"fmt"
	"math"
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
	fmt.Print(a.internal_mean_std)
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

func (anomaly *Anomaly) ShapiroWilk(series []float64) {

	N := float64(len(series))
	n := int(N)
	mean,_ := stats.Mean(series)
	Y := series

	sort.Float64s(Y)

	m := make([]float64,N)

	for i := 0 ; i < int(N);i++{
		num := float64((i+1) - 3.0/8.0)
		den := float64(N+1.0/4.0)
		m[i] = InvCumulativeNormalDistribution(num/den)
	}
	mm := dot(m,m)
	c := divide(m, math.Sqrt(N))
	u := 1.0/math.Sqrt(N)

	u2 := u*u
	u3 := u2*u
	u4 := u2*u2
	u5 := u4*u
	a := make([]float64,int(N))

	a[n-1] = -2.706056 * u5 + 4.434685 * u4 - 2.071190 * u3 - 0.147981 * u2 + 0.221157 * u + c[n - 1]
	an := a[n-1]
	mn := m[n-1]
	a[0] = -an

	if (n <= 5) {

		phi := (mm - 2*mn*mn) / (1.0 - 2*an*an)
		sqrt := math.Sqrt(phi)

		for i := 1; i < n-1; i++ {
			a[i] = m[i] / sqrt
		}
	}else{
		a[n - 2] = -3.582633 * u5 + 5.682633 * u4 - 1.752461 * u3 - 0.293762 * u2 + 0.042981 * u + c[n - 2];
		anm1 :=a[n-2]
		mnm1 :=m[n-2]
		a[1] = -anm1
		phi := (mm - 2 * mn * mn - 2 * mnm1 * mnm1) / (1.0 - 2 * an * an - 2 * anm1 * anm1)
		sqrt := math.sqrt(phi)
		for i:= 2 ;i < n-2;i++ {
			a[i] = m[i]/sqrt
		}
	}

	Wnum := 0.0
	Wden := 0.0

	for i:=0; i < len(series);i++ {
		w := (Y[i] -mean)
		Wnum += a[i] * Y[i]
		Wden += w*w
	}

	W := Wnum*Wnum / Wden
	return W,
}

func divide(s []float64,n float64) []float64{
	for i:=0; i < len(s);i++{
		s[i] /= n
	}
	return s
}

func InvCumulativeNormalDistribution(p float64) float64 {
	if p <= 0 || p >= 1 {
		panic("Argument to ltqnorm %f must be in open interval (0,1)")
	}

	//Coefficients in rational approximations.
	a := []float64{-3.969683028665376e+01, 2.209460984245205e+02,
				   -2.759285104469687e+02, 1.383577518672690e+02,
				   -3.066479806614716e+01, 2.506628277459239e+00}
	b := []float64{-5.447609879822406e+01, 1.615858368580409e+02,
				   -1.556989798598866e+02, 6.680131188771972e+01,
				   -1.328068155288572e+01}
	c := []float64{-7.784894002430293e-03, -3.223964580411365e-01,
				   -2.400758277161838e+00, -2.549732539343734e+00,
				   4.374664141464968e+00, 2.938163982698783e+00}
	d := []float64{7.784695709041462e-03, 3.224671290700398e-01,
				   2.445134137142996e+00, 3.754408661907416e+00}

	// Define break-points.
	plow := 0.02425
	phigh := 1 - plow

	// Rational approximation for lower region:
	if p < plow {
		q := math.Sqrt(-2 * math.Log(p))
		return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q + c[5]) / ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q + 1)
	}

	// Rational approximation for upper region:
	if phigh < p {
		q := math.Sqrt(-2 * math.Log(1-p))
		return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q + c[5]) / ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q + 1)
	}

	//Rational approximation for central region:
	q := p - 0.5
	r := q * q
	return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r + a[5]) * q / (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r + 1)
}


func dot(s1 []float64, s2[]float64) float64{
	var tot float64;
	for i:=0; i < len(s1);i++ {
		tot += s1[i]*s2[i]
	}
	return tot
}