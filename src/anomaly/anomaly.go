package anomaly

import (
  "github.com/montanaflynn/stats"
	"fmt"
	"math"
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
	fmt.Print(values)
	return values
}