package anomaly

import (
	"testing"
	"math/rand"
)
func TestAnomaly(t *testing.T) {

	series := make(map[string][]float64)
	rand.Seed(10000)
	series["a"] = make([]float64,10000)

	for i := 0; i < 10000; i++ {
		series["a"][i] = rand.NormFloat64()*1*1

	}

	features := []string{"a"}

	n := New().setSeries(series).setFeatures(features).build()
	n.train()
	n.fit([]float64{series["a"][2]})
}
