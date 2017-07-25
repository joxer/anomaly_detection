package anomaly

import (
	"testing"
	"math/rand"
)
func TestAnomaly(t *testing.T) {

	series :=  make(map[string][]float64)
	rand.Seed(100)
	series["a"] = []float64{rand.NormFloat64()*2+3,
		rand.NormFloat64()*2+3,
		rand.NormFloat64()*2+3,
		rand.NormFloat64()*2+3,rand.NormFloat64()*2+3,
	}




	features := []string{"a"}

	n := New().setSeries(series).setFeatures(features).build()
	n.train()
	n.ShapiroWilk(series["a"])
	//f.Save("new golang file")
}
