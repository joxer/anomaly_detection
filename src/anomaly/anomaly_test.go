package anomaly

import (
	"testing"
	"math/rand"
	"github.com/joxer/golang-api/plotly"
)
func TestAnomaly(t *testing.T) {

	series := make(map[string][]float64)
	rand.Seed(100)
	series["a"] = make([]float64,200)


	for i := 0; i < 200; i++ {
		series["a"][i] = rand.NormFloat64()+5

	}

	features := []string{"a"}

	n := New().setSeries(series).setFeatures(features).build()
	n.train()

	x := make([]interface{}, 0)
	y := make([]interface{},0)
	for i := 0; i < 200; i++ {
		x = append(x, ([]float64{series["a"][i]}))
		y = append(y,i)
	}

	f := plotly.Figure{
		Data: []plotly.Trace{
			plotly.Trace{
				Type: "bar",
				X: x,
				Y: y,
			},
		},
	}


	f.Save("new golang file")
}
