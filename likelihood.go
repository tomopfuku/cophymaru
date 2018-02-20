package cophymaru

//LL is a struct for likelihood functions/calculations
type LL struct {
	MULTI bool
}

//Calc will calculate the log-likelihood
func (ll *LL) Calc(tree *Node, startFresh bool, weights []float64) float64 {
	if ll.MULTI == true {
		return WeightedUnrootedLogLikeParallel(tree, startFresh, weights)
	}
	return WeightedUnrootedLogLike(tree, startFresh, weights)
}

//InitLL will initialize the likelihood struct
func InitLL(multithread bool) *LL {
	ll := new(LL)
	if multithread == true {
		ll.MULTI = true
	} else {
		ll.MULTI = false
	}
	return ll
}
