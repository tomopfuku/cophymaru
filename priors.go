package cophymaru

import (
	"math"
)

func normalPDF(p float64, mean float64, sd float64) float64 {
	prob := (1.0 / math.Sqrt(2.0*math.Pi*sd)) * math.Exp(-(math.Pow(p-mean, 2.) / (2 * sd)))
	return prob
}

func gammaPDF(p, alpha, beta float64) float64 {
	prob := (math.Pow(alpha, beta) / math.Gamma(alpha)) * math.Exp(-beta*p) * math.Pow(p, (alpha-1))
	return prob
}

//GammaTreeLengthPrior will calculate the prior probability of the current tree length
func GammaTreeLengthPrior(nodels []*Node, alpha, beta float64) float64 {
	tl := TreeLength(nodels)
	return gammaPDF(tl, alpha, beta)
}

//DirchletBranchLengthLogPrior will return the log prior probability of all branch lengths in a trejke
func DirchletBranchLengthLogPrior(nodels []*Node, alpha, beta float64) float64 {
	T := TreeLength(nodels)
	ntips := 0
	for _, n := range nodels {
		if len(n.CHLD) == 0 {
			ntips++
		}
	}
	x := (2. * float64(ntips)) - 4.
	prob := (math.Pow(beta, alpha) / math.Gamma(alpha)) * math.Exp(-beta*T) * math.Pow(T, (alpha-1.-x))
	prob = prob * factorial(ntips)
	return math.Log(prob)
}

func factorial(val int) float64 {
	fact := float64(val)
	for i := val; i > 0; i-- {
		fact = fact * float64(i)
	}
	return float64(fact)
}

//this calculates the probability of drawing parameter p from an exponential distribution with shape parameter == lambda
func exponentialPDF(p float64, beta float64) float64 {
	prob := (1.0 / beta) * math.Exp(-(p / beta))
	return prob
}

//ExponentialBranchLengthLogPrior will give the exponential prior probability of observing a particular set of branch lengths on a tree, given shape lambda
func ExponentialBranchLengthLogPrior(nodels []*Node, beta float64) (prob float64) {
	var tempp float64
	for _, n := range nodels {
		tempp = exponentialPDF(n.LEN, beta)
		prob += math.Log(tempp)
	}
	return
}

//EBExponentialBranchLengthLogPrior will give the exponential probability of each branch length when exponential mean == ML branch length
func EBExponentialBranchLengthLogPrior(nodels []*Node, means []float64) (prob float64) {
	var tempp float64
	for i, n := range nodels {
		tempp = exponentialPDF(n.LEN, means[i])
		prob += math.Log(tempp)
	}
	return
}

//EBNormalBranchLengthLogPrior will calculate the prior probability of a branch length from a normal distribution centered on the ML branch length estimate
func EBNormalBranchLengthLogPrior(nodels []*Node, means []float64) (prob float64) {
	var tempp float64
	for i, n := range nodels {
		tempp = normalPDF(n.LEN, means[i], 0.05)
		prob += math.Log(tempp)
	}
	return
}
