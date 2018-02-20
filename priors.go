package cophymaru

import (
	"fmt"
	"math"
	"os"
)

//BranchLengthPrior is a struct for specifying and calculating different priors
type BranchLengthPrior struct {
	TYPE   string
	ALPHA  float64
	BETA   float64
	NTIPS  int
	SGAMMA float64
	SFACT  float64
}

//Calc will return the log prior probability of the branch length prior
func (pr *BranchLengthPrior) Calc(nodes []*Node) float64 {
	if pr.TYPE == "1" {
		return ExponentialBranchLengthLogPrior(nodes, pr.BETA)
	} else if pr.TYPE == "2" {
		return DirichletBranchLengthLogPrior(nodes, pr)
	}
	return 0.0
}

//InitializePrior will set up a new instance of a prior
func InitializePrior(priorType string, nodes []*Node) *BranchLengthPrior {
	pr := new(BranchLengthPrior)
	pr.TYPE = priorType
	if pr.TYPE == "0" {
		pr.BETA = 0.
	} else if pr.TYPE == "1" {
		pr.BETA = 1.0
	} else if pr.TYPE == "2" {
		pr.ALPHA = 1.0
		pr.BETA = 10.0
		ntips := 0
		for _, n := range nodes {
			if len(n.CHLD) == 0 {
				ntips++
			}
		}
		pr.NTIPS = ntips
		pr.SFACT = factorial(ntips)
		pr.SGAMMA = math.Gamma(pr.ALPHA)

	} else {
		fmt.Println("PLEASE SPECIFY VALID OPTION FOR PRIOR. TYPE maru -h TO SEE OPTIONS")
		os.Exit(0)
	}
	return pr
}

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

//DirichletBranchLengthLogPrior will return the log prior probability of all branch lengths in a trejke
func DirichletBranchLengthLogPrior(nodes []*Node, pr *BranchLengthPrior) float64 {
	T := TreeLength(nodes)
	x := (2. * float64(pr.NTIPS)) - 4.
	prob := math.Log(math.Pow(pr.BETA, pr.ALPHA)/pr.SGAMMA) + (-pr.BETA * T) + math.Log(math.Pow(T, (pr.ALPHA-1.-x)))
	prob = prob + math.Log(pr.SFACT)
	return prob
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
