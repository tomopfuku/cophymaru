package cophymaru

import (
	"fmt"

	"gonum.org/v1/gonum/optimize"
)

func OptimizeMorphStratHeights(tree *Node) {
	fcn := func(heights []float64) float64 {
		large := 100000000000.0
		for _, i := range heights {
			if i <= 0.0 {
				return large
			}
		}
		preNodes := tree.PreorderArray()
		bad := AssignInternalNodeHeights(preNodes, heights)
		if bad {
			return large
		}
		morphLL := RootedLogLikeParallel(tree, true, 4)
		stratLL := PoissonTreeLoglike(preNodes)
		lnl := morphLL + stratLL
		return -lnl
	}
	settings := optimize.Settings{} //DefaultSettings()
	settings.MajorIterations = 10
	settings.Concurrent = 0
	settings.FuncEvaluations = 100
	//settings.FunctionThreshold = 0.1
	settings.GradientThreshold = 0.1
	settings.Recorder = nil
	//FC := optimize.FunctionConverge{}
	//FC.Absolute = 10
	//FC.Relative = 10
	//FC.Iterations = 10
	//settings.FunctionConverge = &FC
	p := optimize.Problem{Func: fcn, Grad: nil, Hess: nil}
	var p0 []float64
	for _, n := range tree.PreorderArray() {
		if len(n.CHLD) > 0 {
			p0 = append(p0, n.HEIGHT)
		}
	}
	meth := &optimize.NelderMead{}
	_, err := optimize.Minimize(p, p0, nil, meth)
	if err != nil {
		fmt.Println(err)
	}
	//fmt.Println(res)
}

func OptimizeGlobalRateHeights(tree *Node) {
	fcn := func(params []float64) float64 {
		rate := params[0]
		heights := params[1:]
		large := 100000000000.0
		for _, i := range heights {
			if i <= 0.0 {
				return large
			}
		}
		preNodes := tree.PreorderArray()
		AssignGlobalRate(preNodes, rate)
		bad := AssignInternalNodeHeights(preNodes, heights)
		if bad {
			return large
		}
		morphLL := RootedLogLikeParallel(tree, true, 4)
		stratLL := PoissonTreeLoglike(preNodes)
		lnl := morphLL + stratLL
		return -lnl
	}
	settings := optimize.Settings{} //DefaultSettings()
	settings.MajorIterations = 10
	settings.Concurrent = 0
	settings.FuncEvaluations = 100
	//settings.FunctionThreshold = 0.1
	settings.GradientThreshold = 0.1
	settings.Recorder = nil
	//FC := optimize.FunctionConverge{}
	//FC.Absolute = 10
	//FC.Relative = 10
	//FC.Iterations = 10
	//settings.FunctionConverge = &FC
	p := optimize.Problem{Func: fcn, Grad: nil, Hess: nil}
	var p0 []float64
	p0 = append(p0, 0.5)
	for _, n := range tree.PreorderArray() {
		if len(n.CHLD) > 0 {
			p0 = append(p0, n.HEIGHT)
		}
	}
	meth := &optimize.NelderMead{}
	_, err := optimize.Minimize(p, p0, nil, meth)
	if err != nil {
		fmt.Println(err)
	}
	//fmt.Println(res)
}

func OptimizeBranchRates(tree *Node) {
	fcn := func(rates []float64) float64 {
		preNodes := tree.PreorderArray()
		large := 100000000000.0
		for _, i := range rates {
			if i <= 0.0 {
				return large
			}
		}
		var bad bool
		bad = AssignBranchRates(preNodes, rates)
		if bad {
			fmt.Println("HERE")
			return large
		}
		morphLL := RootedLogLikeParallel(tree, true, 4)
		//stratLL := PoissonTreeLoglike(preNodes)
		lnl := morphLL //+ stratLL
		return -lnl
	}
	settings := optimize.Settings{} //DefaultSettings()
	settings.MajorIterations = 10
	settings.Concurrent = 0
	settings.FuncEvaluations = 100
	//settings.FunctionThreshold = 0.1
	settings.GradientThreshold = 0.1
	settings.Recorder = nil
	//FC := optimize.FunctionConverge{}
	//FC.Absolute = 10
	//FC.Relative = 10
	//FC.Iterations = 10
	//settings.FunctionConverge = &FC
	p := optimize.Problem{Func: fcn, Grad: nil, Hess: nil}
	var p0 []float64
	preNodes := tree.PreorderArray()
	for _, node := range preNodes[1:] {
		p0 = append(p0, node.RATE)
	}
	meth := &optimize.NelderMead{}
	_, err := optimize.Minimize(p, p0, nil, meth)
	if err != nil {
		fmt.Println(err)
	}
	//fmt.Println(res)
}
