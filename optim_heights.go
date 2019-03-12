package cophymaru

import (
	"fmt"

	"gonum.org/v1/gonum/optimize"
)

func OptimizeMorphStratHeights(tree *Node) (float64, float64, []float64) {
	fcn := func(heights []float64) float64 {
		lam := 2.4
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
		stratLL := ADPoissonTreeLoglike(preNodes, lam)
		lnl := morphLL + stratLL
		return -lnl
	}
	settings := optimize.Settings{} //DefaultSettings()
	settings.MajorIterations = 10
	settings.Concurrent = 0
	settings.FuncEvaluations = 100
	settings.GradientThreshold = 0.1
	settings.Recorder = nil
	p := optimize.Problem{Func: fcn, Grad: nil, Hess: nil}
	var p0 []float64
	preNodes := tree.PreorderArray()
	for _, n := range preNodes {
		if n.ISTIP == false && n.ANC == false {
			p0 = append(p0, n.HEIGHT)
		}
	}
	meth := &optimize.NelderMead{}
	res, err := optimize.Minimize(p, p0, nil, meth)
	if err != nil {
		fmt.Println(err)
	}
	AssignInternalNodeHeights(preNodes, res.X)
	var retparams []float64
	retparams = append(retparams, 1.0)
	for _, bl := range res.X {
		retparams = append(retparams, bl)
	}
	return -res.F, float64(len(res.X)), retparams //res.X

}

func OptimizeLamMorphStratHeights(tree *Node) (float64, float64, []float64) {
	fcn := func(params []float64) float64 {
		lam := params[0]
		heights := params[1:]
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
		stratLL := ADPoissonTreeLoglike(preNodes, lam)
		lnl := morphLL + stratLL
		return -lnl
	}
	settings := optimize.Settings{} //DefaultSettings()
	settings.MajorIterations = 10
	settings.Concurrent = 0
	settings.FuncEvaluations = 100
	settings.GradientThreshold = 0.1
	settings.Recorder = nil
	p := optimize.Problem{Func: fcn, Grad: nil, Hess: nil}
	var p0 []float64
	p0 = append(p0, 1.0)
	preNodes := tree.PreorderArray()
	for _, n := range preNodes {
		if n.ISTIP == false && n.ANC == false {
			p0 = append(p0, n.HEIGHT)
		}
	}
	meth := &optimize.NelderMead{}
	res, err := optimize.Minimize(p, p0, nil, meth)
	if err != nil {
		fmt.Println(err)
	}
	AssignInternalNodeHeights(preNodes, res.X[1:])
	return -res.F, float64(len(res.X)), res.X

}

func OptimizeGlobalRateHeights(tree *Node) (float64, []float64) {
	fcn := func(params []float64) float64 {
		rate := params[0]
		lam := params[1]
		heights := params[2:]
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
		stratLL := ADPoissonTreeLoglike(preNodes, lam)
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
	p0 = append(p0, 1.0)
	preNodes := tree.PreorderArray()
	for _, n := range preNodes {
		if len(n.CHLD) > 0 && n.ANC == false {
			p0 = append(p0, n.HEIGHT)
		}
	}
	meth := &optimize.NelderMead{}
	res, err := optimize.Minimize(p, p0, nil, meth)
	if err != nil {
		fmt.Println(err)
	}
	AssignGlobalRate(preNodes, res.X[0])
	AssignInternalNodeHeights(preNodes, res.X[2:])
	return -res.F, res.X
}

func OptimizeBranchRates(tree *Node) (float64, float64, []float64) {
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
	settings.GradientThreshold = 0.1
	settings.Recorder = nil
	p := optimize.Problem{Func: fcn, Grad: nil, Hess: nil}
	var p0 []float64
	preNodes := tree.PreorderArray()
	for _, node := range preNodes[1:] {
		if node.ANC == true && node.ISTIP == true {
			continue
		}
		p0 = append(p0, node.RATE)
	}
	meth := &optimize.NelderMead{}
	res, err := optimize.Minimize(p, p0, nil, meth)
	if err != nil {
		fmt.Println(err)
	}
	AssignBranchRates(preNodes, res.X)
	return -res.F, float64(len(res.X)), res.X
}
