package cophymaru

import (
	"fmt"
	"math"
	"math/rand"
)

func ADStratTreeSearch(tree *Node) {
	curbest := checkAllADLL(tree)
	fmt.Println(curbest, tree.Newick(true))
	fmt.Println(tree.Rateogram(true))
}

func checkAllAD(tree *Node) (bestAIC float64) {
	var morphlnl, stratlnl, morphK, stratK, comblnl float64
	var strat_res []float64
	for _, node := range tree.PreorderArray() {
		node.RATE = rand.Float64() //0.9
	}
	morphlnl, morphK, _ = OptimizeBranchRates(tree)
	comblnl, stratK, strat_res = OptimizeMorphStratHeights(tree)
	morphlnl, morphK, _ = OptimizeBranchRates(tree)
	fmt.Println(tree.Newick(true))
	K := morphK + stratK
	lnl := comblnl //morphlnl + stratlnl
	tips := tree.PreorderTips()
	testNodes := candidateAncestors(tips)
	bestAIC = AIC(lnl, K)
	curAIC := AIC(lnl, K)
	fmt.Println(RootedLogLikeParallel(tree, true, 4), curAIC)
	stratlnl = ADPoissonTreeLoglike(tree.PreorderArray(), strat_res[0])
	for _, n := range testNodes {
		bad := MakeAncestor(n)
		if bad {
			continue
		}
		morphlnl, morphK, _ = OptimizeBranchRates(tree)
		//stratlnl, stratK = OptimizeMorphStratHeights(tree)
		stratlnl = ADPoissonTreeLoglike(tree.PreorderArray(), strat_res[0])
		K = morphK + stratK
		lnl = morphlnl + stratlnl
		curAIC = AIC(lnl, K)
		fmt.Println(n.NAME, RootedLogLikeParallel(tree, true, 4), morphlnl, stratlnl, lnl, curAIC, bestAIC)
		if (curAIC + 0.5) < bestAIC {
			bestAIC = curAIC
		} else {
			UnmakeAncestor(n)
		}
	}
	return
}

func checkAllADLL(tree *Node) (bestLL float64) {
	var morphlnl, stratlnl, comblnl float64
	var strat_res []float64
	for _, node := range tree.PreorderArray() {
		node.RATE = rand.Float64() //0.9
	}
	morphlnl, _, _ = OptimizeBranchRates(tree)
	comblnl, _, strat_res = OptimizeMorphStratHeights(tree)
	morphlnl, _, _ = OptimizeBranchRates(tree)
	fmt.Println(tree.Newick(true))
	lnl := comblnl //morphlnl + stratlnl
	tips := tree.PreorderTips()
	testNodes := candidateAncestors(tips)
	bestLL = lnl
	curLL := lnl
	fmt.Println(RootedLogLikeParallel(tree, true, 4), curLL)
	stratlnl = ADPoissonTreeLoglike(tree.PreorderArray(), strat_res[0])
	for _, n := range testNodes {
		bad := MakeAncestor(n)
		if bad {
			continue
		}
		morphlnl, _, _ = OptimizeBranchRates(tree)
		//stratlnl, stratK = OptimizeMorphStratHeights(tree)
		stratlnl = ADPoissonTreeLoglike(tree.PreorderArray(), strat_res[0])
		lnl = morphlnl + stratlnl
		curLL = lnl
		fmt.Println(n.NAME, morphlnl, stratlnl, lnl, curLL, bestLL, math.Exp((curLL-bestLL)/2))
		if curLL-0.5 > bestLL {
			bestLL = curLL
		} else {
			UnmakeAncestor(n)
		}
	}
	return
}

func candidateAncestors(tips []*Node) (anc []*Node) {
	for _, n := range tips {
		if n.FAD > n.GetSib().FAD && n.PAR.NAME != "root" {
			anc = append(anc, n)
		}
	}
	return
}
