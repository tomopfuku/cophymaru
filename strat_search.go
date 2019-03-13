package cophymaru

import (
	"fmt"
	"math"
	"math/rand"
)

func ADStratTreeSearch(tree *Node) {
	curbest := checkAllADLL(tree)
	fmt.Println(curbest, tree.Newick(true))
	fmt.Println(tree.Rateogram())
}

func checkAllAD(tree *Node) (bestAIC float64) {
	var morphlnl, stratlnl, morphK, stratK, comblnl float64
	for _, node := range tree.PreorderArray() {
		node.RATE = rand.Float64() //0.9
	}
	morphlnl, morphK, _ = OptimizeBranchRates(tree)
	_, lam := OptimizePreservationLam(tree)
	comblnl, stratK, _ = OptimizeMorphStratHeights(tree, lam)
	morphlnl, morphK, _ = OptimizeBranchRates(tree)
	fmt.Println(tree.Newick(true))
	K := morphK + stratK
	lnl := comblnl //morphlnl + stratlnl
	tips := tree.PreorderTips()
	testNodes := candidateAncestors(tips)
	bestAIC = AIC(lnl, K)
	curAIC := AIC(lnl, K)
	bifAIC := AIC(lnl, K)
	stratlnl = ADPoissonTreeLoglike(tree.PreorderArray(), lam)
	fmt.Println(RootedLogLikeParallel(tree, true, 4), stratlnl, comblnl, curAIC, lam)
	anclikes := make(map[string]float64)
	for _, n := range testNodes {
		bad := MakeAncestor(n)
		if bad {
			continue
		}
		morphlnl, _, _ = OptimizeBranchRates(tree)
		//morphlnl, _, _ = OptimizeLocalRatesAD(n)
		morphlnl = RootedLogLikeParallel(tree, true, 4)

		//stratlnl, stratK = OptimizeMorphStratHeights(tree)
		stratlnl = ADPoissonTreeLoglike(tree.PreorderArray(), lam)
		K = (morphK - 1) + stratK
		lnl = morphlnl + stratlnl
		curAIC = AIC(lnl, K)
		//fmt.Println(n.NAME, RootedLogLikeParallel(tree, true, 4), morphlnl, stratlnl, lnl, curAIC, bestAIC)
		aicW := 0.0
		if curAIC < bifAIC {
			rellike := math.Exp((bifAIC - curAIC) / 2.0)
			//bifrellike := math.Exp((bifAIC - bifAIC) / 2.0)
			sum := rellike + 1.0 //bifrellike
			aicW = rellike / sum
		}
		fmt.Println(n.NAME, morphlnl, stratlnl, curAIC, bifAIC, aicW)
		if (curAIC) < bestAIC {
			bestAIC = curAIC
			anclikes[n.NAME] = curAIC
			UnmakeAncestor(n)
		} else {
			UnmakeAncestor(n)
		}
	}
	for _, n := range testNodes {
		UnmakeAncestor(n)
	}

	/*
		rellikes := make(map[string]float64)
		sum := 0.0
		for k, _ := range anclikes {
			_ = MakeAncestorLabel(k, tree.PreorderArray())
			morphlnl, morphK, _ = OptimizeBranchRates(tree)
			//stratlnl, stratK = OptimizeMorphStratHeights(tree)
			stratlnl = ADPoissonTreeLoglike(tree.PreorderArray(), strat_res[0])
			K = morphK + stratK
			lnl = morphlnl + stratlnl
			curAIC = AIC(lnl, K)
			//curAIC = v
			rellike := math.Exp((bestAIC - curAIC) / 2.0)

			rellikes[k] = rellike
			sum += rellike
			//if rellike > 2000000 {
			UnmakeAncestorLabel(k, tree.PreorderArray())
			//	}
			fmt.Println(k, curAIC, bifAIC, bestAIC, lnl, rellike)
		}
		for k, v := range rellikes {
			aicW := (v / sum) * 100
			fmt.Println(k, aicW)
			if aicW > 5.0 {
				MakeAncestorLabel(k, tree.PreorderArray())
			}
		}
	*/
	return
}

func checkAllADLL(tree *Node) (bestLL float64) {
	var morphlnl, stratlnl, comblnl float64
	for _, node := range tree.PreorderArray() {
		node.RATE = rand.Float64() //0.9
	}
	_, lam := OptimizePreservationLam(tree)
	fmt.Println(lam)
	morphlnl, _ = OptimizeGlobalRateHeights(tree, lam)
	comblnl, _, _ = OptimizeMorphStratHeights(tree, lam)
	morphlnl = RootedLogLikeParallel(tree, true, 4)
	//morphlnl, _, _ = OptimizeBranchRates(tree)
	fmt.Println(tree.Newick(true))
	lnl := comblnl //morphlnl + stratlnl
	tips := tree.PreorderTips()
	testNodes := candidateAncestors(tips)
	bestLL = lnl
	curLL := lnl
	bifLL := lnl
	bifMorphLL := morphlnl
	//fmt.Println(RootedLogLikeParallel(tree, true, 4), curLL)
	nodes := tree.PreorderArray()
	stratlnl = ADPoissonTreeLoglike(nodes, lam)
	var rellike float64
	ancsupport := make(map[string]float64)
	for _, n := range testNodes {
		bad := MakeAncestor(n)
		if bad {
			continue
		}
		//for _, node := range tree.PreorderArray() {
		//	node.RATE = rand.Float64() //0.9
		//}
		morphlnl, _, _ = OptimizeBranchRates(tree)
		//morphlnl, _, _ = OptimizeLocalRatesAD(n)
		morphlnl = RootedLogLikeParallel(tree, true, 4)
		//stratlnl, stratK = OptimizeMorphStratHeights(tree)
		stratlnl = ADPoissonTreeLoglike(nodes, lam)
		lnl = morphlnl + stratlnl
		curLL = lnl
		//rellike = math.Exp(curLL) / (math.Exp(curLL) + math.Exp(bifLL))
		rellike = math.Exp(curLL - (curLL + math.Log1p(math.Exp(bifLL-curLL))))
		morphrellike := math.Exp(morphlnl - (morphlnl + math.Log1p(math.Exp(bifMorphLL-morphlnl))))
		ancsupport[n.NAME] = rellike
		fmt.Println(n.NAME, rellike, morphlnl, bifMorphLL, morphrellike)
		if curLL > bestLL {
			bestLL = curLL
			UnmakeAncestor(n)
		} else {
			UnmakeAncestor(n)
		}
	}
	for k, v := range ancsupport {
		if v > 0.7 {
			MakeAncestorLabel(k, nodes)
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
