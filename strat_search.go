package cophymaru

import (
	"fmt"
	"math"
)

func ADStratTreeSearch(tree *Node) {
	//curbest := checkAllADLL(tree)
	curbest := checkAllADAIC(tree)
	//Make2BudAncestor(NodeFromLabel("H_erg", tree.PreorderArray()).GetSib())
	fmt.Println(curbest, tree.Newick(true))

	/*sublen := tree.Unroot()
	morphK := AncMissingTraitsEM(tree, 100)
	tree.Root(sublen)
	morphlnl := RootedLogLikeParallel(tree, true, 4)
	fmt.Println(AIC(morphlnl, morphK))
	*/
	//fmt.Println(tree.Rateogram())
}

func checkAllADAIC(tree *Node) (bestLL float64) {
	var morphlnl, stratlnl, comblnl float64
	_, lam := OptimizePreservationLam(tree)
	morphlnl, _ = OptimizeGlobalRateHeights(tree, lam)
	sublen := tree.Unroot()
	bifMorphK := AncMissingTraitsEM(tree, 100)
	tree.Root(sublen)
	tree.CalcBranchRates()
	var stratK float64
	comblnl, stratK, _ = OptimizeMorphStratHeights(tree, lam)
	stratK += 1.0
	morphlnl = RootedLogLikeParallel(tree, true, 4)
	//morphlnl, _, _ = OptimizeBranchRates(tree)
	fmt.Println(tree.Newick(true))
	lnl := comblnl //morphlnl + stratlnl
	tips := tree.PreorderTips()
	testNodes := candidateAncestors(tips)
	curLL := lnl
	bestLL = AIC(curLL, bifMorphK+stratK)
	bifLL := bestLL
	bifMorphLL := morphlnl
	bifMorphAIC := AIC(bifMorphLL, bifMorphK)
	//fmt.Println(RootedLogLikeParallel(tree, true, 4), curLL)
	nodes := tree.PreorderArray()
	stratlnl = ADPoissonTreeLoglike(nodes, lam)
	fmt.Println(morphlnl, stratlnl, comblnl)
	var rellike float64
	ancsupport := make(map[string]float64)
	for _, n := range testNodes {
		bad := MakeAncestor(n)
		if bad {
			continue
		}
		sublen := tree.Unroot()
		morphK := AncMissingTraitsEM(tree, 100)
		tree.Root(sublen)
		tree.CalcBranchRates()
		//fmt.Println(tree.Phylogram())
		//comblnl, _, _ = OptimizeMorphStratHeights(tree, lam)
		morphlnl = RootedLogLikeParallel(tree, true, 4)
		stratlnl = ADPoissonTreeLoglike(tree.PreorderArray(), lam)
		lnl = morphlnl + stratlnl
		curLL = AIC(lnl, morphK+stratK)
		//rellike = math.Exp(curLL) / (math.Exp(curLL) + math.Exp(bifLL))
		rellike = math.Exp(curLL - (curLL + math.Log1p(math.Exp(bifLL-curLL))))
		morphrellike := math.Exp(morphlnl - (morphlnl + math.Log1p(math.Exp(bifMorphLL-morphlnl))))
		ancsupport[n.NAME] = rellike
		fmt.Println(n.NAME, curLL, bifLL, AIC(morphlnl, morphK), bifMorphAIC, morphrellike)
		if curLL < bestLL {
			bestLL = curLL
			//UnmakeAncestor(n)

		} else {
			UnmakeAncestor(n)
		}
	}
	return
}

func checkAllADLL(tree *Node) (bestLL float64) {
	var morphlnl, stratlnl, comblnl float64
	_, lam := OptimizePreservationLam(tree)
	morphlnl, _ = OptimizeGlobalRateHeights(tree, lam)
	sublen := tree.Unroot()
	bifMorphK := AncMissingTraitsEM(tree, 100)
	tree.Root(sublen)
	tree.CalcBranchRates()
	var stratK float64
	comblnl, stratK, _ = OptimizeMorphStratHeights(tree, lam)
	stratK += 1.0
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
	fmt.Println(morphlnl, stratlnl, comblnl)
	var rellike float64
	ancsupport := make(map[string]float64)
	for _, n := range testNodes {
		bad := MakeAncestor(n)
		if bad {
			continue
		}
		sublen := tree.Unroot()
		morphK := AncMissingTraitsEM(tree, 100)
		tree.Root(sublen)
		tree.CalcBranchRates()
		//fmt.Println(tree.Phylogram())
		comblnl, _, _ = OptimizeMorphStratHeights(tree, lam)
		morphlnl = RootedLogLikeParallel(tree, true, 4)
		stratlnl = ADPoissonTreeLoglike(tree.PreorderArray(), lam)
		lnl = morphlnl + stratlnl
		curLL = lnl
		//rellike = math.Exp(curLL) / (math.Exp(curLL) + math.Exp(bifLL))
		rellike = math.Exp(curLL - (curLL + math.Log1p(math.Exp(bifLL-curLL))))
		morphrellike := math.Exp(morphlnl - (morphlnl + math.Log1p(math.Exp(bifMorphLL-morphlnl))))
		ancsupport[n.NAME] = rellike
		fmt.Println(n.NAME, lnl, bifLL, morphlnl, bifMorphLL, AIC(curLL, morphK+stratK), AIC(bifLL, bifMorphK+stratK), morphrellike)
		if curLL > bestLL {
			bestLL = curLL
			//UnmakeAncestor(n)

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
