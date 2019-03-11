package main

import (
	"cophymaru"
	"flag"
	"fmt"
	"math/rand"
	"runtime"
	"time"
)

func main() {
	treeArg := flag.String("t", "", "input tree")
	stratArg := flag.String("r", "", "temporal ranges")
	traitArg := flag.String("m", "", "continuous traits")
	threadArg := flag.Int("T", 1, "maximum number of cores to use during run")
	flag.Parse()
	runtime.GOMAXPROCS(*threadArg)
	rand.Seed(time.Now().UTC().UnixNano())

	nwk := cophymaru.ReadLine(*treeArg)[0]
	tree := cophymaru.ReadTree(nwk)
	nodels := tree.PreorderArray()
	traits, _, ntraits := cophymaru.ReadContinuous(*traitArg)
	cophymaru.MapContinuous(tree, traits, ntraits)
	cophymaru.ReadStrat(*stratArg, nodels)
	cophymaru.MakeStratHeights(tree)
	fmt.Println(tree.Newick(true))
	for _, node := range nodels {
		node.RATE = 0.4
	}
	cophymaru.InitParallelPRNLEN(nodels)
	start := time.Now()
	var morphlnl, stratlnl, morphK, stratK float64
	morphlnl, morphK = cophymaru.OptimizeBranchRates(tree)
	stratlnl, stratK = cophymaru.OptimizeMorphStratHeights(tree)
	fmt.Println(tree.Newick(true))
	K := morphK + stratK
	lnl := morphlnl + stratlnl
	fmt.Println("bifurcating AIC:", cophymaru.AIC(lnl, K))
	cophymaru.MakeAncestorLabel("H_hei", tree.PreorderArray())
	cophymaru.MakeAncestorLabel("A_africanus", tree.PreorderArray())
	cophymaru.MakeAncestorLabel("H_erg", tree.PreorderArray())
	cophymaru.MakeAncestorLabel("H_erectus", tree.PreorderArray())
	cophymaru.MakeAncestorLabel("P_aet", tree.PreorderArray())
	cophymaru.MakeAncestorLabel("H_rud", tree.PreorderArray())
	//cophymaru.MakeAncestorLabel("Pan_M", tree.PreorderArray())
	fmt.Println(tree.Newick(true))
	morphlnl, morphK = cophymaru.OptimizeBranchRates(tree)
	//stratlnl, stratK = cophymaru.OptimizeMorphStratHeights(tree)
	K = morphK + stratK
	lnl = morphlnl + stratlnl
	fmt.Println("AD AIC:", cophymaru.AIC(lnl, K))
	elapsed := time.Since(start)
	fmt.Println(tree.Newick(true))
	fmt.Println(elapsed)
	//fmt.Println(tree.Rateogram(true))
}
