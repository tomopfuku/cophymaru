package main

import (
	"cophymaru"
	"flag"
	"fmt"
	"math/rand"
	"time"
)

func main() {
	treeArg := flag.String("t", "", "input tree")
	stratArg := flag.String("r", "", "temporal ranges")
	traitArg := flag.String("m", "", "continuous traits")
	flag.Parse()
	rand.Seed(time.Now().UTC().UnixNano())

	nwk := cophymaru.ReadLine(*treeArg)[0]
	tree := cophymaru.ReadTree(nwk)
	nodels := tree.PreorderArray()
	traits, _, ntraits := cophymaru.ReadContinuous(*traitArg)
	cophymaru.MapContinuous(tree, traits, ntraits)
	cophymaru.ReadStrat(*stratArg, nodels)
	cophymaru.MakeStratHeights(tree)
	fmt.Println(tree.Newick(true))
	stratLL := cophymaru.PoissonTreeLoglike(tree.PreorderArray())
	fmt.Println(cophymaru.PoissonTreeLoglike(tree.PreorderArray()))
	for _, node := range nodels {
		node.RATE = 1.0
	}
	cophymaru.InitParallelPRNLEN(nodels)
	fmt.Println(cophymaru.RootedLogLikeParallel(tree, true, 2))
	traitLL := cophymaru.RootedLogLikeParallel(tree, true, 2)
	fmt.Println(stratLL + traitLL)
}
