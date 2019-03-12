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

	cophymaru.InitParallelPRNLEN(nodels)

	cophymaru.ADStratTreeSearch(tree)

}
