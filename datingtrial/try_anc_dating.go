package main

import (
	"cophymaru"
	"flag"
	"fmt"
	"math/rand"
	"os"
	"runtime"
	"strings"
	"time"
)

func main() {
	treeArg := flag.String("t", "", "input tree")
	stratArg := flag.String("r", "", "temporal ranges")
	traitArg := flag.String("m", "", "continuous traits")
	threadArg := flag.Int("T", 1, "maximum number of cores to use during run")
	outgrpArg := flag.String("g", "", "Specify the outgroup taxon/taxa (separated by commas if there are multiple (example: Alouatta,Cercopithecus)")
	ancArg := flag.String("a", "", "temporal ranges")

	flag.Parse()
	runtime.GOMAXPROCS(*threadArg)
	rand.Seed(time.Now().UTC().UnixNano())
	nwk := cophymaru.ReadLine(*treeArg)[0]
	tree := cophymaru.ReadTree(nwk)
	tree.SetOutgroup(strings.Split(*outgrpArg, ","))
	for _, n := range tree.PreorderArray() {
		n.Len = rand.Float64()
		n.LSLen = rand.Float64()
	}
	nodels := tree.PreorderArray()
	traits, _, ntraits := cophymaru.ReadContinuous(*traitArg)
	cophymaru.MapContinuous(tree, traits, ntraits)
	cophymaru.InitMissingValues(tree.PreorderArray())
	cophymaru.ReadStrat(*stratArg, nodels)
	cophymaru.MakeStratHeights(tree)
	cophymaru.MakeAncestorLabel(*ancArg, tree.PreorderArray())
	cophymaru.MakeAncestorLabel("H_rud", tree.PreorderArray())
	cophymaru.MakeAncestorLabel("H_hei", tree.PreorderArray())
	//cophymaru.MakeAncestorLabel("H_erg", tree.PreorderArray())
	//cophymaru.MakeAncestorLabel("H_erectus", tree.PreorderArray())
	cophymaru.MakeAncestorLabel("A_africanus", tree.PreorderArray())
	//cophymaru.MakeAncestorLabel("H_sap", tree.PreorderArray())
	_, _ = cophymaru.OptimizeGlobalRateHeights(tree, 2.8)

	fmt.Println(tree.Newick(true))
	sublen := tree.Unroot()

	cophymaru.AncMissingTraitsEM(tree, 1000)

	tree.Root(sublen)
	tree.CalcBranchRates()
	//cophymaru.InitParallelPRNLen(nodels)
	fmt.Println(tree.Phylogram())

	_, _, _ = cophymaru.OptimizeMorphStratHeights(tree, 2.8)
	fmt.Println(tree.Newick(true))
	os.Exit(0)

	fmt.Println("rateogram:", tree.Rateogram())
	fmt.Println(cophymaru.RootedLogLikeParallel(tree, true, 4))
	fmt.Println(tree.Phylogram())
	fmt.Println("rateogram:", tree.Rateogram())
	fmt.Println(tree.Newick(true))

	os.Exit(0)

	cophymaru.ADStratTreeSearch(tree)

}
