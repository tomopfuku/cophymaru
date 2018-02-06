package main

/*
this is going to test BM pruning functions
*/

import (
	"cophymaru"
	"flag"
	"fmt"
	"strings"
)

func postorder(curnode *cophymaru.Node) {
	for _, chld := range curnode.CHLD {
		postorder(chld)
	}
	fmt.Println(curnode.NAME, curnode.CONTRT, curnode.LEN)
}

func main() {
	treeArg := flag.String("t", "", "input tree")
	traitArg := flag.String("m", "", "continuous traits")
	iterArg := flag.Int("i", 1, "num iterations for branch length calculation")
	genArg := flag.Int("gen", 5000, "number of MCMC generations to run")
	brPrior := flag.String("blpr", "0", "\tspecifies the prior to place on branch lengths. \n\tOptions:\n\n\t\t0    Flat\n\t\t1    Exponential (mean = 10)\n")
	fosArg := flag.String("fos", "", "file containing names of fossil tips\n\ttips should be comma-separated")
	startArg := flag.String("st", "0", "\tSpecify whether to use ML or random starting tree and branch lengths\n\t\t0    ML\n\t\t1    random\n")
	flag.Parse()
	//var ntax,ntraits int
	nwk := cophymaru.ReadLine(*treeArg)[0]
	tree := cophymaru.ReadTree(nwk)
	fmt.Println(tree.Newick(true))
	traits, ntax, ntraits := cophymaru.ReadContinuous(*traitArg)
	fmt.Println(ntax)
	cophymaru.MapContinuous(tree, traits, ntraits)
	//cophymaru.IterateBMLengths(tree,*iterArg)
	var fosSlice []string // read in fossil names from command line
	for _, i := range strings.Split(*fosArg, ",") {
		fosSlice = append(fosSlice, i)
	}
	if *startArg == "0" {
		starttr, startll := cophymaru.InsertFossilTaxa(tree, traits, fosSlice, *iterArg)
		fmt.Println("STARTING ML TREE:\n", starttr, "\n\nSTARTING MCMC WITH LOG-LIKELIHOOD ", startll)
	} else if *startArg == "1" {
		cophymaru.InsertFossilTaxaRandom(tree, traits, fosSlice, *iterArg)
	}
	cophymaru.MCMC(tree, *genArg, fosSlice, "tmp/test.t", "tmp/test.mcmc", *brPrior)
}
