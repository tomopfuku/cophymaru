package main

/*
this is going to test BM pruning functions
*/

import (
	"cophy"
	"flag"
	"fmt"
	"strings"
)

func postorder(curnode *cophy.Node) {
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
	brPrior := flag.String("blpr", "0", "specifies the prior to place on branch lengths. \nOptions:\n\n0    Flat\n1    Exponential (mean = 10)\n")
	fosArg := flag.String("fos", "", "file containing names of fossil tips\ntips should be comma-separated")
	startArg := flag.String("st", "0", "Specify whether to use ML or random starting tree and branch lengths\n0    ML\n1    random\n")
	flag.Parse()
	//var ntax,ntraits int
	nwk := cophy.ReadLine(*treeArg)[0]
	tree := cophy.ReadTree(nwk)
	fmt.Println(tree.Newick(true))
	traits, ntax, ntraits := cophy.ReadContinuous(*traitArg)
	fmt.Println(ntax)
	cophy.MapContinuous(tree, traits, ntraits)
	//cophy.IterateBMLengths(tree,*iterArg)
	var fosSlice []string // read in fossil names from command line
	for _, i := range strings.Split(*fosArg, ",") {
		fosSlice = append(fosSlice, i)
	}
	if *startArg == "0" {
		starttr, startll := cophy.InsertFossilTaxa(tree, traits, fosSlice, *iterArg)
		fmt.Println("STARTING ML TREE:\n", starttr, "\n\nSTARTING MCMC WITH LOG-LIKELIHOOD ", startll)
	} else if *startArg == "1" {
		cophy.InsertFossilTaxaRandom(tree, traits, fosSlice, *iterArg)
	}
	cophy.MCMC(tree, *genArg, fosSlice, "tmp/test.t", "tmp/test.mcmc", *brPrior)
}
