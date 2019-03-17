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
	traitArg := flag.String("m", "", "discrete traits")

	flag.Parse()
	rand.Seed(time.Now().UTC().UnixNano())
	nwk := cophymaru.ReadLine(*treeArg)[0]
	tree := cophymaru.ReadTree(nwk)
	for _, n := range tree.PreorderArray() {
		n.LSLEN = n.LEN //rand.Float64()
	}
	//nodels := tree.PreorderArray()
	traits := cophymaru.ReadDiscrete(*traitArg)
	cophymaru.MapDiscrete(tree, traits)
	tree.SetOutgroup([]string{"taxon_1"})
	addlen := tree.Unroot()
	_ = cophymaru.IterateDiscBL(tree, 1)
	tree.Root(addlen)
	fmt.Println(tree.Phylogram())

	fmt.Println(cophymaru.DiscLogLikeParallel(tree))
}
