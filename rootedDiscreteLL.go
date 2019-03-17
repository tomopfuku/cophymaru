package cophymaru

import (
	"fmt"
	"math"
	"os"
)

func rootedDiscSiteLL(n *Node, nlikes *float64, startFresh bool, site int) {
	for _, chld := range n.CHLD {
		rootedDiscSiteLL(chld, nlikes, startFresh, site)
	}
	nchld := len(n.CHLD)
	//n.CONPRNLEN[site] = n.LSLEN
	log2pi := 1.8378770664093453
	if nchld != 0 {
		if nchld != 2 {
			fmt.Println("This should only be perfomed on fully bifurcating trees/subtrees! Check for multifurcations and singletons.")
			os.Exit(0)
		}
		c0 := n.CHLD[0]
		c1 := n.CHLD[1]
		obsDist := calcSeqDistance(c0, c1)
		curlike := 0.0
		curVar := c0.LSLEN + c1.LSLEN //(c0.CONPRNLEN[site]) + (c1.CONPRNLEN[site])
		//fmt.Println(obsDist, curVar)
		contrast := curVar - obsDist
		curVar = obsDist //1.0
		curlike = ((-0.5)*(log2pi) + (math.Log(curVar))) - (math.Pow(contrast, 2) / (2 * math.Pow(curVar, 2)))
		internalStateProbs(n)
		*nlikes += curlike
		//tempBranchLength := n.CONPRNLEN[site] + (((c0.CONPRNLEN[site]) * (c1.CONPRNLEN[site])) / ((c0.CONPRNLEN[site]) + (c1.CONPRNLEN[site]))) // need to calculate the prune length by adding the averaged lengths of the daughter nodes to the length
		//n.CONPRNLEN[site] = tempBranchLength                                                                                                    // need to calculate the "prune length" by adding the length to the uncertainty
		n.LL[site] = curlike
	}
}

func rootedDiscTreeLike(tree *Node, startFresh bool, jobs <-chan int, results chan<- float64) {
	for site := range jobs {
		tmpll := 0.
		rootedDiscSiteLL(tree, &tmpll, startFresh, site)
		tmpll = tree.LL[site]
		results <- tmpll
	}
}

func DiscLogLikeParallel(tree *Node) (sitelikes float64) {
	startFresh := true
	workers := 1
	nsites := len(tree.CHLD[0].DISCTRT)
	jobs := make(chan int, nsites)
	results := make(chan float64, nsites)
	for w := 0; w < workers; w++ {
		go rootedDiscTreeLike(tree, startFresh, jobs, results)
	}

	for site := 0; site < nsites; site++ {
		jobs <- site
	}
	close(jobs)

	for site := 0; site < nsites; site++ {
		sitelikes += <-results
	}
	return
}
