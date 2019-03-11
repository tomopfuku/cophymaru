package cophymaru

import (
	"fmt"
	"math"
	"os"
)

func rootedSiteLL(n *Node, nlikes *float64, startFresh bool, site int) {
	for _, chld := range n.CHLD {
		rootedSiteLL(chld, nlikes, startFresh, site)
	}
	nchld := len(n.CHLD)
	if n.MRK == true {
		if startFresh == false {
			if nchld != 0 {
				*nlikes += n.LL[site]
			}
		}
	}
	if n.MRK == false || startFresh == true {
		n.CONPRNLEN[site] = n.LEN
		log2pi := 1.8378770664093453
		if nchld != 0 {
			if nchld != 2 {
				fmt.Println("This BM pruning algorithm should only be perfomed on fully bifurcating trees/subtrees! Check for multifurcations and singletons.")
				os.Exit(0)
			}
			c0 := n.CHLD[0]
			c1 := n.CHLD[1]
			curlike := float64(0.0)
			var tempChar float64
			curVar := (c0.CONPRNLEN[site] * c0.RATE) + (c1.CONPRNLEN[site] * c1.RATE)
			contrast := c0.CONTRT[site] - c1.CONTRT[site]
			curlike += ((-0.5) * ((log2pi) + (math.Log(curVar)) + (math.Pow(contrast, 2) / (curVar))))
			tempChar = (((c0.CONPRNLEN[site] * c0.RATE) * c1.CONTRT[site]) + ((c1.CONPRNLEN[site] * c1.RATE) * c0.CONTRT[site])) / (curVar)
			n.CONTRT[site] = tempChar
			*nlikes += curlike
			tempBranchLength := n.CONPRNLEN[site] + (((c0.CONPRNLEN[site] * c0.RATE) * (c1.CONPRNLEN[site] * c1.RATE)) / ((c0.RATE * c0.CONPRNLEN[site]) + (c1.RATE * c1.CONPRNLEN[site]))) // need to calculate the prune length by adding the averaged lengths of the daughter nodes to the length
			n.CONPRNLEN[site] = tempBranchLength                                                                                                                                            // need to calculate the "prune length" by adding the length to the uncertainty
			n.LL[site] = curlike
			//n.MRK = true
		}
	}
}

func rootedTreeLike(tree *Node, startFresh bool, jobs <-chan int, results chan<- float64) {
	for site := range jobs {
		tmpll := 0.
		rootedSiteLL(tree, &tmpll, startFresh, site)
		tmpll = tree.LL[site]
		results <- tmpll
	}
}

func RootedLogLikeParallel(tree *Node, startFresh bool, workers int) (sitelikes float64) {
	nsites := len(tree.CHLD[0].CONTRT)
	jobs := make(chan int, nsites)
	results := make(chan float64, nsites)
	for w := 0; w < workers; w++ {
		go rootedTreeLike(tree, startFresh, jobs, results)
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
