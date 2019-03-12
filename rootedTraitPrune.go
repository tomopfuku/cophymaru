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

func rootedMissingSiteLL(n *Node, nlikes *float64, startFresh bool, site int) {
	for _, chld := range n.CHLD {
		rootedMissingSiteLL(chld, nlikes, startFresh, site)
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
			if c0.MIS[site] == false && c1.MIS[site] == false {
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
			} else if c0.MIS[site] == true && c1.MIS[site] == false {
				n.CONPRNLEN[site] += (c1.CONPRNLEN[site])
				lensum := n.CONPRNLEN[site] + c1.CONPRNLEN[site]
				nprop := n.CONPRNLEN[site] / lensum
				descprop := c1.CONPRNLEN[site] / lensum
				n.RATE = (n.RATE * nprop) + (c1.RATE * descprop)
				n.CONTRT[site] = c1.CONTRT[site]
				n.LL[site] = 0.0
			} else if c1.MIS[site] == true && c0.MIS[site] == false {
				n.CONPRNLEN[site] += c0.CONPRNLEN[site]
				lensum := n.CONPRNLEN[site] + c0.CONPRNLEN[site]
				nprop := n.CONPRNLEN[site] / lensum
				descprop := c0.CONPRNLEN[site] / lensum
				n.RATE = (n.RATE * nprop) + (c0.RATE * descprop)
				n.CONTRT[site] = c0.CONTRT[site]
				n.LL[site] = 0.0
			}
		}
	}
}

func rootedTreeLike(tree *Node, startFresh bool, jobs <-chan int, results chan<- float64) {
	for site := range jobs {
		tmpll := 0.
		rootedMissingSiteLL(tree, &tmpll, startFresh, site)
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
