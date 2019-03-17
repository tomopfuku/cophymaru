package cophymaru

import (
	"fmt"
	"math"
	"os"
)

//IterateDiscBL will iteratively calculate the ML branch lengths for a particular topology
func IterateDiscBL(tree *Node, niter int) (nparam float64) {
	AssertUnrootedTree(tree)
	nodes := tree.PreorderArray()
	//InitMissingValues(nodes)
	itercnt := 0
	for {
		discCalcLengths(tree)
		itercnt++
		if itercnt == niter {
			break
		}
	}
	nparam = 0.0
	for _, n := range nodes {
		if n.ANC == true && n.ISTIP == true {
			continue
		}
		nparam += 1.0
	}
	return
}

//discCalcLengths will perform a single pass of the branch length ML estimation
func discCalcLengths(tree *Node) {
	rnodes := tree.PreorderArray()
	lnode := 0
	for _, cn := range tree.CHLD {
		AncDiscPruneRooted(cn)
	}
	DiscTritomyML(tree)
	/*

		fmt.Println(tree.Phylogram())
		//next := NodeFromLabel("taxon_4", rnodes).PAR
		//old := NodeFromLabel("taxon_4", rnodes).PAR
		tree = next.RerootLS(tree)
		for _, cn := range tree.CHLD {
			AncDiscPruneRooted(cn)
			if cn.NAME == "root" || cn.NAME == "FU" {
				fmt.Println(cn.NAME, cn.DISCTRT[0])
			}

		}
		DiscTritomyML(tree)
			next = NodeFromLabel("taxon_1", rnodes).PAR.PAR
			tree = next.RerootLS(old)
			for _, cn := range tree.CHLD {
				AncDiscPruneRooted(cn)
				//if cn.NAME == "root" || cn.NAME == "FU" {
				fmt.Println(cn.NAME, cn.DISCTRT[0])
				//}
			}
			fmt.Println(tree.Phylogram())
			for _, cn := range tree.CHLD {
				//if cn.NAME == "root" || cn.NAME == "FU" {
				fmt.Println(cn.NAME, cn.DISCTRT[0])
				//}
			}

			internalStateProbs(rnodes[0])
			fmt.Println(old.PAR.NAME, old.CHLD[0].NAME, old.CHLD[1].NAME)
			internalStateProbs(old)

			for _, cn := range tree.CHLD {
				//if cn.NAME == "root" || cn.NAME == "FU" {
				fmt.Println(cn.NAME, cn.DISCTRT[0])
				//}
			}
			DiscTritomyML(tree)
			//fmt.Println(tree.Phylogram())
	*/
	for ind, newroot := range rnodes {
		if len(newroot.CHLD) == 0 {
			continue
		} else if newroot != rnodes[0] {
			tree = newroot.RerootLS(rnodes[lnode])

			lnode = ind
		}
		for _, cn := range tree.CHLD {
			AncDiscPruneRooted(cn)
		}
		DiscTritomyML(tree)
	}
	tree = rnodes[0].RerootLS(tree)
	//fmt.Println(tree.Newick(true))
}

func internalStateProbs(n *Node) {
	c0 := n.CHLD[0]
	c1 := n.CHLD[1]
	totalVar := c0.LSPRNLEN + c1.LSPRNLEN
	var c0weight, c1weight float64
	if c0.LSPRNLEN != 0.0 {
		c0weight = (c0.LSPRNLEN / totalVar)
	} else {
		c0weight = 0.0
	}
	if c1.LSPRNLEN != 0.0 {
		c1weight = (c1.LSPRNLEN / totalVar)
	} else {
		c1weight = 0.0
	}
	for site, c0probs := range c0.DISCTRT {
		c1probs := c1.DISCTRT[site]
		for state, prob0 := range c0probs {
			prob1 := c1probs[state]
			n.DISCTRT[site][state] = (prob0 * c1weight) + (prob1 * c0weight)
		}
	}
}

//AncDiscPruneRooted will prune Disc branch lens and PICs down to a rooted node
//root node should be a real (ie. bifurcating) root
func AncDiscPruneRooted(n *Node) {
	for _, chld := range n.CHLD {
		AncDiscPruneRooted(chld)
	}
	n.LSPRNLEN = n.LSLEN
	nchld := len(n.CHLD)
	if n.ISTIP == false { //&& n.MRK == false {
		if nchld != 2 {
			fmt.Println("This should only be perfomed on fully bifurcating trees/subtrees! Check for multifurcations and singletons.")
			os.Exit(0)
		}
		internalStateProbs(n)
		/*
			c0 := n.CHLD[0]
			c1 := n.CHLD[1]
			if c0.NAME == "taxon_4" || c1.NAME == "taxon_4" {
				n.NAME = "FU"
				//fmt.Println(n.NAME, "HERE", n.DISCTRT[0])
			}
			//bot := ((1.0 / c0.LSPRNLEN) + (1.0 / c1.LSPRNLEN))
				if c0.LSPRNLEN != 0.0 && c1.LSPRNLEN != 0.0 {
					n.LSPRNLEN += (c0.LSPRNLEN * c1.LSPRNLEN) / (c0.LSPRNLEN * c1.LSPRNLEN)
				}
		*/
	}
}

func calcSeqDistance(c0, c1 *Node) (dist float64) {
	dist = 0.0
	for site := range c0.DISCTRT {
		for state := range c0.DISCTRT[site] {
			prob0 := c0.DISCTRT[site][state]
			prob1 := c1.DISCTRT[site][state]
			if prob0 > prob1 {
				dist += (prob0 - prob1)
			}
		}
	}
	//dist = dist / float64(len(c0.DISCTRT))
	return
}

//DiscTritomyML will calculate the MLEs for the branch lengths of a tifurcating 3-taxon tree assuming that direct ancestors may be in the tree
func DiscTritomyML(tree *Node) {
	c0 := tree.CHLD[0]
	c1 := tree.CHLD[1]
	c2 := tree.CHLD[2]

	d01 := calcSeqDistance(c0, c1)
	d02 := calcSeqDistance(c0, c2)
	d12 := calcSeqDistance(c1, c2)
	sum0 := 0.0
	sum1 := 0.0
	sum2 := 0.0
	if c0.ANC == true && c0.ISTIP == true {
		sum1 = math.Pow(d01, 2)
		sum2 = math.Pow(d02, 2)
	} else if c1.ANC == true && c1.ISTIP == true {
		sum2 = math.Pow(d12, 2)
		sum0 = math.Pow(d01, 2)
	} else if c2.ANC == true && c2.ISTIP == true {
		sum1 = math.Pow(d12, 2)
		sum0 = math.Pow(d02, 2)
	} else {
		sum0 = d01 * d02 // (d01 + (d02 - d12)) / 2
		sum1 = d01 * d12 //(d01 + (d12 - d02)) / 2
		sum2 = d02 * d12 //(d02 + (d12 - d01)) / 2
	}

	/*
		if sum0 != 0.0 {
			sum0 = sum0 - (tree.CHLD[0].LSPRNLEN - tree.CHLD[0].LSLEN)
		}
		if sum1 != 0.0 {
			sum1 = sum1 - (tree.CHLD[1].LSPRNLEN - tree.CHLD[1].LSLEN)
		}
		if sum2 != 0.0 {
			sum2 = sum2 - (tree.CHLD[2].LSPRNLEN - tree.CHLD[2].LSLEN)
		}
	*/
	if sum0 < 0. {
		sum0 = 0.0
	}
	if sum1 < 0. {
		sum1 = 0.0
	}
	if sum2 < 0. {
		sum2 = 0.0
	}

	tree.CHLD[0].LSLEN = sum0
	tree.CHLD[1].LSLEN = sum1
	tree.CHLD[2].LSLEN = sum2
}
