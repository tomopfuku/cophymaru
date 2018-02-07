package cophymaru

import (
	"math/rand"
	"time"
)

//TODO: this has a bug somewhere. need to track it down at some point.

//InsertFossilTaxa finds the ML placement of a set of fossils via stepwise addition
func InsertFossilTaxa(tree *Node, traits map[string][]float64, fosNms []string, iter int, missing bool) (besttr string, bestll float64) {
	nodes := tree.PreorderArray()
	//TODO: get this missing data stuff up and running
	//trMeans := CalcSiteMeans(nodes)
	//MakeMissingMeansArray(nodes)
	for _, curFos := range fosNms {
		ftip := new(Node)
		ftip.NAME = curFos
		ftip.LEN = 0.01
		ftip.CONTRT = traits[curFos]
		//MakeMissingDataSlice(ftip)
		//MakeMissingMeansTip(ftip, trMeans) //replace missing traits with the mean value across all sample traits
		newpar := new(Node)
		newpar.LEN = 0.01
		newpar.AddChild(ftip)
		for range ftip.CONTRT {
			newpar.CONTRT = append(newpar.CONTRT, float64(0.0))
		}
		bestll = -1000000000000.0
		besttr = ""
		curll := float64(0.0)
		var bestPlace *Node
		for _, n := range nodes[1:] {
			GraftFossilTip(newpar, n)
			//fmt.Println(tree.Newick(true))
			//start := time.Now()
			if missing == false {
				IterateBMLengths(tree, iter)
			} else if missing == true {
				MissingTraitsEM(tree, iter)
			}
			//end := time.Now()
			//fmt.Println(end.Sub(start))
			curll = CalcUnrootedLogLike(tree)
			if curll > bestll {
				bestll = curll
				besttr = tree.Newick(true)
				bestPlace = n
			}
			//fmt.Println(curll)
			//fmt.Println(tree.Newick(true))
			PruneTip(newpar, n)
		}
		//fmt.Println(bestll,besttr)
		GraftFossilTip(newpar, bestPlace)
		IterateBMLengths(tree, iter)
		//nodes = tree.PreorderArray() // need to reinitialize the node list to include the now-placed fossil
		//fmt.Println(tree.Newick(true))
	}
	return
}

//InsertFossilTaxaRandom will randomly insert all of the fossil taxa in a dadaset
func InsertFossilTaxaRandom(tree *Node, traits map[string][]float64, fosNms []string, iter int, missing bool) {
	nodes := tree.PreorderArray()
	for _, curFos := range fosNms {
		ftip := new(Node)
		ftip.NAME = curFos
		ftip.LEN = 0.01
		ftip.CONTRT = traits[curFos]
		MakeMissingDataSlice(ftip)
		//MakeMissingMeansTip(ftip, trMeans) //replace missing traits with the mean value across all sample traits
		newpar := new(Node)
		newpar.LEN = 0.01
		newpar.AddChild(ftip)
		for range ftip.CONTRT {
			newpar.CONTRT = append(newpar.CONTRT, float64(0.0))
		}
		reattach := randomNode(nodes[1:])
		GraftFossilTip(newpar, reattach)
		nodes = tree.PreorderArray() // need to reinitialize the node list to include the now-placed fossil
	}
	if missing == false {
		IterateBMLengths(tree, iter)
	} else if missing == true {
		MissingTraitsEM(tree, iter)
	}
}

func randomNode(nodes []*Node) (rnode *Node) {
	s1 := rand.NewSource(time.Now().UnixNano())
	r1 := rand.New(s1)
	rnoden := r1.Intn(len(nodes))
	rnode = nodes[rnoden] //choose a random node
	return
}

//PruneTip will prune a fossil (or any) tip from a tree
//newpar is the parent of the tip you wish to prune, n is the other child of newpar
func PruneTip(newpar *Node, n *Node) {
	newpar.PAR.AddChild(n)
	newpar.RemoveChild(n)
	n.PAR = newpar.PAR
	newpar.PAR.RemoveChild(newpar)
	n.LEN = n.LEN + newpar.LEN
}

//this is the inverse of above
func GraftTip(newpar *Node, n *Node) {
	n.PAR.AddChild(newpar)
	n.PAR.RemoveChild(n)
	newpar.AddChild(n)
	newpar.LEN = n.LEN / 2.
	n.LEN = newpar.LEN
}
