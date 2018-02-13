package cophymaru

import (
	"math/rand"
)

//CalibrateSiteWeights will calculate weight vector to increase emphasis on sites concordant with the molecular reference tree during fossil placement procedures
func CalibrateSiteWeights(tree *Node, weightType string) (weights []float64) {
	siteLL := SitewiseLogLike(tree)
	for range siteLL {
		weights = append(weights, 0.)
	}
	for i := 0; i < 100; i++ {
		rtree := RandomUnrootedTree(tree)
		IterateBMLengths(rtree, 5)
		rsiteLL := SitewiseLogLike(rtree)
		for Si := range rsiteLL {
			if siteLL[Si] > rsiteLL[Si] {
				weights[Si] += 1. // add 1 to site weight if site i has a higher probability of occurring on the reference tree than on random tree j
			}
		}
	}
	if weightType == "float" {
		for i := range weights {
			weights[i] = weights[i] / 100.
		}
	}
	return
}

//RandomUnrootedTree will generate a random tree from all the taxa present in tree (the input tree will remain unaltered)
func RandomUnrootedTree(tree *Node) (root *Node) {
	nodes := tree.PreorderArray()
	var labels []string
	for _, n := range nodes { //make slice containing all of the tip labels
		if len(n.CHLD) == 0 {
			labels = append(labels, n.NAME)
		}
	}
	var randNodes []*Node
	for _, name := range labels {
		newtip := new(Node)
		newtip.NAME = name
		newtip.LEN = rand.Float64()
		for _, n := range nodes {
			if n.NAME == "" {
				continue
			} else if n.NAME == name {
				newtip.CONTRT = n.CONTRT
			}
		}
		MakeMissingDataSlice(newtip)
		newpar := new(Node)
		newpar.LEN = 0.1
		newpar.AddChild(newtip)
		for range newtip.CONTRT {
			newpar.CONTRT = append(newpar.CONTRT, float64(0.0))
			newpar.MIS = append(newpar.MIS, false)
		}
		randNodes = append(randNodes, newpar)
	}
	root = new(Node)
	for range randNodes[0].CONTRT {
		root.CONTRT = append(root.CONTRT, float64(0.0))
		root.MIS = append(root.MIS, false)
	}

	for i := 0; i < 3; i++ {
		rsub := RandomNode(randNodes)
		randNodes = popSampled(rsub, randNodes)
		tip := rsub.CHLD[0]
		rsub.RemoveChild(tip)
		tip.LEN = rand.Float64() * 2.
		root.AddChild(tip)
	}
	rootArray := root.PreorderArray()
	for _, node := range randNodes {
		reattach := RandomNode(rootArray[1:])
		node.LEN = rand.Float64()
		GraftFossilTip(node, reattach)
		rootArray = root.PreorderArray()
	}
	return
}

func popSampled(n *Node, nodes []*Node) (newNodes []*Node) {
	for _, node := range nodes {
		if node != n {
			newNodes = append(newNodes, node)
		}
	}
	return
}

func drawRandomLabel(n []string) (rtip string) {
	rnoden := rand.Intn(len(n))
	rtip = n[rnoden]
	return
}
