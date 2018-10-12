package cophymaru

import "math/rand"

//RandomStartingTree will generate a random tree from all the taxa present in tree (the input tree will remain unaltered)
func RandomStartingTree(traits map[string][]float64) (root *Node) {
	var randNodes []*Node
	for name, tr := range traits {
		newtip := new(Node)
		newtip.NAME = name
		newtip.LEN = rand.Float64()
		newtip.CONTRT = tr
		MakeMissingDataSlice(newtip)
		newpar := new(Node)
		newpar.LEN = 0.1
		newpar.AddChild(newtip)
		for range newtip.CONTRT {
			newpar.CONTRT = append(newpar.CONTRT, float64(0.0))
			newpar.MIS = append(newpar.MIS, false)
			newpar.LL = append(newpar.LL, 0.)
		}
		randNodes = append(randNodes, newpar)
	}
	root = new(Node)
	for range randNodes[0].CONTRT {
		root.CONTRT = append(root.CONTRT, float64(0.0))
		root.MIS = append(root.MIS, false)
		root.LL = append(root.LL, 0.)
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
