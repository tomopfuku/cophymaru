package cophymaru

func MakeAncestor(node *Node) bool {
	sib := node.GetSib()
	if node.FAD < sib.FAD {
		return true
	}
	node.PAR.LEN += node.PAR.HEIGHT - sib.FAD //sib.LEN
	sib.LEN -= node.PAR.HEIGHT - sib.FAD
	node.PAR.HEIGHT = node.PAR.PAR.HEIGHT - node.PAR.LEN
	sib.HEIGHT = node.PAR.HEIGHT - sib.LEN
	sib.FAD = sib.HEIGHT
	node.RATE = 0.0
	node.PAR.FAD = node.FAD
	node.LEN = node.PAR.HEIGHT - node.LAD //0.0
	if node.LEN < 0.0 {
		node.LEN = 0.0
	}
	node.ANC = true
	node.PAR.ANC = true
	node.PAR.NAME = node.NAME + "_ancestral"
	return false
}

/*
func MakeAncestor(node *Node) {
	node.RATE = 0.0
	sib := node.GetSib()
	sib.LEN += node.PAR.LEN
	node.LEN += node.PAR.LEN
	node.PAR.HEIGHT += node.PAR.LEN
	node.PAR.LEN = 0.0
	node.ANC = true
	//node.PAR.LEN += node.LEN
	//node.PAR.LEN -= (sib.HEIGHT - node.HEIGHT)
	//node.LEN = 0.0
	//sib.ANC = true
}
*/

func MakeAncestorLabel(label string, nodes []*Node) {
	n := NodeFromLabel(label, nodes)
	MakeAncestor(n)
}
