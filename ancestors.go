package cophymaru

func MakeAncestorSameRate(node *Node) bool {
	if node.ANC == true || node.ISTIP == false {
		return true
	}
	sib := node.GetSib()
	if node.FAD < sib.FAD && sib.ISTIP == true {
		return true
	}
	lensum := node.LEN + node.PAR.LEN
	nprop := node.LEN / lensum
	descprop := node.PAR.LEN / lensum
	node.PAR.RATE = (node.RATE * nprop) + (node.PAR.RATE * descprop)

	adj := node.PAR.HEIGHT - sib.FAD       //constr
	node.PAR.LEN += adj                    //sib.FAD //sib.LEN
	sib.LEN -= adj                         //node.PAR.HEIGHT - sib.FAD
	if sib.LEN < 0.0 && sib.LEN > -0.001 { //correct for any rounding errors
		sib.LEN = 0.0
	}
	node.PAR.HEIGHT = node.PAR.PAR.HEIGHT - node.PAR.LEN
	sib.HEIGHT = node.PAR.HEIGHT - sib.LEN
	if sib.ISTIP == false && sib.ANC == false {
		sib.FAD = sib.HEIGHT
	}

	//sib.RATE += node.RATE
	node.RATE = 0.0
	node.PAR.FAD = node.FAD
	node.LEN -= adj
	node.HEIGHT = node.PAR.HEIGHT - node.LEN
	if node.LEN < 0.0 {
		node.LEN = 0.0
	}
	node.ANC = true
	node.PAR.ANC = true
	node.PAR.NAME = node.NAME + "_ancestral"
	return false
}

func MakeAncestor(node *Node) bool {
	if node.ANC == true || node.ISTIP == false {
		return true
	}
	sib := node.GetSib()
	if node.FAD < sib.FAD && sib.ISTIP == true {
		return true
	}
	adj := node.PAR.HEIGHT - sib.FAD       //constr
	node.PAR.LEN += adj                    //sib.FAD //sib.LEN
	sib.LEN -= adj                         //node.PAR.HEIGHT - sib.FAD
	if sib.LEN < 0.0 && sib.LEN > -0.001 { //correct for any rounding errors
		sib.LEN = 0.0
	}
	node.PAR.HEIGHT = node.PAR.PAR.HEIGHT - node.PAR.LEN
	sib.HEIGHT = node.PAR.HEIGHT - sib.LEN
	if sib.ISTIP == false && sib.ANC == false {
		sib.FAD = sib.HEIGHT
	}
	node.RATE = 0.0
	node.PAR.FAD = node.FAD
	node.LEN -= adj
	node.HEIGHT = node.PAR.HEIGHT - node.LEN
	if node.LEN < 0.0 {
		node.LEN = 0.0
	}
	node.ANC = true
	node.PAR.ANC = true
	node.PAR.NAME = node.NAME + "_ancestral"
	return false
}

func MakeAncestorLabel(label string, nodes []*Node) bool {
	n := NodeFromLabel(label, nodes)
	bad := MakeAncestor(n)
	return bad
}

func UnmakeAncestorLabel(label string, nodes []*Node) {
	n := NodeFromLabel(label, nodes)
	UnmakeAncestor(n)
}

func UnmakeAncestor(node *Node) bool {
	if node.ANC != true || node.ISTIP == false {
		return true
	}
	sib := node.GetSib()
	sublen := (node.PAR.FAD - node.PAR.HEIGHT) //+ 0.01
	node.PAR.LEN -= sublen                     //node.PAR.PAR.HEIGHT - node.PAR.HEIGHT
	node.PAR.HEIGHT = node.PAR.FAD             // + 0.01
	node.PAR.FAD = node.PAR.HEIGHT
	sib.LEN += sublen //node.PAR.HEIGHT - sib.FAD
	//sib.HEIGHT = node.PAR.HEIGHT - sib.LEN
	//if sib.ISTIP == false && sib.ANC == false {
	//	sib.FAD = sib.HEIGHT
	//}
	node.RATE = node.PAR.RATE
	node.LEN += sublen //+ 0.000001
	node.ANC = false
	node.PAR.ANC = false
	node.PAR.NAME = ""
	return false
}
