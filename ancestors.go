package cophymaru

//Make2BudAncestor will change a bifurcating budding clade and switch to two independent budding events from node
func Make2BudAncestor(node *Node) {
	youngest := 100000000000000.0
	var moveChld *Node
	for _, c := range node.CHLD {
		if c.FAD < youngest {
			moveChld = c
			youngest = c.FAD
		}
	}
	oldest := moveChld.GetSib()
	oldest.LEN += node.LEN
	oldest.ANC = false
	node.ANC = false
	node.RemoveChild(oldest)
	moveTo := node.GetSib()
	node.PAR.RemoveChild(moveTo)
	node.PAR.AddChild(oldest)
	node.AddChild(moveTo)
	moveTo.LEN -= node.HEIGHT - moveTo.HEIGHT
	node.NAME = moveTo.NAME + "_ancestor1"
	node.ISTIP = true
	node.ANC = true
}

/*
func Make2BudTrees(ancestor, newdesc *Node) (bool, *Node) {
	par := ancestor.PAR
	if par.ANC == false {
		return true, nil
	}
	newanc := ancestor.DeepCopySingleNode()
	//dummy_anc := newanc.DeepCopySingleNode()
	newanc.NAME += "_ancestral"
	newanc.ISTIP = false
	newdesc.PAR.RemoveChild(newdesc)
	return false, nil
}
*/

func MakeAncestor(node *Node) bool {
	if node.ANC == true || node.ISTIP == false {
		return true
	}
	sib := node.GetSib()
	if node.FAD < sib.FAD {
		return true
	}
	sib.DIRDESC = true
	node.PAR.NAME = node.NAME + "_ancestral"
	node.PAR.MIS = node.MIS
	var newheight float64
	if node.FAD-node.LAD != 0.0 {
		//newheight := node.FAD - 0.05
		newheight = (sib.FAD + node.FAD) / 2.0
	} else {
		newheight = node.FAD
	}
	adjust := node.PAR.HEIGHT - newheight
	node.PAR.LEN += adjust
	node.LEN -= adjust
	sib.LEN -= adjust
	node.PAR.HEIGHT = newheight
	node.PAR.FAD = node.FAD
	node.RATE = 0.0
	if node.LEN < 0.0 {
		node.LEN = 0.0
	}
	node.ANC = true
	node.PAR.ANC = true
	return false
}

/*
func MakeAncestor(node *Node) bool {
	if node.ANC == true || node.ISTIP == false {
		return true
	}
	sib := node.GetSib()
	sib.DIRDESC = true
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
	fmt.Println(node.PAR.NAME, node.PAR.HEIGHT)

	return false
}
*/
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
	sib.DIRDESC = false
	node.PAR.NAME = ""
	return false
}
