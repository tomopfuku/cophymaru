package cophymaru

//Make2BudAncestor will change a bifurcating budding clade and switch to two independent budding events from node
func Make2BudAncestor(node *Node) {
	youngest := 100000000000000.0
	var moveChld *Node
	for _, c := range node.Chs {
		if c.FAD < youngest {
			moveChld = c
			youngest = c.FAD
		}
	}
	oldest := moveChld.GetSib()
	oldest.Len += node.Len
	oldest.ANC = false
	node.ANC = false
	node.RemoveChild(oldest)
	moveTo := node.GetSib()
	node.Par.RemoveChild(moveTo)
	node.Par.AddChild(oldest)
	node.AddChild(moveTo)
	moveTo.Len -= node.Height - moveTo.Height
	node.Nam = moveTo.Nam + "_ancestor1"
	node.ISTIP = true
	node.ANC = true
}

/*
func Make2BudTrees(ancestor, newdesc *Node) (bool, *Node) {
	par := ancestor.Par
	if par.ANC == false {
		return true, nil
	}
	newanc := ancestor.DeepCopySingleNode()
	//dummy_anc := newanc.DeepCopySingleNode()
	newanc.Nam += "_ancestral"
	newanc.ISTIP = false
	newdesc.Par.RemoveChild(newdesc)
	return false, nil
}
*/
func MakeAncestor(node *Node) bool {
	if node.ANC == true {
		return true
	}
	if node.ISTIP == false {
		return true
	}
	sib := node.GetSib()
	if node.FAD < sib.FAD {
		return true
	}
	if len(node.Chs) != 0 {
		return true
	}
	sib.DIRDESC = true
	node.Par.Nam = node.Nam + "_ancestral"
	node.Par.MIS = node.MIS
	var newheight float64
	if node.FAD-node.LAD != 0.0 {
		//newheight := node.FAD - 0.05
		newheight = (sib.FAD + node.FAD) / 2.0
	} else {
		newheight = node.FAD
	}
	adjust := node.Par.Height - newheight
	node.Par.Len += adjust
	node.Len -= adjust
	sib.Len -= adjust
	node.Par.Height = newheight
	node.Par.FAD = node.FAD
	node.RATE = 0.0
	if node.Len < 0.0 {
		node.Len = 0.0
	}
	node.ANC = true
	node.Par.ANC = true
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
	adj := node.Par.Height - sib.FAD       //constr
	node.Par.Len += adj                    //sib.FAD //sib.Len
	sib.Len -= adj                         //node.Par.Height - sib.FAD
	if sib.Len < 0.0 && sib.Len > -0.001 { //correct for any rounding errors
		sib.Len = 0.0
	}
	node.Par.Height = node.Par.Par.Height - node.Par.Len
	sib.Height = node.Par.Height - sib.Len
	if sib.ISTIP == false && sib.ANC == false {
		sib.FAD = sib.Height
	}
	node.RATE = 0.0
	node.Par.FAD = node.FAD
	node.Len -= adj
	node.Height = node.Par.Height - node.Len
	if node.Len < 0.0 {
		node.Len = 0.0
	}
	node.ANC = true
	node.Par.ANC = true
	node.Par.Nam = node.Nam + "_ancestral"
	fmt.Println(node.Par.Nam, node.Par.Height)

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
	sublen := (node.Par.FAD - node.Par.Height) //+ 0.01
	node.Par.Len -= sublen                     //node.Par.Par.Height - node.Par.Height
	node.Par.Height = node.Par.FAD             // + 0.01
	node.Par.FAD = node.Par.Height
	sib.Len += sublen //node.Par.Height - sib.FAD
	//sib.Height = node.Par.Height - sib.Len
	//if sib.ISTIP == false && sib.ANC == false {
	//	sib.FAD = sib.Height
	//}
	node.RATE = node.Par.RATE
	node.Len += sublen //+ 0.000001
	node.ANC = false
	node.Par.ANC = false
	sib.DIRDESC = false
	node.Par.Nam = ""
	return false
}
