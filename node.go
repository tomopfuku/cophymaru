package cophymaru

import (
	"bytes"
	"fmt"
	"os"
	"strconv"
)

//Node is a node struct to represent rooted and unrooted phylogenetic trees
type Node struct {
	PAR    *Node
	CHLD   []*Node
	NAME   string
	LEN    float64
	PRNLEN float64
	CONTRT []float64
	//VISITED bool
	MRK       bool
	MIS       []bool //this is a slice of bools that indicates whether the index is missing from CONTRT
	LL        []float64
	CONPRNLEN []float64
	FAD       float64
	LAD       float64
	HEIGHT    float64
	FINDS     float64
	RATE      float64
	ISTIP     bool
	ANC       bool
	LSLEN     float64
	LSPRNLEN  float64
	OUTGRP    bool
	DIRDESC   bool
	DISCTRT   map[int]map[string]float64
}

//DeepCopySingleNode will create a deep copy of a single node, without including its ancestors or descendants
func (n *Node) DeepCopySingleNode() *Node {
	nn := &Node{nil, nil, n.NAME + "_copy", n.LEN, n.PRNLEN, n.CONTRT, n.MRK, n.MIS, n.LL, n.CONPRNLEN, n.FAD, n.LAD, n.HEIGHT, n.FINDS, n.RATE, n.ISTIP, n.ANC, n.LSLEN, n.LSPRNLEN, n.OUTGRP, n.DIRDESC, n.DISCTRT}
	return nn
}

func (n *Node) GetSib() *Node {
	if n.NAME == "root" {
		fmt.Println("Root node has no sibling")
		os.Exit(0)
	}
	par := n.PAR
	var sib *Node
	if len(par.CHLD) != 2 {
		if len(par.CHLD) == 1 {
			fmt.Println("Singleton encountered in tree")
			os.Exit(0)
		} else {
			fmt.Println("Multifurcation found in tree")
			os.Exit(0)
		}
	}
	for _, c := range par.CHLD {
		if c != n {
			sib = c
		}
	}
	if sib == nil {
		fmt.Println("something is messed up with the tree. can't find a sister node for node", n)
		os.Exit(0)
	}
	return sib
}

//GetContDesc will find the descendent of a direct ancestor that represents the continuation of the lineage after a budding event
func (n *Node) GetContDesc() (cont *Node) {
	if n.ANC == false {
		fmt.Println("can't get the continuation of a hypothetical ancestor")
	}
	for _, nn := range n.CHLD {
		if nn.ANC == true && nn.NAME+"_ancestral" == n.NAME {
			cont = nn
		}
	}
	return
}

//Root will root an unrooted tritomy tree in place
func (n *Node) Root(rootsublen float64) {
	if len(n.CHLD) == 2 {
		fmt.Println("cannot root already rooted tree")
	}
	nn := &Node{n, nil, "", 0., 0.0, nil, false, nil, nil, nil, 0.0, 0.0, 0.0, 0.0, 0.0, false, false, 0.0, 0.0, false, false, nil}
	nn.CONPRNLEN = make([]float64, len(n.CONTRT))
	nn.CONTRT = make([]float64, len(n.CONTRT))
	nn.LL = make([]float64, len(n.CONTRT))
	for range nn.LL {
		nn.MIS = append(nn.MIS, false)
	}
	var ignodes, ognodes []*Node
	for _, c := range n.CHLD {
		if c.ANC == true {
			nn.NAME = c.NAME + "_ancestral"
			nn.ANC = true
			nn.FAD = c.FAD
		}
		if c.OUTGRP == true {
			ognodes = append(ognodes, c)
		} else {
			ignodes = append(ignodes, c)
		}
	}
	if len(ognodes) == 1 {
		n.AddChild(nn)
		n.RemoveChild(ignodes[0])
		n.RemoveChild(ignodes[1])
		nn.AddChild(ignodes[0])
		nn.AddChild(ignodes[1])
		//nn.LEN = ognodes[0].LEN / 2.0
		nn.LEN = rootsublen
		ognodes[0].LEN = ognodes[0].LEN - rootsublen //nn.LEN
		nn.LSLEN = ognodes[0].LSLEN / 2.0
		ognodes[0].LSLEN = nn.LSLEN
	} else if len(ognodes) == 2 {
		n.AddChild(nn)
		n.RemoveChild(ognodes[0])
		n.RemoveChild(ognodes[1])
		nn.AddChild(ognodes[0])
		nn.AddChild(ognodes[1])
		nn.LEN = rootsublen
		ignodes[0].LEN = ignodes[0].LEN - rootsublen //nn.LEN
		nn.LSLEN = ignodes[0].LSLEN / 2.0
		ignodes[0].LSLEN = nn.LSLEN

	} else {
		fmt.Println("there was a problem rooting the tree")
		os.Exit(0)
	}
	nn.HEIGHT = n.HEIGHT - nn.LEN

	if n.DISCTRT != nil {
		nn.LL = make([]float64, len(n.DISCTRT))
		nn.DISCTRT = make(map[int]map[string]float64)
		for site, probs := range n.DISCTRT {
			nn.DISCTRT[site] = make(map[string]float64)
			for state, prob := range probs {
				nn.DISCTRT[site][state] = prob
			}
		}
	}
}

func (n *Node) Unroot() (rootlen float64) {
	if len(n.CHLD) != 2 {
		fmt.Println("cannot unroot already unrooted tree")
		os.Exit(0)
	}
	var igPar, ogPar *Node
	for _, c := range n.CHLD {
		if c.OUTGRP == true {
			ogPar = c
		} else {
			igPar = c
		}
	}
	if len(ogPar.CHLD) > 0 {
		rootlen = ogPar.LEN
		igPar.LEN += ogPar.LEN
		igPar.LSLEN += ogPar.LSLEN
		for _, cc := range ogPar.CHLD {
			ogPar.RemoveChild(cc)
			n.AddChild(cc)
			n.RemoveChild(ogPar)
		}
	} else {
		rootlen = igPar.LEN
		ogPar.LEN += igPar.LEN
		ogPar.LSLEN += igPar.LSLEN
		for _, cc := range igPar.CHLD {
			igPar.RemoveChild(cc)
			n.AddChild(cc)
			n.RemoveChild(igPar)
		}
	}
	return
}

func (n *Node) SetOutgroup(outgroup []string) {
	if n.PAR != nil {
		fmt.Println("can't set outgroup from node that isn't the root")
	}
	tipdic := make(map[string]bool)
	for _, n := range n.PreorderTips() {
		tipdic[n.NAME] = true
	}
	outdic := make(map[string]bool)
	for _, tax := range outgroup {
		outdic[tax] = true
		if _, ok := tipdic[tax]; !ok {
			fmt.Println("one of the taxa named in the outgroup is not in the tree. please fix.")
			os.Exit(0)
		}
	}
	var ogPars []*Node
	for _, c := range n.CHLD {
		for _, cc := range c.PreorderTips() {
			if _, ok := outdic[cc.NAME]; ok {
				ogPars = append(ogPars, c)
			}
		}
	}
	for _, ogPar := range ogPars {
		ogPar.OUTGRP = true
		for _, n := range ogPar.PreorderArray() {
			n.OUTGRP = true
		}
	}
}

//PostorderArray will return an array of all the nodes in the tree in Postorder
func (n *Node) PostorderArray() (ret []*Node) {
	var buffer []*Node
	for _, cn := range n.CHLD {
		for _, cret := range cn.PreorderArray() {
			buffer = append(buffer, cret)
		}
	}
	buffer = append(buffer, n)
	ret = buffer
	return
}

//PreorderTips will return a preorder array of all the tips in a tree
func (n *Node) PreorderTips() (ret []*Node) {
	var buffer []*Node
	if n.ISTIP {
		buffer = append(buffer, n)
	}
	for _, cn := range n.CHLD {
		for _, cret := range cn.PreorderTips() {
			if cret.ISTIP {
				buffer = append(buffer, cret)
			}
		}
	}
	ret = buffer
	return
}

//PreorderArray will return a preordered array of all the nodes in a tree
func (n *Node) PreorderArray() (ret []*Node) {
	var buffer []*Node
	buffer = append(buffer, n)
	for _, cn := range n.CHLD {
		for _, cret := range cn.PreorderArray() {
			buffer = append(buffer, cret)
		}
	}
	ret = buffer
	return
}

//Newick will return a newick string representation of the tree
func (n Node) Newick(bl bool) (ret string) {
	var buffer bytes.Buffer
	for in, cn := range n.CHLD {
		if in == 0 {
			buffer.WriteString("(")
		}
		buffer.WriteString(cn.Newick(bl))
		if bl == true {
			s := strconv.FormatFloat(cn.LEN, 'f', -1, 64)
			buffer.WriteString(":")
			buffer.WriteString(s)
		}
		if in == len(n.CHLD)-1 {
			buffer.WriteString(")")
		} else {
			buffer.WriteString(",")
		}
	}
	buffer.WriteString(n.NAME)
	ret = buffer.String()
	return
}

//Phylogram will return a newick string representation of the tree with least squares estimates as the branch lengths
func (n Node) Phylogram() (ret string) {
	var buffer bytes.Buffer
	for in, cn := range n.CHLD {
		if in == 0 {
			buffer.WriteString("(")
		}
		buffer.WriteString(cn.Phylogram())
		s := strconv.FormatFloat(cn.LSLEN, 'f', -1, 64)
		buffer.WriteString(":")
		buffer.WriteString(s)
		if in == len(n.CHLD)-1 {
			buffer.WriteString(")")
		} else {
			buffer.WriteString(",")
		}
	}
	buffer.WriteString(n.NAME)
	ret = buffer.String()
	return
}

//Rateogram will return a newick string representation of the tree with rates as the branch lengths
func (n Node) Rateogram() (ret string) {
	var buffer bytes.Buffer
	for in, cn := range n.CHLD {
		if in == 0 {
			buffer.WriteString("(")
		}
		buffer.WriteString(cn.Rateogram())
		s := strconv.FormatFloat(cn.RATE, 'f', -1, 64)
		buffer.WriteString(":")
		buffer.WriteString(s)
		if in == len(n.CHLD)-1 {
			buffer.WriteString(")")
		} else {
			buffer.WriteString(",")
		}
	}
	buffer.WriteString(n.NAME)
	ret = buffer.String()
	return
}

//CalcBranchRates plugs in branch-specific rates using the disparity length / time length
func (n *Node) CalcBranchRates() {
	for _, nn := range n.PreorderArray() {
		if nn == n {
			continue
		}
		if nn.LSLEN != 0.0 && nn.LEN != 0.0 {
			nn.RATE = nn.LSLEN / nn.LEN

		} else if nn.LSLEN == 0.0 {
			nn.RATE = 0.0
		} else if nn.LEN == 0.0 {
			nn.RATE = nn.LSLEN / 0.1
		}
	}
}

//AddChild will add a child to a node
func (n *Node) AddChild(c *Node) error {
	for _, ch := range n.CHLD {
		if c == ch {
			fmt.Println("You are trying to reroot on the current node!")
		}
	}
	n.CHLD = append(n.CHLD, c)
	c.PAR = n
	return nil
}

//RemoveChild will remove a chidl from the slice of children associated with a node
func (n *Node) RemoveChild(c *Node) {
	var newCHLD []*Node
	for _, ch := range n.CHLD {
		if ch != c {
			newCHLD = append(newCHLD, ch)
		}
		n.CHLD = nil
		n.CHLD = newCHLD
		c.PAR = nil
	}
}

//NNodes is a helper method that will return the number of internal nodes descending from n (including n)
func (n *Node) NNodes(count *int) {
	*count++
	for _, ch := range n.CHLD {
		ch.NNodes(count)
	}
}

//RerootLS reroots all the nodes represented in a graph on n
func (n *Node) RerootLS(oldroot *Node) *Node {
	if n == oldroot {
		fmt.Println("you are trying to reroot on the current root!")
	}
	nnodes := 0
	oldroot.NNodes(&nnodes)
	var pathnodes = make([]*Node, nnodes)
	//var pathnodes []*Node
	curnode := n
	pathlen := 0 //this will count the number of nodes between the newroot and the oldroot
	for ind := range pathnodes {
		pathnodes[ind] = curnode
		//pathnodes = append(pathnodes,curnode)
		if curnode == oldroot {
			break
		}
		pathlen++
		curnode = curnode.PAR
	}
	var newpar *Node
	for i := pathlen; i >= 1; i-- {
		newpar = pathnodes[i-1]
		curnode = pathnodes[i]
		curnode.RemoveChild(newpar)
		newpar.AddChild(curnode)
		curnode.LSLEN = newpar.LSLEN
	}
	//curnode = nil
	//newpar = nil
	n.LSLEN = 0.0
	return n
}

//Reroot reroots all the nodes represented in a graph on n
func (n *Node) Reroot(oldroot *Node) *Node {
	if n == oldroot {
		fmt.Println("you are trying to reroot on the current root!")
	}
	nnodes := 0
	oldroot.NNodes(&nnodes)
	var pathnodes = make([]*Node, nnodes)
	//var pathnodes []*Node
	curnode := n
	pathlen := 0 //this will count the number of nodes between the newroot and the oldroot
	for ind := range pathnodes {
		pathnodes[ind] = curnode
		//pathnodes = append(pathnodes,curnode)
		if curnode == oldroot {
			break
		}
		pathlen++
		curnode = curnode.PAR
	}
	var newpar *Node
	for i := pathlen; i >= 1; i-- {
		newpar = pathnodes[i-1]
		curnode = pathnodes[i]
		curnode.RemoveChild(newpar)
		newpar.AddChild(curnode)
		curnode.LEN = newpar.LEN
	}
	//curnode = nil
	//newpar = nil
	n.LEN = 0.0
	return n
}

//UnmarkToRoot will mark all of the nodes betwen the current node and the root to be recalculated
func (n *Node) UnmarkToRoot(oldroot *Node) {
	if n == oldroot {
		fmt.Println("you are trying to reroot on the current root!")
	}
	curnode := n
	pathlen := int(0)
	for {
		pathlen++
		curnode = curnode.PAR
		curnode.MRK = false
		if curnode == oldroot {
			break
		}
		if pathlen > 2000000 {
			fmt.Println("something is wrong with the tree. could not walk back to root.")
			os.Exit(0)
		}
	}
}

//UnmarkAll will unmark all of the nodes on a tree
func (n *Node) UnmarkAll() {
	nodes := n.PreorderArray()
	for _, node := range nodes {
		node.MRK = false
	}
}

//MarkAll will unmark all of the nodes on a tree
func (n *Node) MarkAll() {
	nodes := n.PreorderArray()
	for _, node := range nodes {
		node.MRK = true
	}
}
