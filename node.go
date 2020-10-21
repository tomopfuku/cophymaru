package cophymaru

import (
	"bytes"
	"fmt"
	"os"
	"strconv"
)

//Node is a node struct to represent rooted and unrooted phylogenetic trees
type Node struct {
	Par    *Node
	Chs    []*Node
	Nam    string
	Len    float64
	PRNLen float64
	CONTRT []float64
	//VISITED bool
	Marked    bool
	MIS       []bool //this is a slice of bools that indicates whether the index is missing from CONTRT
	LL        []float64
	CONPRNLen []float64
	FAD       float64
	LAD       float64
	Height    float64
	FINDS     float64
	RATE      float64
	ISTIP     bool
	ANC       bool
	LSLen     float64
	LSPRNLen  float64
	OUTGRP    bool
	DIRDESC   bool
	DISCTRT   map[int]map[string]float64
	Num       int
}

func (n *Node) OldestDescendantAge() float64 {
	oldest := 0.0
	for _, c := range n.PreorderArray() {
		if c.ISTIP == false {
			continue
		}
		if c.FAD > oldest {
			oldest = c.FAD
		}
	}
	return oldest
}

func (n *Node) GetSib() *Node {
	if n.Nam == "root" {
		fmt.Println("Root node has no sibling")
		os.Exit(0)
	}
	par := n.Par
	var sib *Node
	if len(par.Chs) != 2 {
		if len(par.Chs) == 1 {
			fmt.Println("Singleton encountered in tree")
			os.Exit(0)
		} else {
			fmt.Println("Multifurcation found in tree")
			os.Exit(0)
		}
	}
	for _, c := range par.Chs {
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
	for _, nn := range n.Chs {
		if nn.ANC == true && nn.Nam+"_ancestral" == n.Nam {
			cont = nn
		}
	}
	return
}

//Root will root an unrooted tritomy tree in place
func (n *Node) Root(rootsublen float64) {
	if len(n.Chs) == 2 {
		fmt.Println("cannot root already rooted tree")
	}
	nn := &Node{n, nil, "", 0., 0.0, nil, false, nil, nil, nil, 0.0, 0.0, 0.0, 0.0, 0.0, false, false, 0.0, 0.0, false, false, nil, 0}
	nn.CONPRNLen = make([]float64, len(n.CONTRT))
	nn.CONTRT = make([]float64, len(n.CONTRT))
	nn.LL = make([]float64, len(n.CONTRT))
	for range nn.LL {
		nn.MIS = append(nn.MIS, false)
	}
	var ignodes, ognodes []*Node
	for _, c := range n.Chs {
		if c.ANC == true && c.ISTIP == true {
			nn.Nam = c.Nam + "_ancestral"
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
		//nn.Len = ognodes[0].Len / 2.0
		nn.Len = rootsublen
		ognodes[0].Len = ognodes[0].Len - rootsublen //nn.Len
		nn.LSLen = ognodes[0].LSLen / 2.0
		ognodes[0].LSLen = nn.LSLen
	} else if len(ognodes) == 2 {
		n.AddChild(nn)
		n.RemoveChild(ognodes[0])
		n.RemoveChild(ognodes[1])
		nn.AddChild(ognodes[0])
		nn.AddChild(ognodes[1])
		nn.Len = rootsublen
		ignodes[0].Len = ignodes[0].Len - rootsublen //nn.Len
		nn.LSLen = ignodes[0].LSLen / 2.0
		ignodes[0].LSLen = nn.LSLen

	} else {
		fmt.Println("there was a problem rooting the tree")
		os.Exit(0)
	}
	nn.Height = n.Height - nn.Len

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
	if len(n.Chs) != 2 {
		fmt.Println("cannot unroot already unrooted tree")
		os.Exit(0)
	}
	var igPar, ogPar *Node
	for _, c := range n.Chs {
		if c.OUTGRP == true {
			ogPar = c
		} else {
			igPar = c
		}
	}
	if len(ogPar.Chs) > 0 {
		rootlen = ogPar.Len
		igPar.Len += ogPar.Len
		igPar.LSLen += ogPar.LSLen
		for _, cc := range ogPar.Chs {
			ogPar.RemoveChild(cc)
			n.AddChild(cc)
			n.RemoveChild(ogPar)
		}
	} else {
		rootlen = igPar.Len
		ogPar.Len += igPar.Len
		ogPar.LSLen += igPar.LSLen
		for _, cc := range igPar.Chs {
			igPar.RemoveChild(cc)
			n.AddChild(cc)
			n.RemoveChild(igPar)
		}
	}
	return
}

func (n *Node) SetOutgroup(outgroup []string) {
	if n.Par != nil {
		fmt.Println("can't set outgroup from node that isn't the root")
	}
	tipdic := make(map[string]bool)
	for _, n := range n.PreorderTips() {
		tipdic[n.Nam] = true
	}
	outdic := make(map[string]bool)
	for _, tax := range outgroup {
		outdic[tax] = true
		if _, ok := tipdic[tax]; !ok {
			fmt.Println(tax)
			fmt.Println(tipdic)
			fmt.Println("one of the taxa named in the outgroup is not in the tree. please fix.")
			os.Exit(0)
		}
	}
	var ogPars []*Node
	for _, c := range n.Chs {
		for _, cc := range c.PreorderTips() {
			if _, ok := outdic[cc.Nam]; ok {
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
	for _, cn := range n.Chs {
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
	for _, cn := range n.Chs {
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
	for _, cn := range n.Chs {
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
	for in, cn := range n.Chs {
		if in == 0 {
			buffer.WriteString("(")
		}
		buffer.WriteString(cn.Newick(bl))
		if bl == true {
			s := strconv.FormatFloat(cn.Len, 'f', -1, 64)
			buffer.WriteString(":")
			buffer.WriteString(s)
		}
		if in == len(n.Chs)-1 {
			buffer.WriteString(")")
		} else {
			buffer.WriteString(",")
		}
	}
	buffer.WriteString(n.Nam)
	ret = buffer.String()
	return
}

//Phylogram will return a newick string representation of the tree with least squares estimates as the branch lengths
func (n Node) Phylogram() (ret string) {
	var buffer bytes.Buffer
	for in, cn := range n.Chs {
		if in == 0 {
			buffer.WriteString("(")
		}
		buffer.WriteString(cn.Phylogram())
		s := strconv.FormatFloat(cn.LSLen, 'f', -1, 64)
		buffer.WriteString(":")
		buffer.WriteString(s)
		if in == len(n.Chs)-1 {
			buffer.WriteString(")")
		} else {
			buffer.WriteString(",")
		}
	}
	buffer.WriteString(n.Nam)
	ret = buffer.String()
	return
}

//Rateogram will return a newick string representation of the tree with rates as the branch lengths
func (n Node) Rateogram() (ret string) {
	var buffer bytes.Buffer
	for in, cn := range n.Chs {
		if in == 0 {
			buffer.WriteString("(")
		}
		buffer.WriteString(cn.Rateogram())
		s := strconv.FormatFloat(cn.RATE, 'f', -1, 64)
		buffer.WriteString(":")
		buffer.WriteString(s)
		if in == len(n.Chs)-1 {
			buffer.WriteString(")")
		} else {
			buffer.WriteString(",")
		}
	}
	buffer.WriteString(n.Nam)
	ret = buffer.String()
	return
}

//CalcBranchRates plugs in branch-specific rates using the disparity length / time length
func (n *Node) CalcBranchRates() {
	for _, nn := range n.PreorderArray() {
		if nn == n {
			continue
		}
		if nn.LSLen != 0.0 && nn.Len != 0.0 {
			nn.RATE = nn.LSLen / nn.Len

		} else if nn.LSLen == 0.0 {
			nn.RATE = 0.0
		} else if nn.Len == 0.0 {
			nn.RATE = nn.LSLen / 0.1
		}
	}
}

//AddChild will add a child to a node
func (n *Node) AddChild(c *Node) error {
	for _, ch := range n.Chs {
		if c == ch {
			fmt.Println("You are trying to reroot on the current node!")
		}
	}
	n.Chs = append(n.Chs, c)
	c.Par = n
	return nil
}

//RemoveChild will remove a chidl from the slice of children associated with a node
func (n *Node) RemoveChild(c *Node) {
	var newChs []*Node
	for _, ch := range n.Chs {
		if ch != c {
			newChs = append(newChs, ch)
		}
		n.Chs = nil
		n.Chs = newChs
		c.Par = nil
	}
}

//NNodes is a helper method that will return the number of internal nodes descending from n (including n)
func (n *Node) NNodes(count *int) {
	*count++
	for _, ch := range n.Chs {
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
		curnode = curnode.Par
	}
	var newpar *Node
	for i := pathlen; i >= 1; i-- {
		newpar = pathnodes[i-1]
		curnode = pathnodes[i]
		curnode.RemoveChild(newpar)
		newpar.AddChild(curnode)
		curnode.LSLen = newpar.LSLen
	}
	//curnode = nil
	//newpar = nil
	n.LSLen = 0.0
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
		curnode = curnode.Par
	}
	var newpar *Node
	for i := pathlen; i >= 1; i-- {
		newpar = pathnodes[i-1]
		curnode = pathnodes[i]
		curnode.RemoveChild(newpar)
		newpar.AddChild(curnode)
		curnode.Len = newpar.Len
	}
	//curnode = nil
	//newpar = nil
	n.Len = 0.0
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
		curnode = curnode.Par
		curnode.Marked = false
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
		node.Marked = false
	}
}

//MarkAll will unmark all of the nodes on a tree
func (n *Node) MarkAll() {
	nodes := n.PreorderArray()
	for _, node := range nodes {
		node.Marked = true
	}
}
