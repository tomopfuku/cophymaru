package cophy

import (
	"fmt"
	"math"
	"os"
)

func postorder(curnode *Node) {
	for _, chld := range curnode.CHLD {
		postorder(chld)
	}
	fmt.Println(curnode.NAME, curnode.CONTRT)
}

//BMPruneRooted will prune BM branch lens and PICs down to a rooted node
//root node should be a real (ie. bifurcating) root
func BMPruneRooted(n *Node) {
	for _, chld := range n.CHLD {
		BMPruneRooted(chld)
	}
	n.PRNLEN = n.LEN
	nchld := len(n.CHLD)
	if nchld != 0 { //&& n.MRK == false {
		var temp_charst float64
		if nchld != 2 {
			fmt.Println("This BM pruning algorithm should only be perfomed on fully bifurcating trees/subtrees! Check for multifurcations and singletons.")
		}
		c0 := n.CHLD[0]
		c1 := n.CHLD[1]
		bot := ((1.0 / c0.PRNLEN) + (1.0 / c1.PRNLEN))
		n.PRNLEN += 1.0 / bot
		for i, _ := range n.CHLD[0].CONTRT {
			temp_charst = (((1 / c0.PRNLEN) * c0.CONTRT[i]) + ((1 / c1.PRNLEN) * c1.CONTRT[i])) / bot
			n.CONTRT[i] = temp_charst
		}
	}
}

func mark_children(n *Node) {
	n.MRK = true
	for _, ch := range n.CHLD {
		mark_children(ch)
	}
}

func debugParChld(tree *Node) {
	for _, ch := range tree.CHLD {
		fmt.Println(ch.NAME, "\t")
	}
}

//this is a quick check to make sure the tree passed is unrooted
func assertUnrootedTree(tree *Node) {
	if len(tree.CHLD) != 3 {
		fmt.Print("BRANCH LENGTHS MUST BE ITERATED ON AN UNROOTED TREE. THIS TREE IS ROOTED.")
		os.Exit(0)
	}
}

//IterateBMLengths will iteratively calculate the ML branch lengths for a particular topology
func IterateBMLengths(tree *Node, niter int) {
	assertUnrootedTree(tree)
	itercnt := 0
	for {
		calcBMLengths(tree)
		itercnt++
		if itercnt == niter {
			break
		}
	}
}

func calcBMLengths(tree *Node) {
	rnodes := tree.PreorderArray()
	lnode := 0
	for ind, newroot := range rnodes {
		if len(newroot.CHLD) == 0 {
			continue
		} else if newroot != rnodes[0] {
			tree = newroot.Reroot(rnodes[lnode])
			lnode = ind
		}
		for _, cn := range tree.CHLD {
			BMPruneRooted(cn)
		}
		TritomyML(tree)
	}
	tree = rnodes[0].Reroot(tree)
	//fmt.Println(tree.Newick(true))
}

//PruneToStar will prune brlens and traits to a root
func PruneToStar(tree *Node) {
	for _, cn := range tree.CHLD {
		BMPruneRooted(cn)
	}
}

/*
   some potentially sketchy math:
   var bot float64
   twopi := math.Pow((2.*math.Pi),float64(len(tree.CHLD[0].CONTRT)))
   fmt.Println("pi",twopi)
   bot = math.Pow(twopi,float64(len(tree.CHLD[0].CONTRT)))*math.Pow(((tree.CHLD[0].LEN*tree.CHLD[1].LEN)+(tree.CHLD[0].LEN*tree.CHLD[2].LEN)+(tree.CHLD[1].LEN*tree.CHLD[2].LEN)),float64(len(tree.CHLD[0].CONTRT))/2.0)
   bot = 1./bot
   D23 := float64(0.0)
   D13 := float64(0.0)
   D12 := float64(0.0)
   for i,_:= range tree.CHLD[0].CONTRT{
       D23 += math.Pow((tree.CHLD[1].CONTRT[i]-tree.CHLD[2].CONTRT[i]),2)
       D13 += math.Pow((tree.CHLD[0].CONTRT[i]-tree.CHLD[2].CONTRT[i]),2)
       D12 += math.Pow((tree.CHLD[1].CONTRT[i]-tree.CHLD[0].CONTRT[i]),2)
   }
   tright := (tree.CHLD[0].LEN*D23)+(tree.CHLD[1].LEN*D13)+(tree.CHLD[2].LEN*D12)
   bright := ((tree.CHLD[0].LEN*tree.CHLD[1].LEN)+(tree.CHLD[0].LEN*tree.CHLD[2].LEN)+(tree.CHLD[1].LEN*tree.CHLD[2].LEN))*2
   right := -(tright/bright)
   L := math.Log(bot)+right
*/
//fmt.Println("LL",L)
func CalcUnrootedLogLike(tree *Node) (chll float64) {
	chll = 0.0
	for _, ch := range tree.CHLD {
		curlike := 0.0
		CalcRootedLogLike(ch, &curlike)
		chll += curlike
	}
	cdiv := (float64(len(tree.CHLD[0].CONTRT)) / 2.)
	nconst := -(math.Log(2*math.Pi) * cdiv)
	v1 := tree.CHLD[0].PRNLEN
	v2 := tree.CHLD[1].PRNLEN
	v3 := tree.CHLD[2].PRNLEN
	v12 := v1 + v2
	var12 := (cdiv * math.Log(v12))
	csum := float64(0.0)
	lsum := float64(0.0)
	for i, _ := range tree.CHLD[0].CONTRT {
		csum += math.Pow((tree.CHLD[0].CONTRT[i]-tree.CHLD[1].CONTRT[i]), 2) / v12
		lsum += tree.CHLD[2].CONTRT[i] - ((tree.CHLD[1].LEN*tree.CHLD[0].CONTRT[i])+(tree.CHLD[0].LEN*tree.CHLD[1].CONTRT[i]))/(v12)
	}
	//csum = csum/v12
	csum = csum / 2.
	v1timesv2 := v1 * v2
	v32 := v3 + (v1timesv2 / v12)
	d := cdiv * math.Log(v32)
	//fmt.Println(lsum)
	lsum = math.Pow(lsum, 2)
	last := lsum / v32 / 2.
	lbot := nconst - var12 - csum - nconst - d - last
	//fmt.Println(lbot)
	chll += lbot
	return
}

func CalcRootedLogLike(n *Node, nlikes *float64) {
	for _, chld := range n.CHLD {
		CalcRootedLogLike(chld, nlikes)
	}
	n.PRNLEN = n.LEN
	nchld := len(n.CHLD)
	if nchld != 0 { //&& n.MRK == false {
		if nchld != 2 {
			fmt.Println("This BM pruning algorithm should only be perfomed on fully bifurcating trees/subtrees! Check for multifurcations and singletons.")
		}
		c0 := n.CHLD[0]
		c1 := n.CHLD[1]
		curlike := float64(0.0)
		var temp_charst float64
		temp_brlen := n.PRNLEN + ((c0.PRNLEN * c1.PRNLEN) / (c0.PRNLEN + c1.PRNLEN))
		for i, _ := range n.CHLD[0].CONTRT {
			cur_var := c0.PRNLEN + c1.PRNLEN
			contrast := c0.CONTRT[i] - c1.CONTRT[i]
			curlike += ((-0.5) * ((math.Log(2 * math.Pi)) + (math.Log(cur_var)) + (math.Pow(contrast, 2) / (cur_var))))
			temp_charst = ((c0.PRNLEN * c1.CONTRT[i]) + (c1.PRNLEN * c0.CONTRT[i])) / (cur_var)
			n.CONTRT[i] = temp_charst
		}
		*nlikes += curlike
		n.PRNLEN = temp_brlen
	}
}

func TritomyML(tree *Node) {
	ntraits := len(tree.CHLD[0].CONTRT)
	fntraits := float64(ntraits)
	var temp_v1, temp_v2, temp_v3, x1, x2, x3 float64
	temp_v1 = 0.0
	temp_v2 = 0.0
	temp_v3 = 0.0
	for i := range tree.CHLD[0].CONTRT {
		x1 = tree.CHLD[0].CONTRT[i]
		x2 = tree.CHLD[1].CONTRT[i]
		x3 = tree.CHLD[2].CONTRT[i]
		temp_v1 += ((x1 - x2) * (x1 - x3))
		temp_v2 += ((x2 - x1) * (x2 - x3))
		temp_v3 += ((x3 - x1) * (x3 - x2))
	}
	if temp_v1 < 0.0 {
		temp_v1 = 0.000001
		temp_v2 = 0.0
		temp_v3 = 0.0
		for i, _ := range tree.CHLD[0].CONTRT {
			x1 = tree.CHLD[0].CONTRT[i]
			x2 = tree.CHLD[1].CONTRT[i]
			x3 = tree.CHLD[2].CONTRT[i]
			temp_v2 += (x1 - x2) * (x1 - x2)
			temp_v3 += (x1 - x3) * (x1 - x3)
		}
	} else if temp_v2 < 0.0 {
		temp_v1 = 0.0
		temp_v2 = 0.00001
		temp_v3 = 0.0
		for i, _ := range tree.CHLD[0].CONTRT {
			x1 = tree.CHLD[0].CONTRT[i]
			x2 = tree.CHLD[1].CONTRT[i]
			x3 = tree.CHLD[2].CONTRT[i]
			temp_v1 += (x2 - x1) * (x2 - x1)
			temp_v3 += (x2 - x3) * (x2 - x3)
		}
	} else if temp_v3 < 0.0 {
		temp_v1 = 0.0
		temp_v2 = 0.0
		temp_v3 = 0.0001
		for i, _ := range tree.CHLD[0].CONTRT {
			x1 = tree.CHLD[0].CONTRT[i]
			x2 = tree.CHLD[1].CONTRT[i]
			x3 = tree.CHLD[2].CONTRT[i]
			temp_v1 += (x3 - x1) * (x3 - x1)
			temp_v2 += (x3 - x2) * (x3 - x2)
		}
	}
	temp_v1 = temp_v1 / fntraits
	temp_v2 = temp_v2 / fntraits
	temp_v3 = temp_v3 / fntraits
	temp_v1 = temp_v1 - (tree.CHLD[0].PRNLEN - tree.CHLD[0].LEN)
	temp_v2 = temp_v2 - (tree.CHLD[1].PRNLEN - tree.CHLD[1].LEN)
	temp_v3 = temp_v3 - (tree.CHLD[2].PRNLEN - tree.CHLD[2].LEN)
	if temp_v1 <= float64(0.0001) {
		temp_v1 = 0.0001
	} else if temp_v2 <= float64(0.0001) {
		temp_v2 = 0.0001
	} else if temp_v3 <= float64(0.0001) {
		temp_v3 = 0.0001
	}
	tree.CHLD[0].LEN = temp_v1
	tree.CHLD[1].LEN = temp_v2
	tree.CHLD[2].LEN = temp_v3
}
