package cophymaru

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
			temp_charst = (((1 / c0.PRNLEN) * c1.CONTRT[i]) + ((1 / c1.PRNLEN) * c0.CONTRT[i])) / bot
			n.CONTRT[i] = temp_charst
		}
	}
}

//BMPruneRootedSingle will prune BM branch lens and calculate PIC of a single trait down to a rooted node
//root node should be a real (ie. bifurcating) root
func BMPruneRootedSingle(n *Node, i int) {
	for _, chld := range n.CHLD {
		BMPruneRootedSingle(chld, i)
	}
	n.PRNLEN = n.LEN
	nchld := len(n.CHLD)
	if nchld != 0 { //&& n.MRK == false {
		if nchld != 2 {
			fmt.Println("This BM pruning algorithm should only be perfomed on fully bifurcating trees/subtrees! Check for multifurcations and singletons.")
		}
		c0 := n.CHLD[0]
		c1 := n.CHLD[1]
		bot := ((1.0 / c0.PRNLEN) + (1.0 / c1.PRNLEN))
		n.PRNLEN += 1.0 / bot
		tempCharacter := (((1 / c0.PRNLEN) * c1.CONTRT[i]) + ((1 / c1.PRNLEN) * c0.CONTRT[i])) / bot
		n.CONTRT[i] = tempCharacter
	}
}

func markChildren(n *Node) {
	n.MRK = true
	for _, ch := range n.CHLD {
		markChildren(ch)
	}
}

func debugParChld(tree *Node) {
	for _, ch := range tree.CHLD {
		fmt.Println(ch.NAME, "\t")
	}
}

//AssertUnrootedTree is a quick check to make sure the tree passed is unrooted
func AssertUnrootedTree(tree *Node) {
	if len(tree.CHLD) != 3 {
		fmt.Print("BRANCH LENGTHS MUST BE ITERATED ON AN UNROOTED TREE. THIS TREE IS ROOTED.")
		os.Exit(0)
	}
}

//IterateBMLengths will iteratively calculate the ML branch lengths for a particular topology, assuming that traits are fully sampled at the tips. should use MissingTraitsEM if there are missing sites.
func IterateBMLengths(tree *Node, niter int) {
	AssertUnrootedTree(tree)
	itercnt := 0
	for {
		calcBMLengths(tree)
		itercnt++
		if itercnt == niter {
			break
		}
	}
}

//MissingTraitsEM will iteratively calculate the ML branch lengths for a particular topology
func MissingTraitsEM(tree *Node, niter int) {
	AssertUnrootedTree(tree)
	nodes := tree.PreorderArray()
	InitMissingValues(nodes)
	itercnt := 0
	for {
		CalcExpectedTraits(tree) //calculate Expected trait values
		calcBMLengths(tree)      //maximize likelihood of branch lengths
		itercnt++
		if itercnt == niter {
			break
		}
	}
}

//calcBMLengths will perform a single pass of the branch length ML estimation
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

//CalcUnrootedLogLike will calculate the log-likelihood of an unrooted tree, while assuming that some sites have missing data. This _can_ be used to calculate the likelihoods of trees that have complete trait sampling, but it will be slower than CalcRootedLogLike.
func CalcUnrootedLogLike(tree *Node) (chll float64) {
	chll = 0.0
	for _, ch := range tree.CHLD {
		curlike := 0.0
		CalcRootedLogLike(ch, &curlike)
		chll += curlike
	}
	sitelikes := 0.0
	var tmpll float64
	var contrast, curVar float64
	for i := range tree.CHLD[0].CONTRT {
		tmpll = 0.
		if tree.CHLD[0].MIS[i] == false && tree.CHLD[1].MIS[i] == false && tree.CHLD[2].MIS[i] == false { //do the standard calculation when no subtrees have missing traits
			contrast = tree.CHLD[0].CONTRT[i] - tree.CHLD[1].CONTRT[i]
			curVar = tree.CHLD[0].LEN
			tmpll = ((-0.5) * ((math.Log(2. * math.Pi)) + (math.Log(curVar)) + (math.Pow(contrast, 2.) / (curVar))))
			tmpPRNLEN := ((tree.CHLD[0].PRNLEN * tree.CHLD[1].PRNLEN) / (tree.CHLD[0].PRNLEN + tree.CHLD[1].PRNLEN))
			tmpChar := ((tree.CHLD[0].PRNLEN * tree.CHLD[1].CONTRT[i]) + (tree.CHLD[1].PRNLEN * tree.CHLD[0].CONTRT[i])) / curVar
			contrast = tmpChar - tree.CHLD[2].CONTRT[i]
			curVar = tree.CHLD[2].PRNLEN + tmpPRNLEN
			tmpll += ((-0.5) * ((math.Log(2. * math.Pi)) + (math.Log(curVar)) + (math.Pow(contrast, 2.) / (curVar))))
		}
		sitelikes += tmpll
	}
	chll += sitelikes
	return
}

/*/CalcUnrootedLogLike will calculate the log-likelihood of an unrooted tree assuming that traits are completely sampled
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
	for i := range tree.CHLD[0].CONTRT {
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
	chll += lbot
	return
}
*/

//CalcRootedLogLike will return the BM likelihood of a tree assuming that no data are missing from the tips.
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
		var tempChar float64
		tempBranchLength := n.PRNLEN + ((c0.PRNLEN * c1.PRNLEN) / (c0.PRNLEN + c1.PRNLEN))
		for i := range n.CHLD[0].CONTRT {
			curVar := c0.PRNLEN + c1.PRNLEN
			contrast := c0.CONTRT[i] - c1.CONTRT[i]
			curlike += ((-0.5) * ((math.Log(2 * math.Pi)) + (math.Log(curVar)) + (math.Pow(contrast, 2) / (curVar))))
			tempChar = ((c0.PRNLEN * c1.CONTRT[i]) + (c1.PRNLEN * c0.CONTRT[i])) / (curVar)
			n.CONTRT[i] = tempChar
		}
		*nlikes += curlike
		n.PRNLEN = tempBranchLength
	}
}

//MissingUnrootedLogLike will calculate the log-likelihood of an unrooted tree, while assuming that some sites have missing data. This _can_ be used to calculate the likelihoods of trees that have complete trait sampling, but it will be slower than CalcRootedLogLike.
func MissingUnrootedLogLike(tree *Node) (chll float64) {
	chll = 0.0
	for _, ch := range tree.CHLD {
		chll += MissingRootedLogLike(ch)
	}
	sitelikes := 0.0
	var tmpll float64
	var contrast, curVar float64
	for i := range tree.CHLD[0].CONTRT {
		tmpll = 0.
		if tree.CHLD[0].MIS[i] == false && tree.CHLD[1].MIS[i] == false && tree.CHLD[2].MIS[i] == false { //do the standard calculation when no subtrees have missing traits
			contrast = tree.CHLD[0].CONTRT[i] - tree.CHLD[1].CONTRT[i]
			curVar = tree.CHLD[0].LEN
			tmpll = ((-0.5) * ((math.Log(2. * math.Pi)) + (math.Log(curVar)) + (math.Pow(contrast, 2.) / (curVar))))
			tmpPRNLEN := ((tree.CHLD[0].PRNLEN * tree.CHLD[1].PRNLEN) / (tree.CHLD[0].PRNLEN + tree.CHLD[1].PRNLEN))
			tmpChar := ((tree.CHLD[0].PRNLEN * tree.CHLD[1].CONTRT[i]) + (tree.CHLD[1].PRNLEN * tree.CHLD[0].CONTRT[i])) / curVar
			contrast = tmpChar - tree.CHLD[2].CONTRT[i]
			curVar = tree.CHLD[2].PRNLEN + tmpPRNLEN
			tmpll += ((-0.5) * ((math.Log(2. * math.Pi)) + (math.Log(curVar)) + (math.Pow(contrast, 2.) / (curVar))))
		} else if tree.CHLD[0].MIS[i] == false && tree.CHLD[1].MIS[i] == false && tree.CHLD[2].MIS[i] == true { // do standard "rooted" calculation on CHLD[0] and CHLD [1] if CHLD[2] is missing
			contrast = tree.CHLD[0].CONTRT[i] - tree.CHLD[1].CONTRT[i]
			curVar = tree.CHLD[0].PRNLEN + tree.CHLD[1].PRNLEN
			tmpll = ((-0.5) * ((math.Log(2. * math.Pi)) + (math.Log(curVar)) + (math.Pow(contrast, 2.) / (curVar))))
		} else if tree.CHLD[0].MIS[i] == false && tree.CHLD[2].MIS[i] == false && tree.CHLD[1].MIS[i] == true { // do standard "rooted" calculation on CHLD[0] and CHLD [2] if CHLD[1] is missing
			contrast = tree.CHLD[0].CONTRT[i] - tree.CHLD[2].CONTRT[i]
			curVar = tree.CHLD[0].PRNLEN + tree.CHLD[2].PRNLEN
			tmpll = ((-0.5) * ((math.Log(2. * math.Pi)) + (math.Log(curVar)) + (math.Pow(contrast, 2.) / (curVar))))
		} else if tree.CHLD[1].MIS[i] == false && tree.CHLD[2].MIS[i] == false && tree.CHLD[0].MIS[i] == true { // do standard "rooted" calculation on CHLD[1] and CHLD [2] if CHLD[0] is missing
			contrast = tree.CHLD[1].CONTRT[i] - tree.CHLD[2].CONTRT[i]
			curVar = tree.CHLD[1].PRNLEN + tree.CHLD[2].PRNLEN
			tmpll = ((-0.5) * ((math.Log(2. * math.Pi)) + (math.Log(curVar)) + (math.Pow(contrast, 2.) / (curVar))))
		}
		sitelikes += tmpll
	}
	chll += sitelikes
	return
}

//MissingRootedLogLike will return the BM log-likelihood of a tree for a single site, pruning tips that have missing data
func MissingRootedLogLike(n *Node) (sitelikes float64) {
	nodes := n.PostorderArray()
	var contrast, curVar float64
	sitelikes = 0.
	nodelikes := 0.
	for site := range nodes[0].CONTRT {
		nodelikes = 0.
		for _, node := range nodes {
			node.PRNLEN = node.LEN
			if len(node.CHLD) == 0 {
				continue
			}
			if node.CHLD[0].MIS[site] == false && node.CHLD[1].MIS[site] == false {
				contrast = node.CHLD[0].CONTRT[site] - node.CHLD[1].CONTRT[site]
				curVar = node.CHLD[0].PRNLEN + node.CHLD[1].PRNLEN
				nodelikes += ((-0.5) * ((math.Log(2. * math.Pi)) + (math.Log(curVar)) + (math.Pow(contrast, 2.) / (curVar))))
				node.CONTRT[site] = ((node.CHLD[0].PRNLEN * node.CHLD[1].CONTRT[site]) + (node.CHLD[1].PRNLEN * node.CHLD[0].CONTRT[site])) / curVar
				node.PRNLEN += ((node.CHLD[0].PRNLEN * node.CHLD[1].PRNLEN) / (node.CHLD[0].PRNLEN + node.CHLD[1].PRNLEN))
			} else if node.CHLD[0].MIS[site] == false && node.CHLD[1].MIS[site] == true {
				node.PRNLEN += node.CHLD[0].PRNLEN
				node.CONTRT[site] = node.CHLD[0].CONTRT[site]
			} else if node.CHLD[1].MIS[site] == false && node.CHLD[0].MIS[site] == true {
				node.PRNLEN += node.CHLD[1].PRNLEN
				node.CONTRT[site] = node.CHLD[1].CONTRT[site]
			} else if node.CHLD[0].MIS[site] == true && node.CHLD[1].MIS[site] == true {
				node.MIS[site] = true
			}
		}
		sitelikes += nodelikes
	}
	return
}

//TritomyML will calculate the MLEs for the branch lengths of a tifurcating 3-taxon tree
func TritomyML(tree *Node) {
	ntraits := len(tree.CHLD[0].CONTRT)
	fntraits := float64(ntraits)
	var x1, x2, x3 float64
	sumV1 := 0.0
	sumV2 := 0.0
	sumV3 := 0.0
	for i := range tree.CHLD[0].CONTRT {
		x1 = tree.CHLD[0].CONTRT[i]
		x2 = tree.CHLD[1].CONTRT[i]
		x3 = tree.CHLD[2].CONTRT[i]
		sumV1 += ((x1 - x2) * (x1 - x3))
		sumV2 += ((x2 - x1) * (x2 - x3))
		sumV3 += ((x3 - x1) * (x3 - x2))
	}
	if sumV1 < 0.0 {
		sumV1 = 0.000001
		sumV2 = 0.0
		sumV3 = 0.0
		for i := range tree.CHLD[0].CONTRT {
			x1 = tree.CHLD[0].CONTRT[i]
			x2 = tree.CHLD[1].CONTRT[i]
			x3 = tree.CHLD[2].CONTRT[i]
			sumV2 += (x1 - x2) * (x1 - x2)
			sumV3 += (x1 - x3) * (x1 - x3)
		}
	} else if sumV2 < 0.0 {
		sumV1 = 0.0
		sumV2 = 0.00001
		sumV3 = 0.0
		for i := range tree.CHLD[0].CONTRT {
			x1 = tree.CHLD[0].CONTRT[i]
			x2 = tree.CHLD[1].CONTRT[i]
			x3 = tree.CHLD[2].CONTRT[i]
			sumV1 += (x2 - x1) * (x2 - x1)
			sumV3 += (x2 - x3) * (x2 - x3)
		}
	} else if sumV3 < 0.0 {
		sumV1 = 0.0
		sumV2 = 0.0
		sumV3 = 0.0001
		for i := range tree.CHLD[0].CONTRT {
			x1 = tree.CHLD[0].CONTRT[i]
			x2 = tree.CHLD[1].CONTRT[i]
			x3 = tree.CHLD[2].CONTRT[i]
			sumV1 += (x3 - x1) * (x3 - x1)
			sumV2 += (x3 - x2) * (x3 - x2)
		}
	}
	sumV1 = sumV1 / fntraits
	sumV2 = sumV2 / fntraits
	sumV3 = sumV3 / fntraits
	sumV1 = sumV1 - (tree.CHLD[0].PRNLEN - tree.CHLD[0].LEN)
	sumV2 = sumV2 - (tree.CHLD[1].PRNLEN - tree.CHLD[1].LEN)
	sumV3 = sumV3 - (tree.CHLD[2].PRNLEN - tree.CHLD[2].LEN)
	if sumV1 < 0. {
		sumV1 = 0.0001
	}
	if sumV2 < 0. {
		sumV2 = 0.0001
	}
	if sumV3 < 0. {
		sumV3 = 0.0001
	}
	tree.CHLD[0].LEN = sumV1
	tree.CHLD[1].LEN = sumV2
	tree.CHLD[2].LEN = sumV3
}
