package cophy

import (
	"bufio"
	"fmt"
	"log"
	"math"
	"math/rand"
	"os"
	"strconv"
	"time"
)

type treeMCMC struct {
	NGEN int
}

func brlenSlidingWindow(theta float64, wsize float64) (thetaStar float64) {
	u := rand.Float64()
	thetaStar = theta - (wsize / 2.) + (wsize * u)
	if thetaStar < 0. {
		thetaStar = -thetaStar
	}
	return
}

func singleBrlenMultiplierProp(theta float64) (thetaStar, propRat float64) {
	min := math.Log(0.001)
	u := rand.Float64()
	epsilon := 0.5
	c := math.Exp(((u - 0.5) * epsilon))
	thetaStar = theta * c
	propRat = c
	if math.Log(thetaStar) < min { //place a lower constraint on brlen
		thetaStar = math.Exp(2*min - math.Log(thetaStar))
		propRat = thetaStar / theta
	}
	return
}

func getProposedBrlens(nodes []*Node) []float64 {
	var oldL []float64
	for _, i := range nodes {
		oldL = append(oldL, i.LEN)
		i.LEN = brlenSlidingWindow(i.LEN, 0.1)
	}
	return oldL
}

func getFossilNodesFromLabel(fnames []string, nodes []*Node) []*Node {
	var ret []*Node
	for _, n := range nodes {
		for _, lab := range fnames {
			if n.NAME == lab {
				ret = append(ret, n)
			}
		}
	}
	return ret
}

func writeTreeFile(newick string, outFl *bufio.Writer) {
	fmt.Fprint(outFl, newick+"\n")
}

func initOutfile(flnm string) *bufio.Writer {
	f, err := os.Create(flnm)
	if err != nil {
		log.Fatal(err)
	}
	defer f.Close()
	w := bufio.NewWriter(f)
	return w
}

//MCMC will run Markov Chain Monte Carlo simulations, adjusting branch lengths and fossil placements
func MCMC(tree *Node, gen int, fnames []string, treeOutFile, logOutFile string, branchPrior string) {
	f, err := os.Create(treeOutFile)
	if err != nil {
		log.Fatal(err)
	}
	w := bufio.NewWriter(f)
	lf, err := os.Create(logOutFile)
	if err != nil {
		log.Fatal(err)
	}
	lw := bufio.NewWriter(lf)

	nodes := tree.PreorderArray()[1:]
	fos := getFossilNodesFromLabel(fnames, nodes)
	var inNodes []*Node
	for _, n := range nodes {
		if len(n.CHLD) == 2 {
			inNodes = append(inNodes, n)
		}
	}
	ll := CalcUnrootedLogLike(tree)
	var lp float64
	if branchPrior == "0" {
		lp = math.Log(1.) //ExponentialBranchLengthLogPrior(nodes,10.)
	} else if branchPrior == "1" {
		lp = ExponentialBranchLengthLogPrior(nodes, 10.)
	}
	for i := 0; i < gen; i++ {
		if i%5 == 0 {
			lp, ll = fossilPlacementUpdate(ll, lp, fos, nodes, tree)
		}
		if i%5 != 0 {
			lp, ll = singleBranchLengthUpdate(ll, lp, nodes, inNodes, tree, branchPrior) //NOTE uncomment to sample BRLENS
		}
		if i%10 == 0 {
			fmt.Fprint(lw, strconv.Itoa(i)+"\t"+strconv.FormatFloat(lp, 'f', -1, 64)+"\t"+strconv.FormatFloat(ll, 'f', -1, 64)+"\t")
			for _, ln := range nodes {
				fmt.Fprint(lw, strconv.FormatFloat(ln.LEN, 'f', -1, 64)+"\t")
			}
			fmt.Fprint(lw, "\n")
			//writeTreeFile(tree.Newick(true),w)
			fmt.Fprint(w, tree.Newick(true)+"\n")
		}
	}
	lw.Flush()
	err = w.Flush()
	if err != nil {
		log.Fatal(err)
	}
	f.Close()
}

//this move prunes and regrafts a fossil, creating random variables for the altered branch lengths
//the procedure is basically the same as the SPR move as described in the Yang (2014) Mol. Evol. book
func fossilPlacementUpdate(ll, lp float64, fnodes, nodes []*Node, tree *Node) (float64, float64) {
	fn := drawRandomNode(fnodes) //draw a random fossil tip
	reattach := drawRandomReattachment(fn, nodes)
	x, p, lastn := pruneFossilTip(fn)
	r := graftFossilTip(fn.PAR, reattach)
	propRat := r / (x * p)
	llstar := CalcUnrootedLogLike(tree)
	lpstar := math.Log(1.)
	alpha := math.Exp(lpstar-lp) * math.Exp(llstar-ll) * propRat
	s1 := rand.NewSource(time.Now().UnixNano())
	r1 := rand.New(s1)
	r = r1.Float64()
	if r < alpha {
		ll = llstar
		lp = lpstar
	} else { //move fossil back to its previous position and restore old branch lengths
		pruneFossilTip(fn)
		graftFossilTip(fn.PAR, lastn)
		lastn.LEN = x
		fn.PAR.LEN = p
	}
	return lp, ll
}

func singleBranchLengthUpdate(ll, lp float64, nodes, inNodes []*Node, tree *Node, branchPrior string) (float64, float64) {
	updateNode := drawRandomNode(nodes)
	soldL := updateNode.LEN
	var propRat float64
	updateNode.LEN, propRat = singleBrlenMultiplierProp(updateNode.LEN)
	llstar := CalcUnrootedLogLike(tree)
	var lpstar float64
	if branchPrior == "0" {
		lpstar = math.Log(1.)
	} else if branchPrior == "1" {
		lpstar = ExponentialBranchLengthLogPrior(nodes, 10.)
	}
	alpha := math.Exp(lpstar-lp) * math.Exp(llstar-ll) * propRat
	s1 := rand.NewSource(time.Now().UnixNano())
	r1 := rand.New(s1)
	r := r1.Float64()
	if r < alpha {
		ll = llstar
		lp = lpstar
	} else {
		updateNode.LEN = soldL
	}
	return lp, ll
}

func drawRandomNode(n []*Node) (rnode *Node) {
	rnoden := rand.Intn(len(n))
	rnode = n[rnoden]
	return
}

func drawRandomReattachment(fn *Node, nodes []*Node) (rnode *Node) {
	s1 := rand.NewSource(time.Now().UnixNano())
	r1 := rand.New(s1)
	rnoden := r1.Intn(len(nodes))
	rnode = nodes[rnoden] //choose a random reattachment point
	if rnode == fn {
		rnode = drawRandomReattachment(fn, nodes)
	} else if rnode == fn.PAR {
		rnode = drawRandomReattachment(fn, nodes)
	}
	return
}

//TODO: should probably combine these two into a single SPR move function
//this removes a fossil tip, along with its parent node from a tree, returning the new branch length left by the gap
//return value is for calculating proposal ratio later
func pruneFossilTip(tip *Node) (x, p float64, lastn *Node) {
	var n *Node
	newpar := tip.PAR
	for _, chd := range newpar.CHLD {
		if chd != tip {
			n = chd
		}
	}
	newpar.PAR.AddChild(n)
	newpar.RemoveChild(n)
	n.PAR = newpar.PAR
	newpar.PAR.RemoveChild(newpar)
	x = n.LEN
	p = newpar.LEN
	lastn = n
	n.LEN = n.LEN + newpar.LEN
	return
}

//this reattaches a fossil tip to branch n at a random point
func graftFossilTip(newpar *Node, n *Node) float64 {
	r := n.LEN
	n.PAR.AddChild(newpar)
	n.PAR.RemoveChild(n)
	newpar.AddChild(n)
	s1 := rand.NewSource(time.Now().UnixNano())
	r1 := rand.New(s1)
	u := r1.Float64()
	newpar.LEN = u * n.LEN
	n.LEN = n.LEN * (1 - u)
	return r
}
