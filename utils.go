package cophymaru

import (
	"fmt"
	"io/ioutil"
	"math/rand"
	"os"
	"strings"
	"time"
)

//InternalNodeSlice will return a slice containing only internal nodes
func InternalNodeSlice(nodes []*Node) (inNodes []*Node) {
	for _, n := range nodes {
		if len(n.Chs) == 2 {
			inNodes = append(inNodes, n)
		}
	}
	return
}

//InitParallelPRNLen will set up empty slices for the prnlens
func InitParallelPRNLen(nodes []*Node) {
	for _, n := range nodes {
		n.CONPRNLen = make([]float64, len(nodes[0].CONTRT))
	}
}

//TreeLength will return the total length of a slice of nodes
func TreeLength(nodes []*Node) float64 {
	len := 0.
	for _, n := range nodes[1:] {
		len += n.Len
	}
	return len
}

//MakeRandomStartingBranchLengths will initialize a tree with a set of random branch lengths
func MakeRandomStartingBranchLengths(tree *Node) {
	nodes := tree.PreorderArray()
	for _, n := range nodes {
		s1 := rand.NewSource(time.Now().UnixNano())
		r1 := rand.New(s1)
		u := r1.Float64()
		n.Len = u
	}
}

//ReadLine is like the Python readline() and readlines()
func ReadLine(path string) (ln []string) {
	b, err := ioutil.ReadFile(path)
	if err != nil {
		fmt.Println(err)
		fmt.Println("There was an error when reading in the file:", path, ". Are you sure that it exists?")
		os.Exit(0)
	}
	ss := string(b)
	ln = strings.Split(ss, "\n")
	return
}

//ReadFossils will read in a list of fossil tips one line at a time into a slice
//TODO: get this working
func ReadFossils(path string) (fos []string) {
	l := ReadLine(path)
	for _, i := range l {
		fmt.Println(i)
	}
	return
}

func maxChildFAD(node *Node) float64 {
	highest := 0.0
	for _, c := range node.Chs {
		if c.FAD > highest {
			highest = c.FAD
		}
	}
	return highest
}

func AssignBranchRates(preTree []*Node, rates []float64) bool {
	/*
		if len(preTree)-1 != len(rates) {
			fmt.Println("You are trying to estimate a different number of rates than there are branches")
			os.Exit(0)
		}
	*/
	i := 0
	for _, node := range preTree[1:] {
		if node.Nam == "root" {
			continue
		} else if node.ISTIP == true && node.ANC == true {
			continue
		} else {
			rate := rates[i]
			if rate < 0.0 {
				return true
			}
			node.RATE = rate
			i++
		}
	}
	return false
}

func AssignGlobalRate(preTree []*Node, rate float64) bool {
	if rate < 0.0 {
		return true
	}
	for _, node := range preTree {
		if node.ANC == true && node.ISTIP == true {
			continue
		}
		node.RATE = rate
	}
	return false
}

func AssignInternalNodeHeights(preTree []*Node, heights []float64) bool {
	count := 0
	for _, node := range preTree {
		if node.ISTIP == false {
			newheight := heights[count]
			var constr1, constr2 float64
			constr1 = OldestChildAge(node)
			if node.ANC == false {
				constr2 = constr1
			} else {
				//if newheight > node.FAD {
				//	fmt.Println(node.Nam, node.FAD, newheight)

				//return true
				//}
				cont := node.GetContDesc()
				constr2 = cont.LAD
			}
			if newheight < constr1 || newheight < constr2 {
				//fmt.Println(node.Nam, node.Height, newheight, node.FAD, constr1, constr2)
				return true
			}
			if node.Nam != "root" {
				if newheight > node.Par.FAD {
					//fmt.Println(node.Nam)
					return true
				}
			}
			node.Height = newheight
			if node.ANC == false {
				node.FAD = node.Height
			}
			count++
		}
		if node.Nam != "root" {
			newlen := node.Par.Height - node.Height
			if newlen < 0.0 {
				return true
			}
			node.Len = newlen
			if node.Height > node.Par.Height {
				return true
			}
		}
	}
	return false
}

func NodeFromLabel(label string, nodes []*Node) *Node {
	for _, node := range nodes {
		if node.Nam == label {
			return node
		}
	}
	return nil
}

func AIC(lnl, k float64) (aic float64) {
	aic = (2 * k) - (2 * lnl)
	return
}
