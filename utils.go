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
		if len(n.CHLD) == 2 {
			inNodes = append(inNodes, n)
		}
	}
	return
}

//InitParallelPRNLEN will set up empty slices for the prnlens
func InitParallelPRNLEN(nodes []*Node) {
	for _, n := range nodes {
		n.CONPRNLEN = make([]float64, len(nodes[0].CONTRT))
	}
}

//TreeLength will return the total length of a slice of nodes
func TreeLength(nodes []*Node) float64 {
	len := 0.
	for _, n := range nodes[1:] {
		len += n.LEN
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
		n.LEN = u
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
	for _, c := range node.CHLD {
		if c.FAD > highest {
			highest = c.FAD
		}
	}
	return highest
}

func AssignBranchRates(preTree []*Node, rates []float64) bool {
	if len(preTree)-1 != len(rates) {
		fmt.Println("You are trying to estimate a different number of rates than there are branches")
		os.Exit(0)
	}
	for i, node := range preTree[1:] {
		if node.NAME == "root" {
			continue
		} else {
			rate := rates[i]
			if rate < 0.0 {
				return true
			}
			node.RATE = rate
		}
	}
	return false
}

func AssignGlobalRate(preTree []*Node, rate float64) bool {
	if rate < 0.0 {
		return true
	}
	for _, node := range preTree {
		node.RATE = rate
	}
	return false
}

func AssignInternalNodeHeights(preTree []*Node, heights []float64) bool {
	count := 0
	for _, node := range preTree {
		if len(node.CHLD) != 0 {
			node.HEIGHT = heights[count]
			constr := maxChildFAD(node)
			if node.HEIGHT <= constr {
				return true
			}
			count++
		}
		if node.NAME != "root" {
			node.LEN = node.PAR.HEIGHT - node.HEIGHT
			if node.HEIGHT >= node.PAR.HEIGHT {
				return true
			}
		}
		if node.LEN < 0.0 {
			return true
		}
	}
	return false
}
