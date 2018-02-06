package cophymaru

import (
	"fmt"
	"io/ioutil"
	"math/rand"
	"strings"
	"time"
)

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
		fmt.Print(err)
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
