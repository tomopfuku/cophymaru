package main

import (
	"flag"
	"fmt"
	"math"
	"math/rand"
	"os"
	"time"

	"github.com/carolinetomo/gophy"
)

func main() {
	rand.Seed(time.Now().UTC().UnixNano())
	//lg := bufio.NewWriter(f)
	tfn := flag.String("t", "", "tree filename")
	afn := flag.String("s", "", "seq filename")
	rfn := flag.String("r", "", "file containing temporal ranges")
	wks := flag.Int("w", 4, "number of threads")
	flag.Parse()
	if len(*tfn) == 0 {
		fmt.Fprintln(os.Stderr, "need a tree filename (-t)")
		os.Exit(1)
	}
	if len(*afn) == 0 {
		fmt.Fprintln(os.Stderr, "need a seq filename (-s)")
		os.Exit(1)
	}

	//read a tree file
	//trees := gophy.ReadTreesFromFile(*tfn)
	nwk := gophy.ReadLine(*tfn)[0]
	//rt := cophymaru.ReadTree(nwk)
	rt := gophy.ReadNewickString(nwk)
	t := gophy.NewTree()
	t.Instantiate(rt)
	//read a file of stratgraphic ranges
	gophy.ReadStrat(*rfn, t)
	gophy.MakeStratHeights(t)
	//fmt.Println(t.Rt.Newick(true))
	stratLnL := gophy.PoissonTreeLoglike(t)
	//read a seq file
	nsites := 0
	seqs := map[string][]string{}
	mseqs, numstates := gophy.ReadMSeqsFromFile(*afn)
	seqnames := make([]string, 0)
	for _, i := range mseqs {
		seqs[i.NM] = i.SQs
		seqnames = append(seqnames, i.NM)
		nsites = len(i.SQ)
	}
	x := gophy.NewMultStateModel()
	x.NumStates = numstates
	x.SetMap()
	bf := gophy.GetEmpiricalBaseFreqsMS(mseqs, x.NumStates)
	x.SetBaseFreqs(bf)
	x.EBF = x.BF
	x.SetEqualBF()
	//patterns, patternsint, gapsites, constant, uninformative, _ := gophy.GetSitePatternsMS(mseqs, x)
	_, patternsint, _, _, _, _ := gophy.GetSitePatternsMS(mseqs, x)
	patternval, _ := gophy.PreparePatternVecsMS(t, patternsint, seqs, x)
	//sv := gophy.NewSortedIdxSlice(patternvec)
	//sort.Sort(sv)
	x.SetupQJC()
	//fmt.Println(x.Q, x.BF)
	l := gophy.PCalcLogLikePatternsMS(t, x, patternval, *wks)
	gophy.PCalcSankParsPatternsMultState(t, x, patternval, 1)
	gophy.EstParsBLMultState(t, x, patternval, nsites)
	for _, n := range t.Post {
		n.Len = math.Max(10e-10, n.Len/(float64(nsites)))
	}
	//fmt.Println(t.Rt.Newick(true))
	l = gophy.PCalcLogLikePatternsMS(t, x, patternval, *wks)

	fmt.Println(l, stratLnL, l+stratLnL)
	fmt.Println(t.Rt.NewickChronogram())
	os.Exit(0)
	gophy.OptimizeBLNRMS(t, x, patternval, *wks)
	fmt.Println(t.Rt.Newick(true))
	l = gophy.PCalcLogLikePatternsMS(t, x, patternval, *wks)
	fmt.Println(l)

}
