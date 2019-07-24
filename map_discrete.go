package cophymaru

import (
	"fmt"
	"os"
	"strings"
)

// ReadDiscrete will read in a PHYLIP-formatted trait file and return a map
func ReadDiscrete(traitfl string) map[string][]string {
	traitlines, ntax, ntraits := readPhylipHeader(traitfl)
	tm := make(map[string][]string, ntax)
	//lnum := 0
	//var ntax,ntraits int
	var curtr string
	for _, line := range traitlines {
		if line == "" {
			continue
		}
		ss := strings.Split(line, "\t")
		curtax := ss[0]
		var curtraits []string
		seq := ss[1]
		if len(seq) != ntraits {
			fmt.Println("the number of traits specified in the header differs from the number present in this line:")
			fmt.Println(line)
		}
		for _, c := range seq {
			v := string(c)
			if v != "?" && v != "-" {
				curtr = v
			} else {
				curtr = "-"
			}
			curtraits = append(curtraits, curtr)
		}
		tm[curtax] = curtraits
	}
	var curtraits []string
	i := 0
	for {
		if i == ntraits {
			break
		}
		curtraits = append(curtraits, "-")
		i++
	}
	tm["internal"] = curtraits
	return tm
}

func makeSiteProbsMap(traits map[string][]string) map[string]map[int]map[string]float64 {
	//maxNumSites := 0
	uniqueTraits := make(map[string]bool)
	for _, v := range traits {
		for _, tr := range v {
			if tr == "-" || tr == "?" {
				continue
			}
			if _, ok := uniqueTraits[tr]; !ok {
				uniqueTraits[tr] = true
			}
		}
	}
	siteProbsMap := make(map[string]map[int]map[string]float64)
	for tax, v := range traits {
		siteProbs := make(map[int]map[string]float64)
		for i, tr := range v {
			probs := make(map[string]float64)
			for k := range uniqueTraits {
				if k == tr {
					probs[k] = 1.0
				} else {
					probs[k] = 0.0
				}
			}
			siteProbs[i] = probs
		}
		siteProbsMap[tax] = siteProbs
	}
	return siteProbsMap
}

//MapDiscrete maps the traits contained with in a (golang) map to the tips of a tree and initializes slices of the same length for the internal nodes
func MapDiscrete(t *Node, traits map[string][]string) {
	siteProbsMap := makeSiteProbsMap(traits)
	for _, n := range t.PreorderArray() {
		if n.ISTIP {
			if _, ok := traits[n.Nam]; !ok {
				fmt.Println("No traits provided for ", n.Nam)
				os.Exit(0)
			}
			n.DISCTRT = siteProbsMap[n.Nam]
		} else {
			n.DISCTRT = make(map[int]map[string]float64)
			for site, probs := range siteProbsMap["internal"] {
				n.DISCTRT[site] = make(map[string]float64)
				for state, prob := range probs {
					n.DISCTRT[site][state] = prob
				}
			}
		}
		for range n.DISCTRT {
			n.LL = append(n.LL, 0.0)
			n.CONPRNLen = append(n.CONPRNLen, 0.0)
		}
	}
}
