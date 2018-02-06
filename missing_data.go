package cophymaru

//MakeMissingMeansArray will find the missing sites in a data matrix and plug in values corresponding to the mean of the remaining sites
func MakeMissingMeansArray(tree []*Node) {
	means := CalcSiteMeans(tree)
	for _, n := range tree {
		MakeMissingMeansTip(n, means)
	}
}

//MakeMissingMeansTip will replace missing values with the mean across all tips for a single tip
func MakeMissingMeansTip(n *Node, means []float64) {
	for i := range n.CONTRT {
		if n.MIS[i] == true {
			n.CONTRT[i] = means[i]
		}
	}
}

//CalcSiteMeans will calculate the mean value for all the sites in the matrix for which the site is not missing
func CalcSiteMeans(nodes []*Node) (siteSum []float64) {
	var ntraits []int
	for range nodes[0].CONTRT {
		siteSum = append(siteSum, 0.0)
		ntraits = append(ntraits, 0)
	}
	for _, n := range nodes {
		for i, tr := range n.CONTRT {
			if n.MIS[i] != true {
				siteSum[i] += tr
				ntraits[i]++
			}
		}
	}
	for i, m := range siteSum {
		m = m / float64(ntraits[i])
	}
	return
}
