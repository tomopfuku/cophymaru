package cophymaru

import (
	"fmt"
	"os"
)

//AssertNumChars will check to make sure that two character vectors are the same lengths (as they always should be)
func AssertNumChars(chars1, chars2 []float64) {
	if len(chars1) != len(chars2) {
		fmt.Print("ERROR. Something went wrong with the mapping of traits to the tips.")
		os.Exit(0)
	}
}
