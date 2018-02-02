package cophy

import (
	"fmt"
	"io/ioutil"
	"strings"
)

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
