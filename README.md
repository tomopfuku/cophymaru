# cophymaru
_cophymaru_ is a phylogenetics package that uses continuous character data to recover the placement of fossil taxa on reference trees comprised of extant species (often inferred from molecular data). 

type:

        maru -h 

to reveal the set of options that specifies the name(s) of fossil(s), branch length prior, MCMC settings, number of threads, and other parameters. 

To install, you will need to have Golang installed on your machine. This is simple, and instructions can be found at: https://golang.org/doc/install. 
To compile, clone this repository in your Golang src directory, cd into the repo root, and type:

        go build mcmct/maru.go

You can then move the binary into your path or run it from there. _cophymaru_ is fairly self-contained, and so there shouldn't be major issues with compiling.

