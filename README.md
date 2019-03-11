# cophymaru
_cophymaru_ is a phylogenetics package that uses continuous character data to 1) recover the placement of fossil taxa on reference trees comprised of extant species (often inferred from molecular data) 2) infer phylogenetic relationships from scratch using continuous traits, and 3) reconstruct phylogeny, including ancestor-descendant relationships, between fossil taxa using continuous traits and stratigraphic ranges. 


To install, you will need to have Golang installed on your machine. This is simple, and instructions can be found at: https://golang.org/doc/install. 
To compile, clone this repository in your Golang src directory, cd into the repo root, and type:

        go build mcmct/maru.go
        go build goMandos/mandos.go

You can then move the binaries into your path or run it from there. _cophymaru_ is fairly self-contained, and so there shouldn't be major issues with compiling.

type:

        maru -h 

or

        mandos -h

to reveal the set of options that specifies data inputs, type of analysis, name(s) of fossil(s), branch length prior, MCMC/tree search settings, number of threads, and other parameters. 

