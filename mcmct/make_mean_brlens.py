from mandos import *
import sys

if len(sys.argv) != 3:
    print "usage: " +sys.argv[0]+"<mcmc output> <tree>"
    sys.exit(0)

tree = tree_utils2.read_tree(sys.argv[2])
for i in tree.iternodes():
    i.length = 0.0

mcmc = open(sys.argv[1],"r")

ldens = []
burn = 100
count = 0
for i in mcmc:
    count += 1
    if count < burn:
        continue
    spls = i.strip().split("\t")
    lens = spls[3:]
    tcount = 0
    for j in tree.iternodes():
        if j == tree:
            continue
        j.length += float(lens[tcount])
        tcount += 1

sgen = count-burn
print sgen 
for n in tree.iternodes():
    if n == tree:
        continue
    n.length = n.length/float(sgen)

print tree.get_newick_repr(True)
