from ete3 import Tree
from random import randint
import argparse
import sys

#root on outgroups
def root_tree(tree,outgroup,verbose=False):
	try:
		if len(outgroup) > 1:
			out_ancestor = tree.get_common_ancestor(outgroup)
			tree.set_outgroup(out_ancestor)
		else:
			tree.set_outgroup(outgroup[0])
	except:
		if verbose:
			print("Could not root tree.")
		tree = "NA"
	return tree

# unroot tree
def unroot_tree(tree):
	for node in tree.traverse():
		if node.is_root():
			node.delete()
	print(tree)
	return tree

# only include taxa in tree
def prune_tree(tree, taxa, preserve_branch_legnths=True):
	tree.prune(taxa, preserve_branch_length=preserve_branch_legnths)
	return tree

def parse_tree(treestring):
	try:
		tree = Tree(treestring)
		return tree
	except:
		tree = "NA"
		print("Parsing error with tree: " + treestring + ". Will call this tree NA.")

def parse_taxafile(taxafile):
	taxa = []
	with open(taxafile) as f:
		for line in f.readlines():
			taxa.append(line.rstrip())
	return taxa

def parse_treefile(treefile):
	trees = []
	with open(treefile) as tf:
		for line in tf.readlines():
			if not line == "":
				trees.append(line.rstrip())
	return trees, len(trees)

def rename_tips(t,samplenames):
	samplenames_dict = {}
	with open(samplenames) as f:
		for line in f.readlines():
			if line.strip() == "":
				continue
			elif line.startswith("#"):
				continue
			oldname,newname = line.strip().split()
			samplenames_dict[oldname] = newname
	for tip in t.traverse():
		o = tip.name
		if tip.name in samplenames_dict.keys():
			tip.name = samplenames_dict[o]
	#for oldname in samplenames_dict.keys():
	#	tip = t.search_nodes(name=oldname)
	#	print(tip.attribute)
	#	#tip.name = samplenames_dict[oldname]
	return t

def getTips(tree):
	tips = []
	for tip in tree.traverse():
		if tip.is_leaf():
			tips.append(tip.name)
	return tips


parser = argparse.ArgumentParser(description="Various small tools to parse/process phylogenetic trees or list of such.")

# arguments for in and output:
parser.add_argument("-i", "--input-tree", help="Input tree(s) in newick format.")
parser.add_argument("--tree-list", help="Input treelist, one tree per line in newick format.")

# some arguments that can be used for both single tree input and treelist
parser.add_argument('--taxa', help="Taxa to keep in output tree, prune out all other.")
parser.add_argument('--remove-tips', help="Taxa to remove from output tree, separated by comma.")
parser.add_argument('--outgroup', help="Outgroup taxon or taxa to root output tree by.")
parser.add_argument('-o', '--outfile', help="File to write output to.", required=True)

# add argument to unrott tree
parser.add_argument('--unroot', action="store_true", help="Removes root node and returns an unrooted tree.")

# 
parser.add_argument('--only-topology', action="store_true", help="Use this to remove branchlengths and support values from the output tree.")
parser.add_argument('--rename-tips', type=str, help="File with two columns, old_name'\t'new_name, given to rename the tip labels of the tree")

#some arguments that can be used on a treelist input
parser.add_argument('--randsample', type=int, help="If specified, N number of trees will be randomly sampled from input treelist.")

args = parser.parse_args()
# parse command line inputs
trees,count = parse_treefile(args.input_tree)

outformat=0
if args.only_topology:
	outformat=9

if args.taxa:
	taxa = parse_taxafile(args.taxa)
	prune = True
	sys.stderr.write('\nSubsetting all trees in input file to keep the following taxa: \n {}.\n'.format(','.join(taxa)))
elif args.remove_tips:
	taxa_rm = args.remove_tips.split(",")
	taxa = [t for t in getTips(parse_tree(trees[0])) if not t in taxa_rm]
	prune = True
	sys.stderr.write('\nSubsetting all trees in input file to remove the following taxa: \n {}.\n'.format(','.join(taxa_rm)))
	sys.stderr.write("\nKeeping the following taxa: \n {}.\n".format(','.join(taxa)))
else:
	prune = False
if args.outgroup:
	outgroup = args.outgroup.split(",")
	root = True
	sys.stderr.write('\nRooting all trees on outgroup: \n {}.\n'.format(','.join(outgroup)))
	# check that, if a set of taxa is given, outgroup is included there
	if args.taxa:
		for o in outgroup:
			if not o in taxa:
				print("\nGiven outgroup must be included in taxa to keep after pruning. Exiting.\n")
				sys.exit()
if args.randsample:
	if count == 1:
		print("\nRandsample can only be used with treelist input, exiting.\n")
		sys.exit()
	else:
		nsample = args.randsample
		randsample = True
else:
	randsample = False

outfile = args.outfile

# check what processings to do and run the functions
if randsample:
	sys.stderr.write('\nWill randomly sample {} trees from input treelist.\n\n'.format(nsample))
	treesample = []
	sampled = []
	for i in range(0,nsample):
		r = randint(0,len(trees))
		while r in sampled:
			r = randint(0,len(trees))
		else:
			t = trees[r]
			treesample.append(t)
			sampled.append(r)
	trees = treesample

# next loop through the treefile and process
processed_trees = []
for tree in trees:
	t = parse_tree(tree)
	if prune:
		t = prune_tree(t,taxa)
	if args.outgroup:
		t = root_tree(t,outgroup)
	if args.unroot:
		t = unroot_tree(t)
	if args.rename_tips:
		t = rename_tips(t,args.rename_tips)
	processed_trees.append(t.write(format=outformat))

# and last write them all to output
sys.stderr.write('\nWriting tree(s) to file {}.\n'.format(args.outfile))
with open(outfile, "w") as of:
	for tree in processed_trees:
		of.write(tree + '\n')
sys.stderr.write('\nDone.\n')





