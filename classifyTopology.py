import ete3
from ete3 import Tree
import os
import argparse
import sys
import time

#classify topology according to a number of different criteria, in a similar manner to that of svardal et al. 2019

#this script will only check three groups


#this function is written to parse a file with window topologies, as output by "windowsPhylogenies.py"
#here, the tree will be in the eigth column of the file (7th zero based), and I want to output a list of the trees only, in newick format
#first line is expected to be a header and will be ignored for now

start_time = time.time()

def parse_treefile(treefile, tree_column=8, header=True, sep="\t"):
	treelist=[]
	with open(treefile) as t:
		for nr,line in enumerate(t.readlines()):
			if header:
				if not nr == 0:
					treelist.append(line.split(sep)[tree_column -1].rstrip())
	return treelist

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


#prune tree to only include target taxa
def prune_tree(tree, taxa):
	tree.prune(taxa)
	return tree

def collapsed_leaf(node):
    if len(node2labels[node]) == 1:
       return True
    else:
       return False

def features_from_popfile(clades):
	features = {}
	for clade,sample_list in clades.items():
		for sample in sample_list:
			features[sample] = clade
	return features


def add_leaf_features(tree,features,feature_name="clade"):
	for leaf in tree.get_leaves():
		value = features.get(leaf.name)
		leaf.add_feature(feature_name,value)
	return tree

def monophyletic_branch(node, feature_name="clade"):
	desc = node.get_leaves()
	descending_clades = list(set([c.clade for c in desc]))
	#if nr descending clades is 1, it's a monophyletic branch which we want to collapse
	if len(descending_clades) == 1:
		node.add_feature("name",descending_clades[0]) #changing the label to the clade
		return True
	else:
		node.add_feature("name",descending_clades[0])
		return False


def collapse_topology(tree):
	#copying the tree to keep a raw copy, since the labels will be altered here
	t = tree.copy()
	coll_tree = Tree(t.write(is_leaf_fn=monophyletic_branch))
	coll_tree.ladderize()
	return coll_tree

def get_children_of_split(tree):
	split = {}
	for i,n in enumerate(tree.get_children()):
		split[i] = n.get_leaf_names()
	return split

def get_tree_dummy(tree,first_split,a,b,c):
	tree_dummy = "NA"
	classified = False
	t_class = None
	if len(first_split[0]) == 1 and c in first_split[0]:
		if not c in first_split[1]:
			if len(first_split[1]) == 2:
				tree_dummy = "((" + a + ", " + b + ")," + c + ")"
				t_class = 1
				classified = True
			elif len(first_split[1]) > 2:
				tree_dummy = "((" + a + "/" + b + ")," + c + ")"	
				t_class = 2
				classified = True
		else:
			tree_dummy = "NA"
			classified = False
			t_class = None
	elif len(set(first_split[0])) == 2 and len(set(first_split[1])) == 2:
		splitsets = [set(list(first_split.values())[0]), set(list(first_split.values())[1])]
		if set([a,b]) in splitsets and set([b,c]) in splitsets:
			tree_dummy = "((" + a + "/" + b + "_1), (" + b + "_2/" + c + "))"
			classified = True
			t_class = 3
		else:
			tree_dummy = "NA"
			classified = False
			t_class = None
	return tree_dummy,classified,t_class

#a nest to get the dummy tree, only the last position (c) will matter, so try all of these:
def classify_topology_2(tree,first_split,p1,p2,p3):
	dummy_tree,classified,t_class = get_tree_dummy(tree,first_split,p1,p2,p3)
	if t_class == 1:
		tree_class = 1
		subclass = "1a"
	elif t_class == 2:
		tree_class = 1
		subclass = "1b"
	elif t_class == 3:
		tree_class = 4
		subclass = "4"
	if not classified:
		dummy_tree,classified,t_class = get_tree_dummy(tree,first_split,p3,p2,p1)
		if t_class == 1:
			tree_class = 2
			subclass = "2a"
		elif t_class == 2:
			tree_class = 2
			subclass = "2b"
		elif t_class == 3:
			tree_class = 4
			subclass = "4"
		if not classified:
			dummy_tree,classified,t_class = get_tree_dummy(tree,first_split,p3,p1,p2)
			if t_class == 1:
				tree_class = 3
				subclass = "3a"
			elif t_class == 2:
				tree_class = 3
				subclass = "3b"
			elif t_class == 3:
				tree_class = 5
				subclass = "5"
	if not classified:
		tree_class = "5"
		subclass = "5"
	if tree_class == 5:
		dummy_tree = "NA"
	return tree_class,subclass,dummy_tree
	
def parse_treeline(line,chromcol=1,startcol=2,endcol=3,sitescol=5,treecol=7):
	chrom = line.split("\t")[chromcol]
	start = line.split("\t")[startcol]
	stop = line.split("\t")[endcol]
	sites = line.split("\t")[sitescol]
	tree = line.split("\t")[treecol].strip()
	if tree == "":
		tree = "NA"
	return chrom,start,stop,sites,tree

def parse_tree(treestring):
	try:
		tree = Tree(treestring)
		return tree
	except:
		tree = "NA"
		print("Parsing error with tree: " + treestring + ". Will call this tree NA.")
		print(tree)
		return tree

def all_clades_represented(tree,clades):
	taxa_counts = {p1: 0, p2: 0, p3: 0}
	for leaf in tree.get_leaves():
		if leaf.name in clades[p1]:
			taxa_counts[p1] += 1
		elif leaf.name in clades[p2]:
			taxa_counts[p2] += 1
		elif leaf.name in clades[p3]:
			taxa_counts[p3] += 1
	#if any of the groups are missing, we set the tree to "NA"
	if 0 in list(taxa_counts.values()):
		return False
	else:
		return True


if __name__ == "__main__":

	parser = argparse.ArgumentParser()

	#arguments for in and output:
	parser.add_argument("-t", "--treelist", help="Input file, containing trees from genomic regions where the trees should be in the 8th column")
	parser.add_argument("-o", "--output", help="Output prefix, will render two files.", required=True, metavar="output")

	parser.add_argument("--clades", help="File with two tab separated columns, sample in the first and clade/other attribute in second.", required=True)
	parser.add_argument("--outgroup", help="Population/taxon to use as outgroup, as specified in the clades file.", default="Outgroup")
	parser.add_argument("--p1", help="Sister to p2 in the species tree.")
	parser.add_argument("--p2", help="Target clade to investigate, will be interpreted as a sister to p1 and introgressed with p3.")
	parser.add_argument("--p3", help="Population/taxon hypothezised to have introgressed with p2.")
	parser.add_argument("--chromcol", help="column of input file containing the chromosome name", default=1)
	parser.add_argument("--startcol", help="column of input file containing the start coordinate name", default=2)
	parser.add_argument("--endcol", help="column of input file containing the end coordinate name", default=3)
	parser.add_argument("--sitescol", help="column of input file containing the number of sites used in analysis", default=6)
	parser.add_argument("--treecol", help="column of input file containing the tree", default=7)

	args = parser.parse_args()

	prefix = args.output
	input_treelist = args.treelist
	p2 = args.p2
	p1 = args.p1
	p3 = args.p3
	popfile = args.clades
	outgroup = args.outgroup
	#parse the popfule
	clades = {}
	with open(popfile) as f:
		for line in f.readlines():
			sample = line.split("\t")[0]
			clade = line.split("\t")[1].rstrip()
			if not clade in list(clades.keys()):
				clades[clade] = [sample]
			else:
				clades[clade].append(sample)

	outgroup = clades[outgroup]

	#open csv files for writing the output
	window_output = open(os.path.join(prefix + "_windowtrees.txt"), "w")
	window_output.write("Chrom\tStart\tStop\tSites\tTreeClass\tTreeSubClass\tCollapsedTree\tFullTree\n")

	#empty list to append the trees
	treelist = []

	#input_treelist="/Users/axeljensen/WINDOW_TOPOLOGIES/testtop.csv"
	groups = [clade for clade in clades.keys()]
	p2_taxa = clades[p2]
	p1_taxa = clades[p1]
	p3_taxa = clades[p3]
	target_taxa = p1_taxa + p2_taxa + p3_taxa
	
	#dictionary for storing the tree-legend
	tree_legend = {'5': 'other'}
	#loop through all the lines in the input treefile to categorize them
	previous_chrom = None #just for outputting some progression
	with open(input_treelist) as qt:
		for line in qt.readlines():
			#skip the header
			if not line.startswith("window"):
				chrom,start,stop,sites,tree = parse_treeline(line,args.chromcol,args.startcol,args.endcol,args.sitescol,args.treecol)
				#get some of the metadata
				current_chrom = chrom
				if previous_chrom:
					if not current_chrom == previous_chrom:
						print("Finished with chromosome " + previous_chrom + ". Now starting with " + chrom + ".")
				else:
					print("Starting classification, now working on chromosome: " + chrom + ".")
				if not tree == "NA" or tree == "" or tree == ".": #if tree is NA already in the raw file, write it directly as NA
					#parse treestring with ete3
					t = parse_tree(tree)
					#root with outgroup(s)
					t = root_tree(t,outgroup)
					#prune to only include target taxa
					if not t == "NA":
						t = prune_tree(t, target_taxa)
						#make a {'sample':'clade'} dictionary from the popfile
						features = features_from_popfile(clades)
						#and add those to the tree
						t = add_leaf_features(t,features,feature_name="clade")
						#collapse monophyletic branches
						t_collapsed = collapse_topology(t)
						#get the children of the deepest split for catogrization
						first_split = get_children_of_split(t_collapsed)
						#then have a go at classifying the topology
						#count the leaves for each group, if any is empty we discard this tree
						tree_class,subclass,dummy_tree =  classify_topology_2(t_collapsed,first_split,p1,p2,p3)
						if not subclass in tree_legend.keys() and not tree_class == 5:
							tree_legend[subclass] = dummy_tree
						#check monophyly
						#tree_class = classify_topologies(t,clades,p2,p1,p3)
						#write it to file
						nwktree_full = t.write(format=9)
						nwktree_collapsed = t_collapsed.write(format=9)
					if tree == "NA" or tree == "":
						tree_class = "NA"
						subclass = "NA"
						nkwtree_full = "NA"
						nwktree_collapsed = "NA"
					try:
						window_output.write(chrom + "\t" + str(start) + "\t" + str(stop) + "\t" + str(sites) + "\t" + str(tree_class) + "\t" + str(subclass) + "\t" + nwktree_collapsed + "\t" + nwktree_full + "\n")
					except:
						sys.exit()
				previous_chrom = chrom

	print("Writing output.")
	window_output.close()

	trees_out = open(os.path.join(prefix + "_tree_ids.txt"), "w")
	trees_out.write("Tree_class\tTopology\n")
	for c,tree in tree_legend.items():
		trees_out.write(c + "\t" + tree + "\n")
	trees_out.close()
	used_minutes = time.time() - start_time
	print("Finished.")
