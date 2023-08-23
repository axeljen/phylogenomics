from ete3 import Tree
import sys
import argparse
import os
sys.path.append("/home/axeljen/phylogenomics/")
import functions as fn


#testtree = "(((((((((((((PD_0028_Allenopithecus_nigroviridis:0.00188498,PD_0126_Allenopithecus_nigroviridis_ssp:0.00135842):0.01032673,(((((PD_0029_Allochrocebus_lhoesti:0.00065189,PD_0127_Allochrocebus_lhoesti_ssp:0.00127078):0.00084073,PD_0031_Allochrocebus_solatus:0.00268690):0.00155568,PD_0030_Allochrocebus_preussi:0.00461430):0.00447330,((PD_0173_Chlorocebus_pygerythrus_ssp:0.00324368,SAMN02692280_Chlorocebus_sabaeus_saintkitts:0.00377091):0.00318324,RD-1790_Cercopithecus_dryas:0.00606934):0.00040403):0.00128171,(PD_0034_Erythrocebus_patas:0.00117352,PD_0062_Erythrocebus_patas:0.00074588):0.00601962):0.00288441):0.00071209,(PD_0408_Cercocebus_torquatus:0.00818847,SAMN03264642_Macaca_mulatta:0.00780022):0.00968107):0.00182014,PD_0124_Miopithecus_ogouensis_ssp:0.01252318):0.00212702,((GP612_C_hamlyni:0.00195005,PD_0060_Cercopithecus_hamlyni:0.00234551):0.00123626,ME408_C_lomamiensis:0.00274149):0.00668384):0.00160789,(((PD_0032_Cercopithecus_ascanius_schmidti:0.00356239,PD_0066_Cercopithecus_ascanius:0.00204417):0.00260447,((((PD_0041_Cercopithecus_cephus_ssp:0.00401567,PD_0095_Cercopithecus_cephus:0.00462199):0.00004710,PD_0094_Cercopithecus_cephus:0.00494403):0.00111339,PD_0406_Cercopithecus_ascanius_ssp:0.00569839):0.00098149,PD_0091_Cercopithecus_petaurista:0.00640799):0.00060767):0.00064834,PD_0410_Cercopithecus_ascanius_ssp:0.00571126):0.00062537):0.00080749,((ME403_C_mitis_stuhlmanni:0.00442156,((PD_0061_Cercopithecus_albogularis:0.00332575,SAMN25565059_C_mitis:0.00294347):0.00103008,PD_0096_Cercopithecus_mitis:0.00439150):0.00023660):0.00107886,((PD_0033_Cercopithecus_nictitans_nictitans:0.00246579,SAMN13555794_Cercopithecus_nictitans:0.00279497):0.00017341,PD_0090_Cercopithecus_nictitans:0.00332951):0.00050423):0.00085588):0.00046582,((PD_0058_Cercopithecus_neglectus:0.00090764,PD_0407_Cercopithecus_neglectus:0.00034714):0.00458227,(PD_0059_Cercopithecus_diana:0.00267265,PD_0093_Cercopithecus_roloway:0.00504714):0.00239151):0.00120431):0.00092066,((((FK104_C_denti:0.00057239,JH010_C_wolfi_elegans:0.00070424):0.00005939,(GP611_C_wolfi_elegans:0.00098728,PPN005_C_wolfi_wolfi:0.00048957):0.00018499):0.00196082,JH026_C_denti:0.00213712):0.00097186,PD_0092_Cercopithecus_pogonias:0.00248385):0.00227256):0.00123399,(PD_0043_Cercopithecus_lowei_ssp:0.00591324,(SAMN13555781_Cercopithecus_mona:0.00302283,SAMN13555792_Cercopithecus_mona:0.00211948):0.00082420):0.00092359):0.00085430,((PD_0042_Cercopithecus_mona_ssp:0.00173420,((SAMN13555780_Cercopithecus_mona:0.00169894,(SAMN13555786_Cercopithecus_mona:0.00155401,(SAMN13555789_Cercopithecus_mona:0.00066727,SAMN13555793_Cercopithecus_mona:0.00175784):0.00018595):0.00063469):0.00054268,SAMN13555790_Cercopithecus_mona:0.00220023):0.00037665):0.00017171,((SAMN13555783_Cercopithecus_mona:0.00161123,(SAMN13555784_Cercopithecus_mona:0.00107329,SAMN13555787_Cercopithecus_mona:0.00061548):0.00051663):0.00060589,SAMN13555791_Cercopithecus_mona:0.00171031):0.00067788):0.00006749):0.00024272,SAMN13555785_Cercopithecus_mona:0.00243443):0.00031268,SAMN13555782_Cercopithecus_mona:0.00181860,SAMN13555788_Cercopithecus_mona:0.00154059);"
#tree = Tree(testtree)

#popfile = "/crex/proj/sllstore2017021/nobackup/DENTI_LESULA_PROJECT/ANALYSES/NJ_TREES_WINDOWS/popfile_species.txt"

#pops = fn.parsePopfile(popfile)
#tree = Tree(testtree)
#outgroup = tree.get_common_ancestor(outgroup_taxa)
#tree.set_outgroup(outgroup)
# this function collapses groups in the popfile into a single branch, if they are monophyletic
def collapseGroups(tree, pops):
	# loop over the groups in the popfile
	for pop in pops:
		# set the group attribute to the samples in the group
		for leaf in tree.get_leaves():
			if leaf.name in pops[pop]:
				leaf.add_feature("group", pop)
	# check if they are monophyletic, collapse those that are
	for pop in pops:
		# if there is only a single sample in the group
		if len(pops[pop]) == 1:
			# just rename the leaf
			tree.search_nodes(name=pops[pop][0])[0].name = pop
		# if there are more than one sample in the group check if they are monophyletic
		elif tree.check_monophyly(values = [pop], target_attr="group")[0]:
			# if they are, get this node
			n = tree.get_common_ancestor(pops[pop], target_attr="name")
			# just keep a single leaf
			n.get_leaves()[0].name = pop
			# and collapse by deleting the rest of children
			if len(n.get_leaves()) > 1:
				for child in n.get_children():
					if child.name != pop:
						child.delete()
		else:
			print("Warning: group {} is not monophyletic.".format(pop))
			#print(tree)
			# if a group is not mono, return NA
			return "NA"
	#print(tree)
	return tree


# prep a parser for the command line arguments
parser = argparse.ArgumentParser(description='Takes a list of trees, and classifies them according to unique tree topologies.')

parser.add_argument('-t', '--trees', type=str, required=False, help='File with trees to classify. In this script Im feeding a table with trees in the 8th column, having some additional metadata in the first 7 columns. Some of it is used here.')
parser.add_argument('-o', '--out', type=str, required=True, help='Output prefix.')

# adding a function to combine multiple runs across different chromosomes or genomic regions, in which case only the treefile and the output prefix are needed
parser.add_argument('-c', '--combine', type=str, default=None, required=False, help='If given, the script will combine multiple runs across different chromosomes or genomic regions. Input file should be a two-column file, with the path to first the summary and then the long input, for each run to combine.')

# popfile with samples and groups for the samples to consider
parser.add_argument('-p', '--popfile', type=str, required=False, help='Popfile with samples and groups for the samples to consider. If more than one sample is in a group, only trees where they are monophyletic will be considered.')

# outgroup sample(s)
parser.add_argument('-g', '--outgroup', type=str, required=False, help='Outgroup sample(s).')

# metacols
parser.add_argument('--meta', type=str, help="Headers of the columns leading up to the tree, which should always be in the last column. If at least one header is given, it will be used in the output file. Number of headers must match the actual number of columns in the input file. If no headers are given, the input file is assumed to be a tree file only, and no metadata will be added to the output file.", default='number,chrom,start,end,lnL,input_tree_taxa,nsites')

# species tree if wanted
parser.add_argument('--species_tree', type=str, help="Species tree to compare the found trees. If given, the script will also calculate the RF distance between the species tree and the found trees.", default=None)

# parse the arguments
args = parser.parse_args()


# if combine is not given, check that the popfile is given
if args.combine == None:
	if args.popfile == None:
		print("Error: popfile is required in a regular run.")
		sys.exit(1)
	# check that treefile is given
	if args.trees == None:
		print("Error: treefile is required in a regular run.")
		sys.exit(1)
	
	# read the popfile
	pops = fn.parsePopfile(args.popfile)

	# get a list of all the samples in the popfile
	samples = []
	for pop in pops:
		samples += pops[pop]

	# list of dictionaries to store the trees and their metadata
	trees = []

	# and a list of unique trees for later on
	unique_trees = []

	counter = 0
	total = len(open(args.trees).readlines())

	print("Starting to process {} trees.".format(total))

	# read the tree file
	with open(args.trees, "r") as f:
		headers = args.meta.split(",")
		headers.append("tree")
		for line in f:
			if 'tree' in line.strip().split():
				continue # skip the header line
			else:
				values = line.strip().split()
				d = dict(zip(headers, values))
				# parse the tree
				if 'tree' in d.keys():
					tree = Tree(d.get('tree'))
					tree.prune(samples)
					outgroup_taxa = args.outgroup.split(",")
					if len(outgroup_taxa) > 1:
						outgroup = tree.get_common_ancestor(outgroup_taxa)
					else:
						outgroup = tree.search_nodes(name=outgroup_taxa[0])[0]
					try:
						tree.set_outgroup(outgroup)
					except:
						#print("Warning: tree {} could not be parsed. Skipping tree.".format(d))
						tree = "NA"
						d['tree'] = 'NA'
				else:
					#print("Warning: tree {} could not be parsed. Skipping tree.".format(d))
					tree = "NA"
					d['tree'] = 'NA'
				#  if the tree is not NA, collapse the groups in the popfile
				if tree != "NA":
					tree = collapseGroups(tree, pops)
					# if this returns NA, print a warning
					#if tree == "NA":
					#	print("Warning: tree {} not monophyletic for all groups in popfile. Skipping tree.".format(d['tree']))
				# now, if tree is NA at this point, just add it to the dictionary
				if tree == "NA":
					d['parsed_tree'] = "NA"
					d['tree_category'] = "NA"
					trees.append(d)
				else:
					# if there are no unique trees so far, add this one
					if len(unique_trees) == 0:
						unique_trees.append(tree)
						d['parsed_tree'] = tree.write(format=9)
						d['tree_category'] = 1
						trees.append(d)
					else:
						# check against the unique trees
						for i in range(len(unique_trees)):
							rf = tree.robinson_foulds(unique_trees[i])[0]
							if rf == 0:
								# this means we found it
								d['parsed_tree'] = tree.write(format=9)
								d['tree_category'] = i+1
								trees.append(d)
								break
							else:
								# if we didn't find it, and we're at the last tree, add this tree to the unique trees
								if i == len(unique_trees)-1:
									unique_trees.append(tree)
									d['parsed_tree'] = tree.write(format=9)
									d['tree_category'] = i+2
									trees.append(d)
									break
			counter += 1
			if counter % 100 == 0:
				print("Processed {} trees.".format(counter))
				print("found {} unique trees so far, out of a total of {} parsed.".format(len(unique_trees), len([i for i in trees if i['parsed_tree'] != "NA"])))

	# remove the original tree from the dictionary
	for tree in trees:
		del tree['tree']

	# write the output file
	with open(args.out + ".full.txt", "w") as f:
		f.write("\t".join(trees[0].keys()) + "\n")
		for tree in trees:
			f.write("\t".join([str(tree[i]) for i in tree]) + "\n")


	# now we make a dictionary with the unique trees and their counts, and if a species tree was given, we calculate the RF distance between the species tree and the unique trees
	unique_trees_dict = {}
	if args.species_tree != None:
		st = Tree(args.species_tree)
		# prune the species tree
		st.prune(samples)
		# root the species tree
		if len(outgroup_taxa) > 1:
			outgroup = st.get_common_ancestor(outgroup_taxa)
		else:
			outgroup = st.search_nodes(name=outgroup_taxa[0])[0]
		st.set_outgroup(outgroup)
		# collapse the groups in the same way as the trees
		st = collapseGroups(st, pops)

	for tree in trees:
		if tree['parsed_tree'] != "NA":
			if tree['tree_category'] in unique_trees_dict.keys():
				unique_trees_dict[tree['tree_category']]['count'] += 1
			else:
				unique_trees_dict[tree['tree_category']] = {'count':1, 'tree':unique_trees[tree['tree_category'] -1].write(format=9)}
	# if a species tree was given, calculate the RF distance
	for tree in unique_trees_dict:
		if args.species_tree != None:
			rf = Tree(unique_trees_dict[tree]['tree']).robinson_foulds(st)[0]
			unique_trees_dict[tree]['rf'] = rf
		else:
			unique_trees_dict[tree]['rf'] = "NA"

	# write the output file
	with open(args.out + ".summary.txt", "w") as f:
		f.write("tree_cat\ttreetop\tcount\trf\n")
		for tree in unique_trees_dict:
			f.write("{}\t{}\t{}\t{}\n".format(tree, unique_trees_dict[tree]['tree'], unique_trees_dict[tree]['count'], unique_trees_dict[tree]['rf']))

else:
	if os.path.isfile(args.combine):
		print("Combining the runs listed in file {}".format(args.combine))
		# first column contains the summary file, second column contains the long input file
		runs = [i.strip().split() for i in open(args.combine).readlines()]
		unique_trees = []
		# a dict for translating trees between runs
		translation = {r[0]: {} for r in runs}
		# fetch all the unique trees from first run and append to unique_trees
		with open(runs[0][0]) as f:
			for i,line in enumerate(f.readlines()):
				if i == 0:
					continue# skip the header line
				else:
					treecat, treetop, count, rf = line.strip().split()
					parsed_tree = Tree(treetop)
					unique_trees.append({'cat':i, 'tree': treetop, 'parsed_tree': parsed_tree, 'count': int(count), 'rf': rf,})
					translation[runs[0][0]][treecat] = i
		# sort the unique trees by count
		unique_trees = sorted(unique_trees, key=lambda k: int(k['count']), reverse=True)
		# now loop over the rest of the runs
		for run in runs[1:]:
			# read the summary file
			with open(run[0]) as f:
				for i,line in enumerate(f.readlines()):
					if i == 0:
						continue
					else:
						found = False
						treecat, treetop, count, rf = line.strip().split()
						parsed_tree = Tree(treetop)
						for tree in unique_trees:
							if tree['parsed_tree'].robinson_foulds(parsed_tree)[0] == 0:
								tree['count'] += int(count)
								translation[run[0]][treecat] = tree['cat']
								found = True
								break
						if not found:
							newcat = len(unique_trees)+1
							unique_trees.append({'cat':newcat, 'tree': treetop, 'parsed_tree': parsed_tree, 'count': int(count), 'rf': rf,})
							translation[run[0]][treecat] = newcat
			unique_trees = sorted(unique_trees, key=lambda k: int(k['count']), reverse=True)
	# now let's make sure they're sorted by frequency and rename the categories with a master translate dict
	unique_trees = sorted(unique_trees, key=lambda k: int(k['count']), reverse=True)
	master_translation = {}
	for i,t in enumerate(unique_trees):
		master_translation[t['cat']] = i+1
			#master_translation[cat] = translation[run[0]][cat]

	# write a combined summary file
	with open(args.out + ".summary.txt", "w") as f:
		f.write("tree_cat\ttreetop\tcount\trf\n")
		for tree in unique_trees:
			f.write("{}\t{}\t{}\t{}\n".format(master_translation[tree['cat']], tree['tree'], tree['count'], tree['rf']))

	# and then go through all the long files and translate the tree categories, and combined to a single file
	merged_trees = []
	# loop over the runs
	for run in runs:
		# read the long file
		with open(run[1], "r") as f:
			for line in f.readlines():
				if 'parsed_tree' in line.strip().split():
					headers = line.strip().split()
					continue
				# parse the line
				values = line.strip().split()
				d = dict(zip(headers, values))
				# use the master translation dict to translate the tree category
				if not d['tree_category'] == "NA":
					d['tree_category'] = master_translation[translation[run[0]][d['tree_category']]]
				else:
					d['tree_category'] = "NA"
				merged_trees.append(d)
	# write the output file
	with open(args.out + ".full.txt", "w") as f:
		f.write("\t".join(merged_trees[0].keys()) + "\n")
		for tree in merged_trees:
			f.write("\t".join([str(tree[i]) for i in tree]) + "\n")




