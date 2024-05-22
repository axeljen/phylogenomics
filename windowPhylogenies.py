#!/usr/bin/python

## Axel Jensen, 2022
## This python script was originally adapted from Simon Martins phyml_sliding_windows.py
## https://github.com/simonhmartin/genomics_general/tree/master/phylo
## but works directly on a vcf file, converts it into window-alignments and run iqtree or phyml phylogenetic estimates

import sys
import os
import argparse
import gzip
from multiprocessing import Process
from multiprocessing import SimpleQueue
from threading import Thread
from time import sleep
import tempfile
#from zipfile import ZipFile
from pysam import VariantFile
import random
#from vcf.functions import getChromLengths, vcfToSeqDict, generateWindows, filterSeqDict
import functions as fn
from datetime import datetime
import glob

def IQTreeCommand(iqtree_path, model, bootstraps, alignment, outgroup, prefix, logfile, bootstrap_flag):
	if outgroup:
		outgroupArg = " -o " + outgroup
	else:
		outgroupArg = ""
	if bootstraps:
		bootstraps_command = " -B " + str(bootstraps)
	else:
		bootstraps_command = ""
	iqtree_command = iqtree_path + "2" + " -s " + alignment + " -pre " + prefix + outgroupArg + bootstraps_command + " -m " + model + " -nt 1 " + ">>" + logfile + " 2>&1 "
	os.system(iqtree_command)

def PHYMLCommand(phyml_path, model, opt, bootstraps, alignment, directory, prefix, logfile):
	if opt:
		opt_command = " -o " + opt
	else:
		opt_command = " -o n "
	if bootstraps:
		bootstraps_command = " -b " + str(bootstraps)
	else:
		bootstraps_command = " -b 0 "
	phyml_command = phyml_path + " -i " + alignment + bootstraps_command + opt_command + " --model " + model + ">>" + logfile + " 2>&1 "
	os.system(phyml_command)

def parseIQTREEOutput(directory, prefix, windowCoords):
	try:
		tree_logfile = open(os.path.join(prefix + ".log"))
		for line in tree_logfile.readlines():
			if "Alignment has" in line:
				info_line = line
		sequences = info_line.split()[2]
		alignment_sites = info_line.split()[5]
		treefile = open(os.path.join(prefix + ".treefile"))
		tree = treefile.readline().rstrip('\n')
		results = True
	except:
		if args.verbosity == "warn":
			sys.stderr.write("Warning: no tree was constructed for window {chrom}:{start}-{end}.\n".format(chrom=windowCoords['chrom'], start=windowCoords['start'], end=windowCoords['end']))
		tree = "NA"
		results = False
		alignment_sites = "NA"
	return tree,alignment_sites

def parsePHYMLOutput(prefix, windowCoords):
	try:
		with open(prefix + "_phyml_tree.txt", "rt") as treeFile: 
			tree = treeFile.readline().strip()
			results = True
	except:
		try:
			with open(prefix + "_phyml_tree", "rt") as treeFile:
				tree = treeFile.readline().strip()
				print("couldn't find {}.phy_phyml_tree.txt".format(prefix))
				results = True
		except:
			sys.stderr.write("Tree not found at " + prefix + "_phyml_tree.txt\n")
			print("couldn't find {}_phyml_tree".format(prefix))
			tree = "NA"
			results = False
			if args.verbosity == "warn":
				sys.stderr.write("Warning: no tree was constructed for window {chrom}:{start}-{end}.\n".format(chrom=windowCoords['chrom'], start=windowCoords['start'], end=windowCoords['end']))
	try:
		with open(prefix + "_phyml_stats.txt", "rt") as statsFile:
			stats = statsFile.read().split()
			lnL = stats[stats.index("Log-likelihood:")+1]
	except:
		try:
			with open(prefix + "_phyml_stats", "rt") as statsFile:
				stats = statsFile.read().split()
				lnL = stats[stats.index("Log-likelihood:")+1]
		except:
			lnL = "NA"
			stats = "NA"
	return tree,stats,lnL

def run_window_phylogeny(windowQueue, resultQueue, input_vcf, software, model, bootstraps, opt, outgroup, directory, output_temp, logfile, bootstrap_flag, max_missing=None, phyml_path=None, iqtree_path=None, test=False, ploidy=2):
	global timePassed
	global windowsTotal
	global windowsQueued
	while True:
		if _FINISH:
			break
		windowsQueued,window,samples,haploidize,heterozygotes = windowQueue.get()
		vcf = VariantFile(input_vcf)
		msa = fn.MSA.fromVcf(vcf, window['chrom'], window['start'], window['end'], samples=samples, haploidize = haploidize, heterozygotes = heterozygotes, index=0, ploidy=int(ploidy))
		if max_missing:
			msa.filterMSA(max_missing)
		seqLen = msa.length
		nseq = len(msa.samples)
		prefix = os.path.join(os.path.basename(directory) + "_" + str(window['start']) + "_" + str(window['end']))
		tempAlignment = tempfile.NamedTemporaryFile(mode="w", prefix=prefix, suffix=".phy", dir = directory, delete=False)
		msa.writePhylip(tempAlignment.name)
		outprefix = os.path.join(directory, prefix)
		tempAlignment.close()
		if software == "phyml":
			PHYMLCommand(phyml_path, model, opt, bootstraps, tempAlignment.name, directory, outprefix, logfile)
			tree,stats,lnL = parsePHYMLOutput(tempAlignment.name, window)
			treeDict = {'window_number':str(window['window_number']), 'chrom': str(window['chrom']), 'start':str(int(window['start']) + 1), 'stop':str(window['end']), 'lnL': str(lnL), 'samples': str(nseq), 'sites':str(seqLen), 'tree':tree}
		elif software == "iqtree":
			IQTreeCommand(iqtree_path, model, bootstraps, tempAlignment.name, outgroup, outprefix, logfile, bootstrap_flag)
			tree,nsites = parseIQTREEOutput(directory, outprefix, window)
			treeDict = {'window_number':str(window['window_number']), 'chrom': str(window['chrom']), 'start':str(int(window['start']) + 1), 'stop':str(window['end']), 'sequences':str(nseq), 'sites':str(nsites), 'tree':tree}
		open(output_temp, "a").write(str('\t'.join(treeDict.values())) + "\n")
		#global resultsReceived
		# remove all files from this run
		if not test:
			filelist = glob.glob(os.path.join(directory, prefix) + "*")
			for f in filelist:
				os.remove(f)

##indefinite loop to write some stats every tenth second as the program runs:
def checkStats(output_temp):
	global resultsReceived
	#resultsReceived = len(results)
	time_passed = 0
	while True:
		sleep(int(args.check_time))
		try:
			finished = open(output_temp).read().split("\n")
			resultsReceived = len(finished) - 1
		except:
			resultsReceived = 0
		time_passed = time_passed + 30
		percentage = (resultsReceived / windowsTotal) * 100
		sys.stderr.write("{time}: Received results for {received}/{total} windows ({perc} %)\n".format(time = datetime.now().strftime("%H:%M:%S"), received = resultsReceived, total = windowsTotal, perc=round(percentage, 2)))
		#print("Finished with " + str(resultsReceived) + "/" + str(windowsTotal) + " windows (" + str(percentage) + "%) in " + str(time_passed) + " seconds.")

#just a function to check whether to use gzip or not when opening input
def openFile(filepath, mode='r'):
	if filepath.endswith(".gz"):
		return gzip.open(filepath, 'rt')
	else:
		return open(filepath, mode)
    
##########################################################################################################################

#if __name__ == "__main__":

parser = argparse.ArgumentParser()

#arguments for in and output:
parser.add_argument("-i", "--input", help="Input vcf file (.vcf or .vcf.gz).", required=True, metavar="genotypes.vcf | genotypes.vcf.gz")
parser.add_argument("-o", "--output", help="Output file.", required=True, metavar="output.tsv")

#arguments for subsetting and splitting alignment:
parser.add_argument("--sequences", help="Headers of sequences to include, either as a file with one header per line, or in the command separated by commas.", metavar="<seq1,seq2,seq3>|seqs.file")
parser.add_argument("--max-missing", type=float, help="Maximum fraction of missing bases ('N' or '-') allowed for a posititon to be kept.", metavar="<int:nBases>")

parser.add_argument("-r","--region", type=str, help="Restrict analyses to this region. If given in combination with window size, windows will be generated within this region. Otherwise, the region will be parsed as a whole and a tree will be output.")

parser.add_argument("-w", "--window-size", help="Window size in bp.", type=int, metavar="bases")
parser.add_argument("--step-size", help="Step size in bp.", type=int, metavar="bases")

parser.add_argument("--ploidy", help="Ploidy of samples in the vcf file (1 or 2).", default=2)

parser.add_argument("--haploidize", help="Combine/subsample diploid loci to one haploid sequence", action="store_true")
parser.add_argument("--handle-heterozygotes", help="Either randomly sample an allele (random), use IUPAC codes (IUPAC), always take the first (first) or second (second) allele.", type=str, default="IUPAC") # changed this to heterozygotes only, keeping to avoid breaking
parser.add_argument("--heterozygotes", help="Either randomly sample an allele (random), use IUPAC codes (IUPAC), always take the first (first) or second (second) allele.", type=str, default="IUPAC")

#arguments for running iqtree:
parser.add_argument("--iqtree-path", help="If iqtree executable isn't in PATH, you can specify the absolute path here", metavar="path/to/iqtree", default="iqtree")
parser.add_argument("--model", help="iqtree substitution model", default="GTR")
parser.add_argument("--outgroup", help="Outgroup(s) for the tree construction, if desired", default=None, metavar="outgroup_seq1,outgroup_seq2")
parser.add_argument("--bootstraps", help="Number of ultrarapid bootstraps to run in iqtree, defaults to none.", type=int, metavar="number_of_bootstraps", default = 0)

#number of threads we're running on:
parser.add_argument("-T", "--threads", help="Number of threads for parallel processing", type=int, default=1, metavar="threads")

#can give a file with predefined coordinates/windows to build trees on
parser.add_argument("--windows-file", help="File with predefined windows/regions to analyze. 'chrom\tstart\tstop'. Only works when input is a vcf as of now.")

#optionally run phyml instead of iqtree, which is still the default for now
parser.add_argument("--phyml", help="Specify to run phyml instead of iqtree.", action="store_true")
parser.add_argument("--phyml-path", type=str, help="Path to phyml executable, if not in path", default="phyml")
#a few options specific to phyml
parser.add_argument("--phyml-model", type=str, default="HKY85", help="Substitution model for running phyml.")
parser.add_argument("--optimize", type=str, default="n", help="Parameter optimization, see phyml manual for options.")

# an option to run analyses in a separate directory from wd
parser.add_argument("--analysis-directory", type=str, default=".")

parser.add_argument("--test", help="skip deleting some tmp files etc., for debugging.", action="store_true")
parser.add_argument("--verbosity", choices = ('warn','err'), help="A bit more verbose output for debugging.", default='err')
parser.add_argument("--check-time", help="Sleeptime when checking results, used for testing.", default=5)
args = parser.parse_args()

# print date etc
sys.stderr.write("#" * 100 + "\n\n" + "{datetime}: starting window phylogenies analyses.\n\n".format(datetime = datetime.now().strftime("%Y-%m-%d, %H:%M:%S")) + "#" * 100 + "\n\n")

### open files
vcf = VariantFile(args.input)
output_file = open(args.output, "w")

#set step length same as window if not given:
if args.window_size:
	window_size = args.window_size
	if args.step_size:
		step_size=args.step_size
	else:
		step_size=int(args.window_size)
else:
	window_size = None
sys.stderr.write("Parameters: \n")
sys.stderr.write("Input vcf read from: {inputfile}, writing final window output to {outputfile}.\n".format(inputfile=args.input, outputfile=args.output))
#if regions file is provided, parse this one to a dictionary
if args.windows_file:
	regions = []
	with open(args.windows_file) as r:
		for number,line in enumerate(r.readlines()):
			if not line == "":
				regions.append({'window_number':number, 'chrom':line.split()[0], 'start':int(line.split()[1]), 'end':int(line.split()[2])})
	sys.stderr.write("Intervals were supplied from file, using these intervals.\n")			
# otherwise check if a region is given
elif args.region:
	chrom = args.region.split(":")[0]
	start,end = [int(i) for i in args.region.split(":")[1].split("-")]
	sys.stderr.write("Region to analyze from command argument: \n{chrom}\n{start}\n{end}\n".format(chrom=chrom,start=start,end=end))
	# check if a window size was given:
	if window_size:
		# generate windows within this region
		sys.stderr.write("Chopping the specified region up into windows:\nWindow size: {}, step size:{}\n".format(window_size,step_size))
		regions = fn.generateWindows({chrom: (start -1,end)}, window_size,step_size)
		#print(regions)
		#sys.exit()
	else:
		# make a single region of this interval for downstream analysis
		sys.stderr.write("Constructing a phylogeny for the entire specified region.\n")
		regions = [{'window_number': 1, 'chrom': chrom, 'start': start -1, 'end': end}]
else:
	sys.stderr.write("No interval file was supplied, will generate windows from the contigs in the vcf input.\n")
	sys.stderr.write("Window size: {wsize}, step size: {stepsize}\n".format(wsize=str(args.window_size), stepsize=str(step_size)))
	chroms = fn.getChromLengths(vcf)
	regions = fn.generateWindows(chroms, int(args.window_size), step_size)
windowsTotal = len(regions)
if not args.bootstraps:
	bootstraps = None
else:
	bootstraps = args.bootstraps
# set bootstrap flag to "-bb" as default
bootstrap_flag = "-bb" 

#create a temporary directory to work in
analysis_directory = args.analysis_directory
tmpdir = tempfile.mkdtemp(suffix=None, prefix=None, dir=analysis_directory)

###create tempfile to store unsorted trees in:
output_temp = os.path.join(tmpdir, os.path.basename(args.input) + "_out_temp")

#write headers to output file and check that the software needed is there
if not args.phyml:
	if args.iqtree_path:
		iqtree = args.iqtree_path
	else:
		iqtree = "iqtree"
	t = os.system(iqtree + " >> /dev/null 2>&1")
	if t > 0:
		t = os.system("iqtree2  >> /dev/null 2>&1")
		if t > 0:
			sys.stderr.write("Error: Please make sure you have iqtree or iqtree2 in your PATH, or give full path to executable. Will exit now.\n")
			sys.exit(1)
		else:
			sys.stderr.write("Running iqtree2.\n")
			iqtree = "iqtree2"
			bootstrap_flag = "-B"
	logfile = args.output[:-4] + "_iqtree.log.txt"
	sys.stderr.write("Running phylogenies with IQTree.\n")
	output_file.write("window\tchrom\tstart\tend\tsequences\tsites\ttree\n")
else:
	# check so that phyml is in path
	if args.phyml_path:
		phyml = args.phyml_path
	else:
		phyml = "phyml"
	#t = os.system("phyml >> /dev/null 2>&1")
	#if t > 0:
	#	sys.stderr.write("Error: Please make sure you have phyml in your PATH, or give full path to executable. Will exit now.\n")
	#	sys.exit(1)
	logfile = args.output[:-4] + "_phyml.log.txt"
	sys.stderr.write("Running phylogenies with PHYML.\n")
	output_file.write("window\tchrom\tstart\tend\tLnl\tsequences\tsites\ttree\n")

if os.path.isfile(logfile):
	sys.stderr.write("Logfile with the given prefix already exist, this will be overwritten.\n")
	os.remove(logfile)

#If sequences to subsample are given to input, extract those:
if args.sequences:
	if os.path.isfile(args.sequences):
		headers_to_keep = open(args.sequences).read().split("\n")
		headers_to_keep = [sample for sample in headers_to_keep if not sample == ""]
	else:
		headers_to_keep = args.sequences
	#sys.stderr.write("Keeping the specified samples for analysis:\n " + ','.join([sample for sample in headers_to_keep if not sample == ""]) + "\n")
else:
	headers_to_keep = "all"
	#sys.stderr.write("Keeping all samples in the vcf file for analyses.\n") 
#outgroup to none if not given
if not args.outgroup:
	outgroup=None
	sys.stderr.write("No outgroup was specified.\n") 
else:
	outgroup=args.outgroup
	sys.stderr.write("Using the following sample(s) as outgroup: {}.\n".format(outgroup)) 

# are we running phyml or iqtree?
if args.phyml:
	model = args.phyml_model
	software = "phyml"
	opt = args.optimize
	sys.stderr.write("Will run PHYML using the {phymlmodel} model, and optimization set to {optim}\n".format(phymlmodel = model, optim = opt))
else:
	software = "iqtree"
	model = args.model
	opt = None
	sys.stderr.write("Will run IQTree using the {iqtreemodel} model.\n".format(iqtreemodel = model))
sys.stderr.write("Bootstraps set to {}\n".format(bootstraps))

# should we haploidize the sequences?
if args.haploidize:
	haploidize = True
	heterozygotes = args.heterozygotes
	sys.stderr.write("Heterozygous positions will be haploidized, and replacing heterozygous using setting: {}.\n".format(heterozygotes))
else:
	haploidize = False
	heterozygotes = None
	#heterozygotes = args.heterozygotes
	sys.stderr.write("Will generate a diploid sequence for each samples, suffixed __1 and __2. Phasing will not be taken into account here.\n")

if not args.max_missing:
	max_missing = None
else:
	max_missing = args.max_missing

#### Global variables for keeping track of stats:
global _FINISH
_FINISH = False
global windowsQueued
windowsQueued = 0
resultsReceived = 0
resultsWritten = 0
resultsHandled = 0
windowsTotal = windowsTotal
global analysisStarted
### Create queues
#queue for holding alignment windows:
windowQueue = SimpleQueue()
#queue for holding results:
resultQueue = SimpleQueue()
#que for holding sorted results:
writeQueue = SimpleQueue()


#list to put all processes in:
workerprocesses = []
if args.test:
	test = True
else:
	test = False

for thread in range(args.threads):
	#vcfToSeqDict(vcf, window['chrom'], haploidize = haploidize, heterozygotes = heterozygotes, samples=samples, start=window['start'], end=window['end'])
	worker =Process(target=run_window_phylogeny, args=(windowQueue, resultQueue, args.input, software, model, 
	bootstraps, opt, outgroup, tmpdir, output_temp, logfile, bootstrap_flag, max_missing, args.phyml_path, 
	args.iqtree_path, test, args.ploidy,))
	worker.daemon = True
	worker.start()
analysisStarted = True

stats = Thread(target=checkStats, args=(output_temp,))
stats.daemon = True
stats.start()
sys.stderr.write("\n" + "#" * 100 + "\n\n")
sys.stderr.write("Starting analysis on " + str(args.threads) + " processes in parallell.\n")

# start generating the window alignments and trees
if headers_to_keep == "all":
	samples = list(vcf.header.samples)
else:
	samples = headers_to_keep

# loop through the windows and feed to the queue
for window in regions:
	while windowsQueued - resultsReceived >= 50:
		sleep(int(args.check_time))
	windowQueue.put((windowsQueued,window,samples,haploidize,heterozygotes))
	windowsQueued += 1
analysisStarted = True
while resultsReceived < windowsQueued:
	sleep(int(args.check_time))
#	print("Received: " + str(resultsReceived) + ", queued: " + str(windowsQueued))
print("Done.")

_FINISH=True

sys.stderr.write("\n\nDone with analysis, sorting and writing final output.\n")
#grab unsorted outtemp, loop through it to sort outfiles
results_unsorted = []
with open(output_temp) as file:
	for line in file.readlines():
		results_unsorted.append(line)
for i in range(1, len(results_unsorted) + 1):
	for result in results_unsorted:
		if result.split()[0] == str(i):
			output_file.write(result)
output_file.close()
if not args.test:
	for i in os.listdir(tmpdir):
		os.remove(os.path.join(tmpdir, i))
	os.rmdir(tmpdir)
#	os.remove(output_temp)
sys.stderr.write("All done.\n")
sys.exit(0)

