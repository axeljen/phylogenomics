import functions as fn
import sys
import argparse
import os
from pysam import VariantFile

parser = argparse.ArgumentParser()

parser.add_argument('-i', '--input', type=str, help='Input vcf file.', required=True)
parser.add_argument('-o', '--output', type=str, help="Output vcf file.", required=True)

parser.add_argument('-r', '--region', type=str, help='Region to subset as chrom:start-end, one-based inclusive.')

parser.add_argument('-P1', type=str, help="Pop 1 of two to look for fixed differences between", required=True)
parser.add_argument('-P2', type=str, help="Pop 2 of two to look for fixed differences between", required=True)

parser.add_argument('-p', '--popfile', type=str, help="Popfile assigning samples to the two populations to look for fixed differences between.", required=True)

parser.add_argument('-m','--max-missing',type=float,help="Max fraction of missing genotypes allowed.", default=1)

args = parser.parse_args()

# parse popfile
pops = fn.parsePopfile(args.popfile)
# check so that both P1 and P2 are in it and have samples assigned
if not args.P1 in pops.keys() or not args.P2 in pops.keys():
	sys.stderr.write("Please make sure that the two pops are present with assigned samples in the popfile.")
	sys.exit(1)

# write some info on what we're doing
P1_samples = pops[args.P1]
P2_samples = pops[args.P2]
sys.stderr.write("Will look for sites with differentially fixed alleles between specified populations:\n")
sys.stderr.write("{pop}:\n{samples}\n".format(pop=args.P1,samples=','.join(pops[args.P1])))
sys.stderr.write("{pop}:\n{samples}\n".format(pop=args.P2,samples=','.join(pops[args.P2])))
sys.stderr.write("\nMax missingness set to {}.\n\nOutput vcf will be written to:\n{}\n".format(args.max_missing,args.output))

# check if region was given, else set to none
chrom = None
start = None
end = None
if args.region:
	try:
		chrom = args.region.split(":")[0]
		start = int(args.region.split(":")[1].split("-")[0]) - 1 #remove one to convert to zero-indexed
		end = int(args.region.split(":")[1].split("-")[1])
		sys.stderr.write("Will parse region: {}:{}-{}".format(chrom,start+1,end))
	except:
		sys.stderr.write("Invalid region format, please check.\n")
		sys.exit(1)

# loop through the vcf records and check them
vcf = VariantFile(args.input)
vcf_out = VariantFile(args.output,"w",header=vcf.header)
records_written = 0
for rec in vcf.fetch(chrom,start,end):
	p1_alleles = []
	p1_missing = 0
	p2_alleles = []
	p2_missing = 0
	for sample in pops[args.P1]:
		gt = rec.samples[sample]['GT']
		if None in gt:
			p1_missing += 1
			continue
		p1_alleles.append(gt[0])
		p1_alleles.append(gt[1])
	for sample in pops[args.P2]:
		gt = rec.samples[sample]['GT']
		if None in gt:
			p2_missing += 1
			continue
		p2_alleles.append(gt[0])
		p2_alleles.append(gt[1])
	# check missingness
	if (p1_missing / len(pops[args.P1])) > args.max_missing:
		continue
	if (p2_missing / len(pops[args.P2])) > args.max_missing:
		continue
	# check if fixed
	p1_set = set(p1_alleles)
	p2_set = set(p2_alleles)
	if len(p1_set) > 1 or len(p2_set) > 1:
		continue
	elif len(list(p1_set)) == 1 and len(list(p2_set)) == 1:
		if list(p1_set)[0] != list(p2_set)[0]:
			vcf_out.write(rec)
			records_written += 1

# write when finished
sys.stderr.write("All done. Found {} sites with fixed differences.\n".format(records_written))