import functions as fn
import sys
import argparse
import os

parser = argparse.ArgumentParser()

parser.add_argument('-i', '--input', type=str, help='Input files to concatenate, either as file with one per line, or as a comma separated list.', required=True)
parser.add_argument('-o', '--output', type=str, help="Output file. If a partitions file is to be written this will be using the same file path only different suffix.")

parser.add_argument('--raxml-partitions', action = "store_true", help="This option will write file specifying the coordinates of each partitions, ready to use with raxml.")

parser.add_argument('--upper-case', action = "store_true", help="Will force nucleotides to be written out in uppercase letters.")

args = parser.parse_args()

# fetch all the paths into a list
seqpaths = []

if os.path.isfile(args.input):
	for line in open(args.input).readlines():
		if line.strip() == "":
			continue
		else:
			seqpaths.append(line.strip())
else:
	seqpaths = [i.strip() for i in args.input.split(",")]
# parse the files into a list of msa objects
seqlist = []

for sp in seqpaths:
	seqlist.append(fn.readSequenceFile(sp))
# use the first in the file as master msa
msa = seqlist[0]
msa.upperCase()

# loop through and concatenate the rest
for m in seqlist[1:]:
	m.upperCase()
	msa.concatenateMSA(m)


# write the output msa
fn.writeSequenceFile(msa,args.output)

# and if partitions are to be written, write those too
if args.raxml_partitions:
	outfile = os.path.join(os.path.dirname(args.output), ''.join(os.path.basename(args.output).split(".")[0:-1]) + "_raxml_partitions.txt")
	msa.writePartitions(outfile)