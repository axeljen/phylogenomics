import functions as fn
import sys
import argparse
import os

parser = argparse.ArgumentParser()

parser.add_argument('-i', '--input', type=str, help='Input alignment to extract partitions from.', required=True)
parser.add_argument('-o', '--output', type=str, help="Output prefix.")

parser.add_argument('--partition-file', type=str, help="File specifying the coordinates of the different partitions (1-based inclusive).")

args = parser.parse_args()

# fetch the input alignment
msa = fn.readSequenceFile(args.input)

# fetch the coordinates of each partition from the partitionsfile
partcords = {}
with open(args.partition_file) as f:
	for i,line in enumerate(f.readlines()):
		partcords['partition_{}'.format(i +1)] = [(int(c.strip().split("-")[0]),int(c.strip().split("-")[1])) for c in line.strip().split("=")[1].strip().split(",")]
# make a subalignment for each partition and concatenate the coordinates and write out
for p in partcords.keys():
	m = msa.copy()
	m.subset(partcords[p][0][0], partcords[p][0][1], index = 1)
	rest = []
	for c in partcords[p][1:]:
		rest.append(msa.copy())
		rest[-1].subset(c[0], c[1], index = 1)
	for i in rest:
		m.concatenateMSA(i)
	m.writePhylip(args.output + "_{}.phy".format(p))
