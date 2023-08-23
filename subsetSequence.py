import functions as fn
import sys
import argparse
import os

parser = argparse.ArgumentParser()

parser.add_argument('-i', '--input', type=str, help='Input sequence.', required=True)
parser.add_argument('-o', '--output', type=str, help="Output file or, if multiple regions are to be extracted, output prefix.")

parser.add_argument('-r', '--region', type=str, help='Region to subset as start-end, one-based inclusive.')

parser.add_argument('--strand', type=str, help="-/+, if -, reverse complement will be output.", default="+")

parser.add_argument('-R', '--regions-file', type=str, help='Regions file with at least two columns, first column should give start coordinate and second end coord, in one-based inclusive format. Third column, if included, should contain strand to extract from (see above), and a fourth column, if included, gives the name to be appended to output file.')

parser.add_argument('--index', type=int, default=1, help="I'm assuming 1-based coordinates given in region or regions-file, can be changed with this option.")

parser.add_argument('-gff', '--gff-file', type=str, help='A gff file from which to extract sequence regions. Specific features to be extracted can be specified with the --extract-types flag.')

parser.add_argument('--extract-types', type=str, help="Comma separated feature types to extract given a gff input file. Leave out to extract all features.")

parser.add_argument('--out-format', type=str, default="phylip")

parser.add_argument('-s', '--sequences', type=str, help="Optional: file listing the sequences to keep, one per line, or comma separated list.")

args = parser.parse_args()

# read input file to msa object
if args.input.endswith((".phy",".phylip",".phy.gz",".phylip.gz")):
	msa = fn.MSA.fromPhylip(args.input)
elif args.input.endswith((".fa",".fasta",".fa.gz",".fasta.gz")):
	msa = fn.MSA.fromFasta(args.input)
else:
	sys.stderr.write("Check input file format, need to be fasta or phylip.")
	sys.exit()

# make a list of the regions to extract
regions = []
if args.region:
	# if a single region is given, add this to list of regiins, along with the output file
	regions.append((int(args.region.split("-")[0] - args.index),int(args.region.split("-")[1]), args.strand, args.output))
else:
	# set suffix based on outformat
	if args.out_format == "phylip":
		suffix = "phy"
	else:
		suffix = "fa"
	if args.regions_file:
		# if a regions file is supplied, go through it and add regions along with an interfix that we'll use to construct output filename
		for line in open(args.regions_file).readlines():
			if len(line.split()) == 4:
				start,end,strand,interfix = line.strip().split()
			elif len(line.split()) == 3:
				start,end,strand = line.strip().split()
				interfix = start + "_" + end
			elif len(line.split()) == 2:
				start,end = line.strip().split()
				interfix = start + "_" + end
				strand = "+"
			regions.append((int(start - args.index),int(end),strand,"{prefix}_{interfix}.{suffix}".format(prefix=args.prefix,interfix=interfix,suffix=suffix)))
	elif args.gff_file:
		gff_all = fn.parseGFF(args.gff_file)
		types = args.extract_types.split(",") if args.extract_types else set([k['feature_type'] for k in gff_all])
		gff_filt = [f for f in gff_all if f['feature_type'] in types]
		for f in gff_filt:
			start = int(f['start'])
			end = int(f['end'])
			strand = f['strand']
			if f['gene_name']:
				filename = "{prefix}_{name}_{start}_{end}.{suffix}".format(prefix=args.output,name=f['gene_name'],start=f['start'],end=f['end'],suffix=suffix)
			else:
				filename = "{prefix}_{start}_{end}.{suffix}".format(prefix=args.output,start=f['start'],end=f['end'],suffix=suffix)
			regions.append((start,end,strand,filename))

# check if samples/sequences to keep were given as an argument, otherwise take all samples from msa object
samples = []
if args.sequences:
	if os.path.isfile(args.sequences):
		for line in open(args.sequences).readlines():
			if line.strip() == "":
				continue
			samples.append(line.strip())
	else:
		samples = [i.strip() for i in args.sequences.split(",")]
	sys.stderr.write("Keeping the following sequences: \n{}\n".format(samples))
else:
	samples = [s for s in msa.samples]

# keep only those samples
msa.subsetSamples(samples)

# if we're extracting regions
if len(regions) > 0:
	sys.stderr.write("Extracting {} regions from input file. \n".format(len(regions)))
	# now just loop through the regions and write them out
	for region in regions:
		# make a copy of complete msa
		msa_copy = msa.copy()
		for seq in msa_copy.sequences.values():
			seq.subsetSequence(region[0],region[1])
			# if minus strand, reverse complement the sequence
			if region[2] == "-" or region[2] == "–":
				seq.reverseComplement()
		# update sequence length
		msa_copy.checkAlnLength()
		if region[3].endswith((".phy",".phylip")):
			msa_copy.writePhylip(region[3])
		elif region[3].endswith((".fa",".fasta")):
			msa_copy.writeFasta(region[3])
else:
	# just write the full msa object
	if args.output.endswith((".phy",".phylip")):
		msa.writePhylip(args.output)
	elif args.output.endswith((".fa",".fasta")):
		msa.writeFasta(args.output)

sys.stderr.write("Done. Output written to {}.\n".format(args.output))