import functions as fn
import sys
import argparse
import os

parser = argparse.ArgumentParser()

parser.add_argument('-i', '--input', type=str, help='Input sequence.', required=True)
parser.add_argument('-o', '--output', type=str, help="Output file.")

parser.add_argument('--max-missing-site', type=float, help="Maximum proportion of missing data allowed at a site. This step takes place after filtering out individuals with too much missing data.", default=1.0)
parser.add_argument('--max-missing-ind', type=float, help="Maximum proportion of missing data allowed for a sequence.", default=1.0)
parser.add_argument('--codon-aware', action="store_true", help="If specified the alignment will be expected to be inframe codon alignment, and filtration will be done on triplets.")
parser.add_argument('--trim-trailing-stop', action="store_true", help="If specified, stop codons will be trimmed from the end of the sequences. Only used together with --codon-aware.")
parser.add_argument('--remove-on-internal-stop', action="store_true", help="If specified, sequences with internal stop codons will be removed.")

args = parser.parse_args()

# read input file to msa object
if args.input.endswith((".phy",".phylip",".phy.gz",".phylip.gz")):
	msa = fn.MSA.fromPhylip(args.input)
elif args.input.endswith((".fa",".fasta",".fa.gz",".fasta.gz")):
	msa = fn.MSA.fromFasta(args.input)
else:
	sys.stderr.write("Check input file format, need to be fasta or phylip.")
	sys.exit()

# check amount of missingness per individual first
samples = []
samples_with_internal_stops = []
samples_with_high_missingness = []
for s in msa.samples:
	if msa.sequences[s].missingness() <= args.max_missing_ind:
		# if we're removing sequences with internal stop codons, do this now
		if args.remove_on_internal_stop:
			if msa.sequences[s].containsInternalStopCodon():
				sys.stderr.write("Removing {} due to internal stop codon.\n".format(s))
				samples_with_internal_stops.append(s)
				continue
		samples.append(s)
	else:
		samples_with_high_missingness.append(s)
		sys.stderr.write("Removing {} due to too much missing data.\n".format(s))

# subset the msa object to only include the samples that passed the filter
msa.subsetSamples(samples)
# filter alignment
if args.codon_aware:
	if args.trim_trailing_stop:
		msa.removeTrailingStopCodon()
	msa.codonAwareFiltration(max_missingness=args.max_missing_site, break_on_stops = False)
else:
	msa.filterMSA(args.max_missing_site)

# write the output
if args.output.endswith((".phy",".phylip")):
	msa.writePhylip(args.output)
elif args.output.endswith((".fa",".fasta")):
	msa.writeFasta(args.output)

sys.stderr.write("Done. Output written to {}.\n".format(args.output))