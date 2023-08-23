import functions as fn
import sys
import argparse
import os

parser = argparse.ArgumentParser()

parser.add_argument('-i', '--input', type=str, help='Input sequence.', required=True)
parser.add_argument('-o', '--output-directory', type=str, help="Output directory.")
parser.add_argument('--outformat', type=str, default = "fa")
parser.add_argument('--frame', type=str, help='Frame start of first codon pos.', default = 0)
parser.add_argument('--strict-triplets', action='store_true', help="Add this to require that sequence that's being chopped be divisible by three or throw an error.")

args = parser.parse_args()

# fetch input
msa = fn.readSequenceFile(args.input)

# prep a new msa for each codon
msa_1 = fn.MSA()
msa_2 = fn.MSA()
msa_3 = fn.MSA()

# chop up each sequence in the msa
for seq in msa.sequences.values():
	prefix = ''.join(os.path.basename(args.input).split(".")[0:-1])
	codons = seq.splitCodons(offset = args.frame, not_divisible_by_three=args.strict_triplets)
	#print(codons[0].sequence)
	msa_1.addSample(name = codons[0].sample, seq = codons[0].sequence, meta = codons[0].meta)
	#print(msa_1.sequences[msa_1.samples[0]].sequence)
	msa_2.addSample(name = codons[1].sample, seq = codons[1].sequence, meta = codons[1].meta)
	#print(msa_2.sequences[msa_2.samples[0]].sequence)
	msa_3.addSample(name = codons[2].sample, seq = codons[2].sequence, meta = codons[2].meta)
	fn.writeSequenceFile(msa_1, os.path.join(args.output_directory,prefix + '_codon_1.' + args.outformat))
	fn.writeSequenceFile(msa_2, os.path.join(args.output_directory,prefix + '_codon_2.' + args.outformat))
	fn.writeSequenceFile(msa_3, os.path.join(args.output_directory,prefix + '_codon_3.' + args.outformat))