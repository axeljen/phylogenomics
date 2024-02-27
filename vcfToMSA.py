import argparse
import sys
import os
import functions as fn
import pysam
#set up and parse input arguments

parser = argparse.ArgumentParser(description="Fetch gentypes from a vcf and make an alignment.")

parser.add_argument('-v','--vcf', type=str, help='Input vcf')

parser.add_argument('-o','--out', type=str, help='Path to write alignment too, if multiple alignments this will be used as a prefix and region coordinates will be appended to this.')
parser.add_argument('-of','--out-format',choices=['fasta','phylip'], type=str, help="If out is given as prefix, use this to specify the output format. Defaults to fasta.", default='fasta')

parser.add_argument('-r','--region', type=str, help='Optional, region given as chr:start-end', default=None)
parser.add_argument('--strand', choices=['+','-'], help="If region is specified to be on minus strand, the reverse complement will be output. Defaults to '+'.", default='+')

parser.add_argument('-R', '--regions-file', type=str, help="Optional: file with regions to convert to MSA from the input vcf file. Three-four tab separated columns: chrom start end (one-based inclusive). The fourth optional column can be included to suppply '+' or '-' strand. If given, regions on the '-' strand will be converted to reverse complement.")

parser.add_argument('--concat', action="store_true", help="Can be used in conjunction with a regions file, to specify that the supplied regions should be concatenated into a single alignment. Without this flag, regions will be written as separate alignments, using the --out argument as a prefix.")

parser.add_argument('--max-missing', type=float, help="Skip or mask sites with more than <float> fraction of missing genotypes.", default=1)

parser.add_argument('--reference', type=str, help="Indexed reference genome in fasta format can be included to extrapolate missing regions in the vcf file, either with 'N's (default) or other missing character of choice or with the reference base (change this with --mask-with flag).")

parser.add_argument('--mask-with', type=str, help="Can be used in conjuction with --reference, to choose what letter to mask sites missing from the vcf file with. Choose 'ref' to use the reference base, or any single letter character. Defaults to 'N'.", default="N")

parser.add_argument('--bypass-filters', action="store_true", help="Set this to ignore any filters in the vcf file. Without it, only sites with 'PASS' in the filter column will be considered.")

parser.add_argument('-s','--samples', type=str, help="File with one sample per line, for all samples to include in the extraction.", default=None)

parser.add_argument('--haploidize', choices=['IUPAC','iupac','random','dont'], help="Haploidize diploid records and change heterozygote genotypes to the IUPAC ambiguity code or sample a random allele. Defaults to iupac. If dont is given, each sample will be output as two haploid sequences.",default='IUPAC')

parser.add_argument('--extrapolate-missing-sites', choices=['ref','N','-'], help="Character to use when extending the alignment to span over sites that are missing in the vcf file. Used only in conjuction with --reference. Defaults to 'N'.", default="N")

# add options to filter on sample DP
parser.add_argument('--min-sample-dp', type=int, help="Minimum sample DP to include a site.", default=None)
parser.add_argument('--max-sample-dp', type=int, help="Maximum sample DP to include a site.", default=None)

args = parser.parse_args()


# parse input vcf to pysam variantfile
sys.stderr.write("Starting to convert {} to alignment.\n\n".format(args.vcf))
vcf = pysam.VariantFile(args.vcf)

# check if region or regionfile is given, and parse into a list
regions = []
if args.region:
	regions.append((args.region.split(":")[0],args.region.split(":")[1].split("-")[0], args.region.split(":")[1].split("-")[1],args.strand))
elif args.regions_file:
	for line in open(args.regions_file).readlines():
		if line.strip() == "":
			continue
		try:
			regions.append((line.split("\t")[0].strip(),int(line.split("\t")[1].strip()),int(line.split("\t")[2].strip()),line.split("\t")[3].strip()))
		except:
			regions.append((line.split("\t")[0].strip(),int(line.split("\t")[1].strip()),int(line.split("\t")[2].strip()),args.strand))


# if a sample file was given, parse it, otherwise grab all the samples from vcf header
samples = []
if args.samples:
	with open(args.samples) as f:
		for line in f.readlines():
			if line == "":
				continue
			samples.append(line.strip())
else:
	samples = list(vcf.header.samples)

sys.stderr.write("\nIncluding the following {} samples in output alignment:\n".format(len(samples)))
sys.stderr.write(','.join(samples))

# check if we should haploidize sequences
if args.haploidize == 'dont':
	haploidize = False
	heterozygotes = None
	sys.stderr.write("\nTwo sequences will be written per individual, as --haploidize was set to 'dont'.\n")
else:
	haploidize = True
	heterozygotes = args.haploidize
	sys.stderr.write("\nHeterozygote genotypes will be haploidized using the {} option.\n".format(args.haploidize))


# if no region is given, we will parse the entire vcf file
if len(regions) == 0:
	sys.stderr.write("No region was specified, will parse all records in the vcf.".format(args.vcf))
	msa = fn.MSA.fromVcf(vcf,None,None,None,samples,haploidize=haploidize,heterozygotes=heterozygotes, padding=args.mask_with, skip_indels = True, index = 1)
	# should strandedness be minus, change to reverse complement
	if args.strand == '-':
		msa.reverseComplement()
	# if max missing is below one, filter it
	if args.max_missing < 1:
		sys.stderr.write("\nFiltering out sites with a higher fraction than {} of missing data.\n".format(args.max_missing))
		msa.filterMSA(args.max_missing)
	# write it to output
	if args.out.endswith(("phy","phylip","phy.gz","phylip.gz")):
		msa.writePhylip(args.out)
	elif args.out.endswith(("fa",".fasta",".fa.gz",".fasta.gz")):
		msa.writeFasta(args.out)
	else:
		sys.stderr.write('Cannot determine output format, please check name of output file and try again.')
		sys.exit(1)
	sys.stderr.write("\nDone. Output alignment ({}) written here:\n{}\n".format(msa.length,args.out))

elif len(regions) == 1:
	# if a single region is given, parse much like the above after setting the region
	chrom,start,end,strand = regions[0][0],int(regions[0][1]),int(regions[0][2]),regions[0][3]
	sys.stderr.write("\nParsing the region {}:{}-{} to an alignment.\n".format(chrom,start,end))
	msa = fn.MSA.fromVcf(vcf,chrom,start,end,samples,haploidize=haploidize,heterozygotes=heterozygotes, padding=args.mask_with, skip_indels = True, index = 1,
	minsampleDP=args.min_sample_dp, maxsampleDP = args.max_sample_dp)
	# should strandedness be minus, change to reverse complement
	if strand == '-':
		msa.reverseComplement()
	# if max missing is below one, filter it
	if args.max_missing < 1:
		sys.stderr.write("\nFiltering out sites with a higher fraction than {} of missing data.\n".format(args.max_missing))
		msa.filterMSA(args.max_missing)
	if args.reference:
		# if a ref genome is given, extrapolate the alignment on top of this now. only works when a region is specified for now!
		ref = pysam.FastaFile(args.reference)
		refseq = ref.fetch(chrom,start -1,end)
		msa.addRefGenome(refseq,start,end,refname='reference', mask_uncalled = args.extrapolate_missing_sites, strand=args.strand)
	# write it to output
	if args.out.endswith(("phy","phylip","phy.gz","phylip.gz")):
		msa.writePhylip(args.out)
	elif args.out.endswith(("fa",".fasta",".fa.gz",".fasta.gz")):
		msa.writeFasta(args.out)
	else:
		sys.stderr.write('Cannot determine output format, please check name of output file and try again.')
		sys.exit(1)
	sys.stderr.write("\nDone. Output alignment ({} bp) written here:\n{}\n".format(msa.length,args.out))

else:
	# otherwise this means that we're dealing with multiple regions, start by storing them all in a list
	msa_list = []
	sys.stderr.write("\nFound {} regions in the supplied regions file. Converting those to alignments.\n.".format(len(regions)))
	for region in regions:
		#print(samples)
		chrom,start,end,strand = region[0],int(region[1]),int(region[2]),region[3]
		sys.stderr.write("\nNow parsing the region {}:{}-{}.".format(chrom,start,end))
		msa = fn.MSA.fromVcf(vcf,chrom,start,end,samples,haploidize=haploidize,heterozygotes=heterozygotes, padding=args.mask_with, skip_indels = True, index = 1,
		minsampleDP=args.min_sample_dp, maxsampleDP=args.max_sample_dp)
		# should strandedness be minus, change to reverse complement
		if strand == '-':
			msa.reverseComplement()
		# if max missing is below one, filter it
		if args.max_missing < 1:
			sys.stderr.write("\nFiltering out sites with a higher fraction than {} of missing data.\n".format(args.max_missing))
			msa.filterMSA(args.max_missing)
		if args.reference:
			# if a ref genome is given, extrapolate the alignment on top of this now. only works when a region is specified for now!
			ref = pysam.FastaFile(args.reference)
			print(chrom,start,end)
			refseq = ref.fetch(chrom,start -1,end)
			msa.addRefGenome(refseq, start, end, refname='reference', mask_uncalled = args.extrapolate_missing_sites, strand = args.strand)
		# add to the list
		msa_list.append(msa)
	#print(msa_list)
	if args.concat:
		# concatenate alignments to a single one
		msa = msa_list[0]
		for m in msa_list[1:]:
			msa.concatenateMSA(m)
		# write 
		if args.out.endswith(("phy","phylip","phy.gz","phylip.gz")):
			msa.writePhylip(args.out)
		elif args.out.endswith(("fa",".fasta",".fa.gz",".fasta.gz")):
			msa.writeFasta(args.out)
		else:
			sys.stderr.write('Cannot determine output format, please check name of output file and try again.')
			sys.exit(1)
		sys.stderr.write("\nDone. Output alignment ({} bp) written here:\n{}\n".format(msa.length,args.out))
	else:
		# write each msa to individual files
		if args.out_format == "fasta":
			suffix = "fa"
		else:
			suffix = "phy"
		for m in msa_list:
			print(m.start)
			outfile = args.out + "{}_{}_{}.{}".format(m.chrom,m.start,m.end,suffix)
			if suffix == 'fa':
				m.writeFasta(outfile)
			else:
				m.writePhylip(outfile)
		sys.stderr.write("\nDone. Output alignments written with prefix:\n{}\n".format(args.out))


	

