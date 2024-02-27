import pysam
import sys
import os
import datetime as dt
# vcf file as first argument
vcf = pysam.VariantFile(sys.argv[1])
# optionally include a second argument for an outfile prefix
if len(sys.argv) > 2:
	prefix = sys.argv[2]
else:
	# otherwise we'll put the output in current directory, named after the input file
	prefix = os.path.splitext(os.path.basename(sys.argv[1]))[0]
	# remove the .vcf extension if it's there
	if prefix.endswith(".vcf"):
		prefix = prefix[:-4]
	# and gz if it's there
	if prefix.endswith(".vcf.gz"):
		prefix = prefix[:-7]

# function to check if the variant is a SNP
def is_snp(variant):
	# if the reference and alternate alleles are both single nucleotides
	if len(variant.ref) == 1:
		altlengths = [len(alt) for alt in variant.alts]
		if all([altlen == 1 for altlen in altlengths]):
			return True
		else:
			return False
	else:
		return False

# function to check if the variant is biallelic
def is_biallelic(variant):
	if len(variant.alts) == 1:
		return True
	else:
		return False

# check if the variant is an indel
def is_indel(variant):
	if len(variant.ref) > 1 or any([len(alt) > 1 for alt in variant.alts]):
		return True

# check if a sample is missing
def is_missing(sample):
	if None in sample["GT"]:
		return True
	else:
		return False
# check if sample is heterozygous
def is_het(sample):
	if sample["GT"][0] != sample["GT"][1]:
		return True
	else:
		return False
# check if sample is homozygous reference
def is_homref(sample):
	if sample["GT"] == (0,0):
		return True
	else:
		return False
# check if sample is homozygous alternate
def is_homalt(sample):
	if sample["GT"] != (0,0) and None not in sample["GT"]:
		return True
	else:
		return False

# lists for storing the results
globalstats = {
	'snps': 0,
	'indels': 0,
	'other': 0,
	'biallelic': 0,
	'multiallelic': 0,
	'invariant': 0,
}
samplestats = {}

processed = 0

# print that we're starting
print("Starting to process variants")

# get timestamp
start = dt.datetime.now()

# loop through the variants
for rec in vcf.fetch():
	# check if there is variation
	if rec.alts is None or len(rec.alts) == 0:
		globalstats['invariant'] += 1
	else:
		# check what kind of variant it is
		if is_snp(rec):
			vartype = "snp"
			globalstats['snps'] += 1
		elif is_indel(rec):
			vartype = "indel"
			globalstats['indels'] += 1
		else:
			vartype = "other"
			globalstats['other'] += 1
		# check if it's biallelic
		if is_biallelic(rec):
			biallelic = "yes"
			globalstats['biallelic'] += 1
		else:
			biallelic = "no"
			globalstats['multiallelic'] += 1
	# if the samples list is empty, populate it with the sample names
	if len(samplestats.keys()) == 0:
		for sample in list(vcf.header.samples):
			samplestats[sample] = {'called': 0, 'het': 0, 'hom': 0, 'hom_alt': 0, 'hom_ref': 0, 'missing': 0}
	# loop through the samples
	for sample in rec.samples:
		# check if the sample is missing
		if is_missing(rec.samples[sample]):
			samplestats[sample]['missing'] += 1
			continue
		# check if the sample is heterozygous
		if is_het(rec.samples[sample]):
			samplestats[sample]['het'] += 1
			samplestats[sample]['called'] += 1
			continue
		# check if the sample is homozygous reference
		if is_homref(rec.samples[sample]):
			samplestats[sample]['hom_ref'] += 1
			samplestats[sample]['hom'] += 1
			samplestats[sample]['called'] += 1
			continue
		# check if the sample is homozygous alternate
		if is_homalt(rec.samples[sample]):
			samplestats[sample]['hom_alt'] += 1
			samplestats[sample]['hom'] += 1
			samplestats[sample]['called'] += 1
			continue
	# report back every 100000 variants
	processed += 1
	if processed % 100000 == 0:
		seconds_elapsed = (dt.datetime.now() - start).total_seconds()
		print("Processed " + str(processed) + " variants in " + str(seconds_elapsed) + "seconds.")

# write the global stats to a file
with open(prefix + ".globalstats.txt", "w") as outfile:
	outfile.write("type\tcount\n")
	for stat in globalstats:
		outfile.write(stat + "\t" + str(globalstats[stat]) + "\n")

# write the sample stats to a file
with open(prefix + ".samplestats.txt", "w") as outfile:
	outfile.write("sample\tcalled\thet\thom\thom_alt\thom_ref\tmissing\n")
	for sample in samplestats:
		outfile.write(sample + "\t" + str(samplestats[sample]['called']) + "\t" + str(samplestats[sample]['het']) + "\t" + str(samplestats[sample]['hom']) + "\t" + str(samplestats[sample]['hom_alt']) + "\t" + str(samplestats[sample]['hom_ref']) + "\t" + str(samplestats[sample]['missing']) + "\n")

	