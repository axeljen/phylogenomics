import pysam
import functions as fn
import argparse

args = argparse.ArgumentParser(description='Calculate heterozygosity in sliding windows along the genome')

args.add_argument('-i', '--input', help='Input file in VCF format', required=True)
args.add_argument('-o', '--output', help='Output file.', required=True)
args.add_argument('--ref-index', help="Reference genome index file, for generating windows", required=True)
args.add_argument('-w', '--window', help='Window size', required=True)
args.add_argument('-s', '--step', help='Step size, defaults to non-overlapping windows', required=False)
args.add_argument('--samples', help='Samples to include in the analysis. Defaults to all samples in vcf file.', required=False, default=None)
args = args.parse_args()

if args.step == None:
	args.step = args.window

def parseSamples(vcf, samples = None):
	if not samples == None:
		samples = []
		for l in samples.readlines():
			samples.append(l.strip())
	else:
		samples = list(vcf.header.samples)
	return samples

def parse_genotype(genotype):
	if None in genotype:
		return "missing"
	else:
		if genotype[0] == genotype[1]:
			return "hom"
		else:
			return "het"

def prep_window_dict(chrom,start,end,samples):
	window_dict = {}
	window_dict['chrom'] = chrom
	window_dict['start'] = start
	window_dict['end'] = end
	window_dict['n_sites'] = 0
	window_dict['samples'] = {}
	for sample in samples:
		window_dict['samples'][sample] = {'het':0,'hom':0,'missing':0}
	return window_dict
	

# open vcf file
vcf = pysam.VariantFile(args.input)

# get chromosome lengths
chromLengths = fn.getChromLengths(vcf, only_contigs_with_records=True, reference_index = args.ref_index)

# make windows
windows = fn.generateWindows(chromLengths, int(args.window), int(args.step))

# get samples
samples = parseSamples(vcf, args.samples)

# list for storing results
results = []

# iterate over windows
done_windows = 0

print("Starting to process {} windows.".format(len(windows)))
for w in windows:
	window_dict = prep_window_dict(w['chrom'],w['start'],w['end'],samples)
	for record in vcf.fetch(w['chrom'],w['start'],w['end']):
		for sample in samples:
			state = parse_genotype(record.samples[sample]['GT'])
			window_dict['samples'][sample][state] += 1
		window_dict['n_sites'] += 1
	results.append(window_dict)
	done_windows += 1
	if done_windows % 100 == 0:
		print("Processed {}/{} windows.".format(done_windows,len(windows)))

# write results to file
with open(args.output,'w') as out:
	out.write('chrom\tstart\tend\tn_sites\tsample\tmissing\thet\thom\n')
	for r in results:
		chrom = r['chrom']
		start = r['start']
		end = r['end']
		n_sites = r['n_sites']
		for sample in samples:
			missing = r['samples'][sample]['missing']
			het = r['samples'][sample]['het']
			hom = r['samples'][sample]['hom']
			out.write(f'{chrom}\t{start}\t{end}\t{n_sites}\t{sample}\t{missing}\t{het}\t{hom}\n')

