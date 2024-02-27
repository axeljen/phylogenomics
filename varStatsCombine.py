# small python script to combine output from varStats.py based on a list of input files
import sys
import os

input_prefixes = sys.argv[1]
output_prefix = sys.argv[2]

# make a list with the samplestats and globalstats files, respectively
samplestats = []
globalstats = []

if os.path.isfile(input_prefixes):
	for line in open(input_prefixes, 'r').readlines():
		samplestats.append(line.strip() + '.samplestats.txt')
		globalstats.append(line.strip() + '.globalstats.txt')
else:
	for prefix in input_prefixes.split(','):
		samplestats.append(prefix + '.samplestats.txt')
		globalstats.append(prefix + '.globalstats.txt')

# process the first of each list to get the header
samplestats_file = open(samplestats[0], 'r')
stats = {}
samplestats_file.readline().strip().split('\t')
for line in samplestats_file.readlines():
	sample = line.strip().split('\t')[0]
	stats[sample] = {
		'called': int(line.strip().split('\t')[1]),
		'het': int(line.strip().split('\t')[2]),
		'hom': int(line.strip().split('\t')[3]),
		'hom_alt': int(line.strip().split('\t')[4]),
		'hom_ref': int(line.strip().split('\t')[5]),
		'missing': int(line.strip().split('\t')[6]),
		}
# loop through the remaining samplestats
for samplestat in samplestats[1:]:
	samplestats_file = open(samplestat, 'r')
	samplestats_file.readline()
	for line in samplestats_file.readlines():
		sample = line.strip().split('\t')[0]
		# print error and exit if sample is not in the header
		if sample not in stats:
			print('ERROR: sample ' + sample + ' not in header of ' + samplestat)
			sys.exit()
		stats[sample]['called'] += int(line.strip().split('\t')[1])
		stats[sample]['het'] += int(line.strip().split('\t')[2])
		stats[sample]['hom'] += int(line.strip().split('\t')[3])
		stats[sample]['hom_alt'] += int(line.strip().split('\t')[4])
		stats[sample]['hom_ref'] += int(line.strip().split('\t')[5])
		stats[sample]['missing'] += int(line.strip().split('\t')[6])
	samplestats_file.close()

# enxt, process the globalstats
globstats = {
	'snps': 0,
	'indels': 0,
	'other': 0,
	'biallelic': 0,
	'multiallelic': 0,
	'invariant': 0,
}

for glob in globalstats:
	with open(glob, 'r') as g:
		for line in g.readlines():
			if line.startswith('type'):
				continue
			elif line.startswith('snps'):
				globstats['snps'] += int(line.strip().split('\t')[1])
				continue
			elif line.startswith('indels'):
				globstats['indels'] += int(line.strip().split('\t')[1])
				continue
			elif line.startswith('other'):
				globstats['other'] += int(line.strip().split('\t')[1])
				continue
			elif line.startswith('biallelic'):
				globstats['biallelic'] += int(line.strip().split('\t')[1])
				continue
			elif line.startswith('multiallelic'):
				globstats['multiallelic'] += int(line.strip().split('\t')[1])
				continue
			elif line.startswith('invariant'):
				globstats['invariant'] += int(line.strip().split('\t')[1])
				continue
			else:
				continue

# write the output files
with open(output_prefix + '.samplestats.txt', 'w') as outfile:
	outfile.write('sample\tcalled\thet\thom\thom_alt\thom_ref\tmissing\n')
	for sample in stats:
		outfile.write(sample + '\t' + str(stats[sample]['called']) + '\t' + str(stats[sample]['het']) + '\t' + str(stats[sample]['hom']) + '\t' + str(stats[sample]['hom_alt']) + '\t' + str(stats[sample]['hom_ref']) + '\t' + str(stats[sample]['missing']) + '\n')
with open(output_prefix + '.globalstats.txt', 'w') as outfile:
	outfile.write('type\tcount\n')
	for stat in globstats:
		outfile.write(stat + '\t' + str(globstats[stat]) + '\n')