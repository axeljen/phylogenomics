# Axel Jensen 2023
# this script takes a vcf file as input and prepares site pattern counts for a dfoil analyses, as specified on the dfoil manual.

import argparse
import sys
import os
import functions as fn
import pysam

class dfoilCounts:
	def __init__(self,vcf,chrom,start,end,P1,P2,P3,P4,Out):
		self.vcf = vcf
		self.chrom = chrom
		self.start = start
		self.end = end
		self.pops = {'P1': P1,'P2':P2,'P3':P3,'P4':P4}
		self.Out = Out
		self.counts = {'AAAAA':0,'AAABA':0,'AABAA':0,'AABBA':0,'ABAAA':0,'ABABA':0,'ABBAA':0,'ABBBA':0,
		'BAAAA':0,'BAABA':0,'BABAA':0,'BABBA':0,'BBAAA':0,'BBABA':0,'BBBAA':0,'BBBBA':0}
		self.goodsites = 0
		self.badsites = 0
	def countpatterns(self):
		# print some info
		if self.chrom:
			sys.stderr.write('Counting site patterns on locus {}:{}-{}.\n'.format(self.chrom,self.start,self.end))
		else:
			sys.stderr.write('Counting site patterns across all sites in vcf file.\n')
		# parse through the vcf
		for rec in vcf.fetch(self.chrom,self.start,self.end):
			# we only want biallelic sites
			if not len(rec.alleles) == 2:
				continue
			# check that the outgroup is fixed, and in this case we use this allele as A, otherwise continue
			allele_a = []
			for sample in self.Out:
				if not None in rec.samples[sample]['GT']:
					allele_a.append(rec.samples[sample]['GT'][0])
					allele_a.append(rec.samples[sample]['GT'][1])
			if not len(set(allele_a)) == 1:
				continue
			else:
				allele_a = allele_a[0]
				if allele_a == 0:
					allele_b = 1
				elif allele_a == 1:
					allele_b = 0
				else:
					continue
				# get the frequencies of derived allele in each pop
				freqs = {'P1': {'A':0,'B':0},'P2': {'A':0,'B':0}, 'P3':{'A':0,'B':0}, 'P4':{'A':0,'B':0}, 'Out': {'A':1,'B':0}}
				counts = {p: [0,0] for p in self.pops} # for each pop a tuple with derived allele count (allele_b) and genotyped alleles
				for pop,samples in self.pops.items():
					#print(pop)
					#print(samples)
					for sample in samples:
						if not None in rec.samples[sample]['GT']:
							counts[pop][1] += 2
							counts[pop][0] += rec.samples[sample]['GT'].count(allele_b)
					#print(counts[pop])
					# go from counts to frequencies
					try:
						freqs[pop]['B'] = counts[pop][0] / counts[pop][1]
						freqs[pop]['A'] = 1 - (counts[pop][0] / counts[pop][1])
					except:
						# if zerodivision error, set to none and skip in next step
						freqs[pop] = None
				#print(counts)
				#print(freqs)
				if not None in freqs.values():
					self.goodsites += 1
					# if none of the pops are missing, add the freqs to the counts
					for sp in self.counts.keys():
						freq = freqs['P1'][sp[0]] * freqs['P2'][sp[1]] * freqs['P3'][sp[2]] * freqs['P4'][sp[3]] * freqs['Out'][sp[4]]
						self.counts[sp] += freq
						if freq > 0:
							#print('{}: {}, {}\n'.format(rec.pos,rec.ref,rec.alts))
							#print(sp)
							# for testing, fetch all alleles in the correct order
							alleles_t = {}
							for pop,samples in self.pops.items():
								alleles_t[pop] = [rec.samples[sample]['GT'] for sample in samples]
							#print("alleles: {}".format(alleles_t))
							#print("{},{},{},{},{}".format(freqs['P1'][sp[0]],freqs['P2'][sp[1]],freqs['P3'][sp[2]],freqs['P4'][sp[3]],freqs['Out'][sp[4]]))
							#print("Sumfreq: {}".format(freqs['P1'][sp[0]] * freqs['P2'][sp[1]] * freqs['P3'][sp[2]] * freqs['P4'][sp[3]] * freqs['Out'][sp[4]]))
					#print("\n\n")
				else:
					self.badsites += 1

#set up and parse input arguments
parser = argparse.ArgumentParser(description="Count site patterns to use as input for dfoil analyses. This script use polymorphic data from vcf file, and round the final counts to integers as dfoil doesn't allow floats.")

parser.add_argument('-v','--vcf', type=str, help='Input vcf')

parser.add_argument('-o','--out', type=str, help='Output file with site pattern counts.')


parser.add_argument('-w','--window-size', type=int, help='Size of window for counting site pattern. If no window is given, counts will be given across entire vcf.')
parser.add_argument('-s','--step-size', type=int, help='Offset between windows (start - next start')

parser.add_argument('-p', '--pop-file', type=str, help='File with population assignments.', required=True)

parser.add_argument('-P1', type=str,help='P1 population')
parser.add_argument('-P2', type=str,help='P2 population')
parser.add_argument('-P3', type=str,help='P3 population')
parser.add_argument('-P4', type=str,help='P4 population')
parser.add_argument('-O', type=str,help='Outgroup population')

args = parser.parse_args()

# parse the popfile
popsparsed = {}
with open(args.pop_file) as of:
	for line in of.readlines():
		try:
			sample,pop = line.strip().split()
			if pop in popsparsed.keys():
				popsparsed[pop].append(sample)
			else:
				popsparsed[pop] = [sample]
		except:
			continue

# make list of P1-Out based on popfile
P1 = popsparsed[args.P1]
P2 = popsparsed[args.P2]
P3 = popsparsed[args.P3]
P4 = popsparsed[args.P4]
Out = popsparsed[args.O]

# open the vcf
vcf = pysam.VariantFile(args.vcf)
if args.window_size:
	# generate windows based on chroms and records in input vcf
	chroms = fn.getChromLengths(vcf)
	windows = fn.generateWindows(chroms,window=args.window_size,step_size=args.step_size)
	sys.stderr.write('Counting sitepatterns in {} windows. Parameters: window size = {}, step size = {}\n'.format(len(windows),args.window_size,args.step_size))
	# start a list dictionary for storing results
	results = []

	# now start parsing through the windows
	for w in windows:
		# initiate dfoil counts object
		dc = dfoilCounts(vcf,w['chrom'],w['start'],w['end'],P1,P2,P3,P4,Out)
		# count the patterns
		dc.countpatterns()
		# append results
		results.append({'chrom':dc.chrom,'start':dc.start,'end':dc.end,'mid':dc.end - (dc.end - dc.start)/2, 'counts':dc.counts})
else:
	# do the counts across the full file
	dc = dfoilCounts(vcf,None,None,None,P1,P2,P3,P4,Out)
	dc.countpatterns()
	results = [{'chrom':'na','start':'na','end':'na','mid':1, 'counts':dc.counts}]

print(results[0]['counts'])

# write the output file
with open(args.out, 'w') as of:
	of.write('#chrom\tmidpos\t{}\n'.format('\t'.join([h for h in results[0]['counts'].keys()])))
	for r in results:
		of.write('{chrom}\t{midpos}\t{counts}\n'.format(chrom=r['chrom'],midpos=str(round(r['mid'])),counts = '\t'.join([str(round(c)) for c in r['counts'].values()])))

sys.stderr.write("Done.\n")

