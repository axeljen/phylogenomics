import sys
from pysam import VariantFile
import Bio
from Bio import pairwise2
import random
import re
import subprocess
import tempfile
import os
# function to translate to complement allele
def complement(allele):
	# table of complementary bases
	comps = {
		'A':'T','T':'A','C':'G','G':'C',
		'a':'t','t':'a','c':'g','g':'c',
		'Y':'R','R':'Y','S':'S','W':'W','K':'M','M':'K','B':'V','D':'H','H':'D','V':'B','N':'N',
		'y':'r','r':'y','s':'s','w':'w','k':'k','m':'m','b':'v','d':'h','h':'d','v':'b','n':'n',
	}
	return comps[allele]

# Function to go from a diploid genotype (two alleles) to a IUPAC code
def IUPAC(alleles):
	IUPAC = {
		'M': ['A','C'],
		'R': ['A','G'],
		'W': ['A','T'],
		'S': ['C','G'],
		'Y': ['C','T'],
		'K': ['G','T']
	}
	hit = False
	for c in IUPAC:
		if set(IUPAC[c]) == set(alleles):
			code = c
			hit = True
	if hit == False:
		code = "N"
	return code

def reverseIUPAC(iupac):
	IUPAC = {
		'M': ['A','C'],
		'R': ['A','G'],
		'W': ['A','T'],
		'S': ['C','G'],
		'Y': ['C','T'],
		'K': ['G','T']
	}
	if iupac in IUPAC.keys():
		return IUPAC[iupac][random.randint(0,1)]
	else:
		sys.stderr.write("Invalid input base for this function: {}. Returning 'N'".format(iupac))
		return "N"

# function to convert a codon to an amino acid
def codon2aa(codon):
	# translation table
	translation = {
		'AAA':'K','AAC':'N','AAG':'K','AAT':'N', 'AAR':'K', 'AAY':'N',
		'ACA':'T','ACC':'T','ACG':'T','ACT':'T', 'ACN':'T', 'ACR':'T', 'ACY':'T', 'ACW':'T', 'ACM':'T', 'ACS':'T', 'ACK':'T', 'ACH':'T', 'ACB':'T', 'ACD':'T', 'ACV':'T',
		'AGA':'R','AGC':'S','AGG':'R','AGT':'S', 'AGR':'R', 'AGY':'S',
		'ATA':'I','ATC':'I','ATG':'M','ATT':'I', 'ATY':'I', 'ATW':'I', 'ATM':'I', 'ATH':'I', 'ATB':'I', 'ATD':'I',
		'CAA':'Q','CAC':'H','CAG':'Q','CAT':'H', 'CAR':'Q', 'CAY':'H',
		'CCA':'P','CCC':'P','CCG':'P','CCT':'P', 'CCN':'P', 'CCR':'P', 'CCY':'P', 'CCW':'P', 'CCM':'P', 'CCS':'P', 'CCK':'P', 'CCH':'P', 'CCB':'P', 'CCD':'P', 'CCV':'P',
		'CGA':'R','CGC':'R','CGG':'R','CGT':'R', 'CGN':'R', 'CGY':'R', 'CGW':'R', 'CGM':'R', 'CGS':'R', 'CGK':'R', 'CGH':'R', 'CGB':'R', 'CGD':'R', 'CGV':'R', 'MGR':'R',
		'CTA':'L','CTC':'L','CTG':'L','CTT':'L', 'CTN':'L', 'CTR':'L', 'CTY':'L', 'CTW':'L', 'CTM':'L', 'CTS':'L', 'CTK':'L', 'CTH':'L', 'CTB':'L', 'CTD':'L', 'CTV':'L',
		'GAA':'E','GAC':'D','GAG':'E','GAT':'D', 'GAR':'D', 'GAY':'D',
		'GCA':'A','GCC':'A','GCG':'A','GCT':'A', 'GCR':'A', 'GCS':'A', 'GCK':'A', 'GCW':'A', 'GCM':'A', 'GCY':'A', 'GCH':'A', 'GCB':'A', 'GCN':'A', 'GCD':'A', 'GCV':'A',
		'GGA':'G','GGC':'G','GGG':'G','GGT':'G', 'GGN':'G', 'GGR':'G', 'GGY':'G', 'GGW':'G', 'GGM':'G', 'GGS':'G', 'GGK':'G', 'GGH':'G', 'GGB':'G', 'GGD':'G', 'GGV':'G',
		'GTA':'V','GTC':'V','GTG':'V','GTT':'V', 'GTR':'V', 'GTY':'V',
		'TAA':'*','TAC':'Y','TAG':'*','TAT':'Y', 'TAR':'*', 'TAY':'Y',
		'TCA':'S','TCC':'S','TCG':'S','TCT':'S', 'TCN':'S',
		'TGA':'*','TGC':'C','TGG':'W','TGT':'C', 'TGY':'C',
		'TTA':'L','TTC':'F','TTG':'L','TTT':'F', 'TTR':'L', 'TTY':'F',
	}
	# capitialize codon
	codon = codon.upper()
	# if codon is not in the translation table, return ?
	if not codon in translation.keys():
		return "?"
	else:
		return translation[codon]

# function to convert a whole sequence to amino acids
def translate(seq, break_on_stop = False):
	# make a list of codons by splitting the sequence in triplets
	codons = [seq[i:i+3] for i in range(0, len(seq), 3)]
	aaseq = ''
	# loop over the codons and translate them
	for i in range(len(codons)):
		aa = codon2aa(codons[i])
		if aa == "*":
			# print a warning if this is not the last codon
			if i < len(codons) - 1:
				print("Warning: internal stop codon encountered at position {}\n".format(len(aaseq) + 1))
				if break_on_stop:
					sys.exit()
		aaseq += aa
	return aaseq

# class for storing single sequence
class Sequence:
	def __init__(self, sample, sequence = "", meta = None):
		self.sample = sample
		self.sequence = sequence
		self.meta = meta
	def addSequence(self,sequence):
		self.sequence = self.sequence + sequence
	def removeBaseByPosition(self,position):
		# remove a base from the sequence by its zero based index
		self.sequence = self.sequence[0:position] + self.sequence[position + 1:]
	def reverseComplement(self):
		revcomp = ""
		for i in range(1,len(self.sequence) + 1):
			revcomp += complement(self.sequence[-i])
		self.sequence = revcomp
	def subsetSequence(self,start,end,index=0):
		start = start - index
		self.sequence = self.sequence[start:end]
	def splitCodons(self,offset=0,not_divisible_by_three='warning'):
		seqs = []
		# check if it's divisible by three or not and act accordingly
		if not len(self.sequence) % 3 == 0:
			if not_divisible_by_three == 'warning':
				sys.stderr.write('Warning: sequence {} of length {} is not divisible by three. Might result in erroneous codon chopping.'.format(self.sample,len(self.sequence)))
			else:
				sys.stderr.write('Error: sequence {} of length {} is not divisible by three. Remove --strict-triplets to proceed anyway.'.format(self.sample,len(self.sequence)))
				sys.exit()
		seqs.append(Sequence(sample = self.sample, meta = {'info':"first_codonpos_concat"}, sequence = self.sequence[offset::3]))
		seqs.append(Sequence(sample = self.sample, meta = {'info':"second_codonpos_concat"}, sequence = self.sequence[offset + 1::3]))
		seqs.append(Sequence(sample = self.sample, meta = {'info':"third_codonpos_concat"}, sequence = self.sequence[offset + 2::3]))
		return seqs
	#function to split a sequence into kmers and return unique kmers and their counts
	def countKmers(self, k):
		kmers = {}
		seqlen = len(self.sequence)
		startpos = 0
		endpos = startpos + k
		while endpos <= seqlen:
			kmer = self.sequence[startpos:endpos]
			if not kmer in kmers.keys():
				kmers[kmer] = 1
			else:
				kmers[kmer] += 1
			startpos += 1
			endpos = startpos + k
		return kmers
	# function to blast a query sequence to this one
	def blastToMe(self,queryseq, blast='blastn'):
		# write this sequence to temporary file
		subject = tempfile.NamedTemporaryFile()
		subject.write(str.encode('>{name}\n{seq}\n'.format(name = self.sample, seq=self.sequence)))
		subject.seek(0)
		# and the query sequence as well
		query = tempfile.NamedTemporaryFile()
		# first, if the query sequence is given as an input file, read it in first
		if os.path.isfile(queryseq):
			msa = readSequenceFile(queryseq)
			seqobject = msa.sequences[list(msa.sequences.keys())[0]]
			name,queryseq = seqobject.sample,seqobject.sequence
		else:
			# otherwise just call it query
			name = 'query'
		# write to tempfile
		query.write(str.encode('>{name}\n{queryseq}\n'.format(name = name, queryseq=queryseq)))
		# reset marker
		query.seek(0)
		# run blast
		blast = subprocess.run([blast, "-subject", subject.name, '-query', query.name, '-outfmt','7'], 
		capture_output = True)
		# grab blast raw output
		results = blast.stdout.decode('ASCII').split("\n")
		# fetch the result fieldsname
		fieldsline = [f for f in results if f.startswith("# Fields:")][0]
		fields = [f.strip() for f in fieldsline.split(":")[1].split(",")]
		resultslines = [i.strip().split() for i in results if not i.startswith("#") and len(i.split("\t")) > 0]
		resDict = {}
		for i,res in enumerate(resultslines):
			if len(res) == 0:
				continue
			resDict["hit_" + str(i + 1)] = {fields[f]: res[f] for f in range(0, len(fields))}
			# if any hits, also append an msa object to the resdict for fun
			#msa_1 = MSA()
			#msa_1.addSample(name = resDict["hit_" + str(i + 1)]['query acc.ver'], seq = queryseq)
			#msa_1.addSample(name = resDict["hit_" + str(i + 1)]['subject acc.ver'], seq = self.sequence[int(resDict["hit_" + str(i + 1)]['s. start']) - 1:int(resDict["hit_" + str(i + 1)]['s. end'])])
			#resDict["hit_" + str(i + 1)]['msa'] = msa_1
		return resDict
	# function to blast this sequence to a query
	def blastMeTo(self,subjectseq, blast='blastn'):
		# write this sequence to temporary file
		query = tempfile.NamedTemporaryFile()
		query.write(str.encode('>{name}\n{seq}\n'.format(name = self.sample, seq=self.sequence)))
		query.seek(0)
		# and the query sequence as well
		subject = tempfile.NamedTemporaryFile()
		# first, if the query sequence is given as an input file, read it in first
		if os.path.isfile(subjectseq):
			msa = readSequenceFile(subjectseq)
			seqobject = msa.sequences[list(msa.sequences.keys())[0]]
			name,subjectseq = seqobject.name,seqobject.sequence
		else:
			# otherwise just call it query
			name = 'subject'
		# write to tempfile
		subject.write(str.encode('>{name}\n{subjectseq}\n'.format(name = name, subjectseq=subjectseq)))
		# reset marker
		subject.seek(0)
		# run blast
		blast = subprocess.run([blast, "-subject", subject.name, '-query', query.name, '-outfmt','7'], 
		capture_output = True)
		# grab blast raw output
		results = blast.stdout.decode('ASCII').split("\n")
		# fetch the result fieldsname
		fieldsline = [f for f in results if f.startswith("# Fields:")][0]
		fields = [f.strip() for f in fieldsline.split(":")[1].split(",")]
		resultslines = [i.strip().split() for i in results if not i.startswith("#") and len(i.split()) > 0]
		resDict = {}
		for i,res in enumerate(resultslines):
			resDict["hit_" + str(i + 1)] = {fields[f]: res[f] for f in range(0, len(fields))}
			# if any hits, also append an msa object to the resdict for fun
			msa_1 = MSA()
			msa_1.addSample(name = resDict["hit_" + str(i + 1)]['query acc.ver'], seq = queryseq)
			msa_1.addSample(name = resDict["hit_" + str(i + 1)]['subject acc.ver'], seq = self.sequence[int(resDict["hit_" + str(i + 1)]['s. start']) - 1:int(resDict["hit_" + str(i + 1)]['s. end'])])
			resDict["hit_" + str(i + 1)]['msa'] = msa_1
		return resDict
	def upperCase(self):
		# will change all nucleotides to be of upper case in sequence
		self.sequence = self.sequence.upper()
	def distance(self, sequence):
		# calculate the distance between to aligned sequences
		# first make sure they're of the same length
		if not len(self.sequence) == len(sequence.sequence):
			sys.stderr.write("Error: sequences are not of the same length. Exiting.")
			sys.exit(1)
		# then calculate the distance
		dist = 0
		comparisons = 0
		for i in range(0,len(self.sequence)):
			(n1,n2) = self.sequence[i], sequence.sequence[i]
			if n1 in ['-','?','N','n'] or n2 in ['-','?','N','n']:
				continue
			else:
				comparisons += 1
				if n1 != n2:
					dist += 1
		return dist, comparisons
	def missingness(self):
		# calculate the proportion of missing data in this sequence
		missing = 0
		for i in self.sequence:
			if i in ['-','?','N','n']:
				missing += 1
		return missing / len(self.sequence)
	def containsInternalStopCodon(self):
		# check if the sequence contains an internal stop codon
		# first translate the sequence
		aa = translate(self.sequence)
		# remove the last one as this could be the real stop 
		aa = aa[:-1]
		if '*' in aa:
			return True
		else:
			return False
# class for storing multi sequence (alignment only for now, but should work fine with non-aligned as well?)
class MSA:
	def __init__(self, samples = [], file = None, chrom=None,start=None,end=None):
		self.samples = samples
		self.sequences = {}
		self.called_positions = []
		self.chrom = chrom
		self.start = start
		self.end = end
		self.length = 0
		self.missingness = []
		self.partitions = []
		self.name = None
		self.reference = None
	# check lengths of sequences
	def checkAlnLength(self):
		# function to quickly check that all sequences are of the same length, just as a sanity check
		if len(self.sequences) == 0:
			# if all sequences have been filtered out
			sys.stderr.write("Error: no sequences in alignment. Exiting.\n")
			sys.exit(1)
		lengths = [len(i.sequence) for i in self.sequences.values()]
		if len(set(lengths)) == 1:
			self.length = lengths[0]
		else:
			# if they're not we should exit and print that something's up
			sys.stderr.write("Sequences differ in length, which they shouldn't. Exiting.\n")
			#sys.exit(1)
		# use this function to update sequence missingness too
		self.missingness = [0 for i in range(0,self.length)]
		for i in range(0,self.length):
			for seq in self.sequences.values():
				if seq.sequence[i] in ['-','N','n','?']:
					self.missingness[i] += 1
			self.missingness[i] = self.missingness[i] / len(self.sequences)
	# add sample
	def addSample(self,name,seq,meta=None):
		seq = Sequence(name,seq,meta = meta)
		self.samples.append(name)
		self.sequences[name] = seq
		self.checkAlnLength()
	@classmethod
	def fromFasta(cls,fastafile):
		# this method fetches a sequence alignment object from a fasta file.
		# initiate an empty MSA object
		msa = MSA()
		msa.samples = []
		# keep track of parsed sequences
		seqcount = 0
		# read file, line by line, and parse to sequence objects
		# start with an empty sample and sequence string
		sample = ""
		sequence = ""
		with open(fastafile) as f:
			for line in f.readlines():
				if line.startswith(">"):
					# start a new sequence/sample
					if sequence != "":
						# add previous sequence to msa
						msa.addSample(sample,sequence)
						# start a new sequence
						sample = line.strip().lstrip(">")
						sequence = ""
					else:
						# if this is our first sequence
						sample = line.strip().lstrip(">")
						sequence = ""
				else:
					sequence += line.strip()
		# add the last sequence
		msa.addSample(sample,sequence)
		# and make sure that the length is as should be
		msa.checkAlnLength()
		return msa
	@classmethod
	def fromPhylip(cls, phylipfile, interleaved=False):
		# initiate msa object
		msa = MSA()
		msa.samples = []
		f = open(phylipfile)
		# first line will contain number and length of sequences, can be split by arbitrary spaces, following solution should get them no matter what
		firstline = [i for i in f.readline().split() if not i == ""]
		nseq,length = int(firstline[0]),int(firstline[1])
		# then continue to read lines and fetch the samples and seqs, much like with the fasta function
		sample = ""
		sequence = ""
		current_sample = 0
		for i,line in enumerate(f.readlines()):
			# if interleaved format, make sure we don't go over the number of samples
			if i <= nseq +1:
				sample = line.split()[0].strip()
				sequence = ''.join(line.split()[1].split()).strip()
				msa.addSample(sample,sequence)
			else:
				# otherwise, go back to the beginning and add the remaining sequence to the sample
				msa.sequences[sample].addSequence(''.join(line.split().strip()))
				if current_sample == nseq - 1:
					# if so, start over from the first sample again
					current_sample = 0
				else:
					current_sample += 1
		# check that length and samples are what we expect
		msa.checkAlnLength()
		if not msa.length == length:
			sys.stderr.write("Error: alignment length doesn't match header.")
		if not len(msa.sequences) == nseq:
			sys.stderr.write("Error: number of sequences doesn't match header.")
		return msa
	# fetch a vcf region and turn into an msa object
	@classmethod
	def fromVcf(cls,vcf,chrom,start,end,samples=None,haploidize=True,heterozygotes="IUPAC", padding="N", skip_indels = True, index = 1, verbosity='err', minsampleDP=None,maxsampleDP=None, ploidy=2):
		# if sample is none take all samples
		#print(samples)
		if not samples:
			samples = list(vcf.header.samples)
		if start:
			start = start - index
		#except:
		#	sys.stderr.write("Coordinates (chrom:start-end) are required for now. Please specify.")
		#	sys.exit()
		msa = MSA(samples = samples,chrom = chrom, start = start, end = end)
		# create a sequence object for all the samples
		msa.sequences = {}
		for sample in samples:
			if haploidize == True or ploidy == 1:
				msa.sequences[sample] = Sequence(sample,meta = {'vcfsample':sample, 'haplotype':0,})
			else:
				msa.sequences[sample + "__1"] = Sequence(sample + "__1", meta = {'vcfsample':sample,'haplotype':0,})
				msa.sequences[sample + "__2"] = Sequence(sample + "__2", meta = {'vcfsample':sample,'haplotype':1,})
		# go through the records and parse to letters
		#print(start)
		for rec in vcf.fetch(chrom,start,end):
			# sanity check to make sure that ref and alt alleles are regular, single bases
			if not rec.ref in ['A','a','T','t','C','c','G','g']:
				if verbosity == 'warn':
					sys.stderr.write("Warning: non-SNP ref allele ({}) encountered at pos {}, skipping this site.\n".format(rec.ref, rec.pos))
				continue
			skip = False
			if rec.alts:
				for a in rec.alts:
					if not a in ['A','a','T','t','C','c','G','g']:
						if verbosity == 'warn':
							sys.stderr.write("Warning: non-SNP alt allele ({}) encountered at pos {}, skipping this site.\n".format(a, rec.pos))
						skip = True
						continue
			if skip == True:
				continue
			msa.called_positions.append(rec.pos) # this one is in case we want to project the alignment onto a contigous refgenome
			missing_calls = 0 # keep track of missingness in case we want to filter on this
			# get all the alleles
			alleles = list(rec.alleles)
			#print(alleles)
			# if the alleles are of different lengths, add padding so that they come out the same ### This is not in use for now. Only allowing SNPS and invars
			longest = 1
			if not len(set(len(i) for i in alleles)) == 1:
				if skip_indels:
					continue
				else:
					for a in alleles:
						if len(a) > longest:
							longest = len(a)
				alleles = [a + padding * (longest - len(a)) for a in alleles]
			# add missing allele to end of this list
			alleles.append(padding * longest)
			# then just go through the sequences and add the base to the sequences
			for seq in msa.sequences.values():
				# first get the samples genotype
				gt = rec.samples[seq.meta['vcfsample']]['GT']
				# check if there are depthfilters to consider
				if minsampleDP is not None:
					if rec.samples[seq.meta['vcfsample']]['DP'] is None or rec.samples[seq.meta['vcfsample']]['DP'] < minsampleDP:
						# set gt to none
						gt = (None,None)
				if maxsampleDP is not None:
					if rec.samples[seq.meta['vcfsample']]['DP'] is None or rec.samples[seq.meta['vcfsample']]['DP'] > maxsampleDP:
						# set gt to none
						gt = (None,None)
				# if heterozygote
				if not haploidize:
					if None in gt:
						msa.sequences[seq.sample].addSequence(alleles[-1])
						missing_calls += 1
					else:
						# in this case output the diploid genotypes to two different sequences
						msa.sequences[seq.sample].addSequence(alleles[gt[seq.meta['haplotype']]])
				else:
					if None in gt:
						missing_calls += 1
						# if genotype is missing, add the specified missing letter
						msa.sequences[seq.sample].addSequence(alleles[-1])
					if ploidy == 2:
						if gt[0] != gt[1]:
							if heterozygotes == "IUPAC":
								# grab the iupac code corresponding to the heterozygot call
								msa.sequences[seq.sample].addSequence(IUPAC((alleles[gt[0]],alleles[gt[1]])))
							elif heterozygotes == "random":
								# otherwise just take a random base
								msa.sequences[seq.sample].addSequence(alleles[gt[random.randint(0,1)]])
							elif heterozygotes == "first":
								# always use the first allele
								msa.sequences[seq.sample].addSequence(alleles[gt[0]])
							elif heterozygotes == "second":
								# always use the second allele
								msa.sequences[seq.sample].addSequence(alleles[gt[1]])
							else:
								sys.stderr.write("Please provide a valid arument on how to handle heterozygote genotypes. Exiting.")
								sys.exit(1)
						else:
							# if homozygous, just take the first genotype
							msa.sequences[seq.sample].addSequence(alleles[gt[0]])
					elif ploidy == 1:
						#print("Genotype for sequence " + seq.sample + ": " + alleles[gt[0]])
						# if haploid, just take the first genotype
						msa.sequences[seq.sample].addSequence(alleles[gt[0]])
					
			# add missingness info for subsequent filter opportunity
			msa.missingness.append(missing_calls / len(msa.samples))
		# check that lengths are ok before proceeding
		msa.checkAlnLength()
		return msa
	# filter msa object based on missingness
	def filterMSA(self, max_missingness):
		# this function takes an msa object and missingness threshold to filter out sites with high missingness
		called_pos_filtered = [] # keep track of the called positions after filtration
		missingness_filtered = [] # and missingness as well
		sequences_filtered = {seq: Sequence(seq) for seq in self.sequences.keys()} # keep filtered sequences in a separate dict for starters
		filtered = 0
		for i in range(self.length):
			try:
				m = self.missingness[i]
			except:
				sys.stderr.write("Failed at {chr}:{start}-{end}. Missingness length: {l}\nfailed at {i}.\nAln length: {y}.".format(chr=self.chrom,start=self.start,end=self.end,l=len(self.missingness),i=i,y=self.length))
			if m > max_missingness:
				filtered += 1
				continue
			else:
				for seq in sequences_filtered.keys():
					sequences_filtered[seq].addSequence(self.sequences[seq].sequence[i])
				called_pos_filtered.append(self.called_positions[i])
				missingness_filtered.append(self.missingness[i])
		# and replace the original values with the new
		self.called_positions = called_pos_filtered
		self.missingness = missingness_filtered
		self.sequences = {sample: sequence for sample,sequence in sequences_filtered.items()}
		# double check so lengths are still ok
		self.checkAlnLength()
		# print how much we filtered out
		sys.stderr.write("Filtered out {} sites with more than {} prop missing data\n".format(filtered,max_missingness))
		sys.stderr.write("Filtered alignment length: {}\n".format(self.length))
	def removePositions(self,positions):
		# this functions takes a list of positions to remove from the alignment, and removes those from all sequences
		for name,seq in self.sequences.items():
			newseq = ""
			for i in range(0,self.length):
				if not i in positions:
					newseq = newseq + seq.sequence[i]
			seq.sequence = newseq
		self.checkAlnLength()
	# function to insert uncalled ref bases to alignment, such that its lenght will be end - start coordinates
	def addRefGenome(self,refseq,start,end,refname='reference', mask_uncalled = True, strand = "+"):
		#print(refseq)
		# will have all the new sequences as a list of same length as reference to begin with, either with N's or refbases depending on mask_uncalled
		if mask_uncalled:
			seqs = {seq.sample: list("N" * len(refseq)) for seq in self.sequences.values()}
		else:
			seqs = {seq.sample: list(refseq) for seq in self.sequences.values()}
		# and then we simply loop through the called positions list and add the called base to each sequence
		if len(self.called_positions) > 0:
			# loop through each position in the reference genome, if strand is minus we need to reverse first
			positions = [i for i in range(start,end + 1)]
			if strand == '-':
				positions.reverse()
				# and then also reverse the called positions
				self.called_positions.reverse()
			for i,pos in enumerate(positions):
				if pos in self.called_positions:
					# find the position of this called position in the sequence
					seqpos = self.called_positions.index(pos)
					for sample,seq in self.sequences.items():
						seqs[sample][i] = self.sequences[sample].sequence[seqpos]
		else:
			for sample,seq in self.sequences.items():
				self.sequences[sample].sequence = "N" * len(refseq)
		# then, join the list to string and replace the previous sequences
		for sample,seq in self.sequences.items():
			self.sequences[sample].sequence = ''.join(seqs[sample])
		#print(refseq)
		# add the reference sequence to the alignment
		self.reference = Sequence(refname,refseq)
		# if minus strand, reverse complement it
		if strand == "-":
			self.reference.reverseComplement()
		# and update the msa length
		self.checkAlnLength()
	# reverse complement all sequences in the alignment
	def reverseComplement(self):
		for seq in self.sequences.values():
			seq.reverseComplement()
	# function to write the alignment to fasta format
	def writeFasta(self,filepath,wraplength = None, names = None):
		if names:
			# rename before writing if prompted
			self.renameSequences(names)
		with open(filepath, "w") as of:
			for seq in self.sequences.values():
				if not wraplength:
					of.write(">{sample}\n{sequence}\n".format(sample = seq.sample, sequence = seq.sequence))
				else:
					of.write(">{sample}\n".format(sample = seq.sample))
					chunks = [seq.sequence[i:i+wraplength] for i in range(0, len(seq.sequence), wraplength)]
					for chunk in chunks:
						of.write(chunk + "\n")
			if self.reference:
				# if there's a ref sequence, we write this too
				if not wraplength:
					of.write(">{}\n{}\n".format(self.reference.sample, self.reference.sequence))
				else:
					of.write(">{}\n".format(self.reference.sample))
					chunks = [self.reference.sequence[i:i+wraplength] for i in range(0, len(self.reference.sequence), wraplength)]
					for chunk in chunks:
						of.write(chunk + "\n")
		of.close()
	# write phylip output
	def writePhylip(self,filepath,strict=False, names = None):
		if names:
			# rename before writing if prompted
			self.renameSequences(names)
		nseq = len(self.sequences)
		if self.reference:
			# if there's a reference sequence too, we add this to nseq
			nseq += 1
		with open(filepath, "w") as of:
			of.write("\t{}\t{}\n".format(nseq, self.length))
			for seq in self.sequences.values():
				if strict:
					name = seq.sample[0:10]
				else:
					name = seq.sample
				of.write(name + " " * 3 + seq.sequence + "\n")
			if self.reference:
				# write reference too if there
				if strict:
					name = self.reference.sample[0:10]
				else:
					name = self.reference.sample
				of.write(name + " " * 3 + self.reference.sequence + "\n")
	def print(self, format = "fasta", positions=None):
		if format == "fasta":
			seqstring = ""
			for seq in self.sequences.values():
				if positions == None:
					start,end = 0,self.length
				else:
					start,end = positions[0],positions[1]
				#print(">{sample}\n{sequence}".format(sample = seq.sample, sequence = ''.join([i for i in seq.sequence[start:end]])))
				seqstring += ">{sample}\n{sequence}\n".format(sample = seq.sample, sequence = ''.join([i for i in seq.sequence[start:end]]))
			return seqstring
		elif format == "phylip":
			print("\t{}\t{}".format(len(self.sequences), self.length))
			for seq in self.sequences.values():
				print(seq.sample + " " * 3 + seq.sequence)
			self.writePhylip(sys.stdout)
		else:
			sys.stderr.write("Unknown format. can only print fasta or phylip.")
	def sortSamples(self, sortorder):
		# sort the samples in the alignment according to a given order
		# first check that the samples are all in the alignment
		if not set(sortorder) == set(self.samples):
			sys.stderr.write("Error: samples in sortorder are not all in the alignment. Exiting.")
			sys.exit(1)
		# then sort the samples
		self.samples = sortorder
		# and sort the sequences dict accordingly
		self.sequences = {sample: self.sequences[sample] for sample in sortorder}
	def concatenateMSA(self,msa):
		# make sure that both msa's contain the same samples
		if not set(self.samples) == set(msa.samples):
			sys.stderr.write("MSA objects must contain exactly the same set of samples for concatenation to be possible. Fix this and try again.")
			sys.exit(1)
		# let's keep track of the concatenated partitions coordinates in here as well
		if len(self.partitions) > 0:
			# then we already have the 'baseline' partition, so just add the coordinates (one-based inclusive) of the one we'll concatenate
			self.partitions.append((self.partitions[-1][1] + 1, self.partitions[-1][1] + msa.length))
		else:
			# add the baseline (first msa) partition and the one we're adding
			self.partitions.append((1,self.length))
			self.partitions.append((self.length + 1, self.length + msa.length))
		for sample in self.sequences.keys():
			self.sequences[sample].addSequence(msa.sequences[sample].sequence)
		# if there's a reference genome, concatenate those too
		if self.reference:
			self.reference.addSequence(msa.reference.sequence)
		self.checkAlnLength()
	def distanceMatrix(self):
		# function to calculate a distance matrix for all pairwise comparisons of sequences in the alignment
		distmat = {s: {} for s in self.samples}
		for i,sample1 in enumerate(self.samples):
			for j,sample2 in enumerate(self.samples):
				distmat[sample1][sample2] = {}
				if i == j:
					distmat[sample1][sample2] = 0
				else:
					# calculate distance between the two sequences
					distmat[sample1][sample2]['count'],distmat[sample1][sample2]['n_comparisons'] = self.sequences[sample1].distance(self.sequences[sample2])
					distmat[sample1][sample2]['distance'] = distmat[sample1][sample2]['count'] / distmat[sample1][sample2]['n_comparisons']
		return distmat
	def privateSubstitutions(self):
		privsubs = {s: 0 for s in self.samples}
		for i in range(self.length):
			for sample in self.samples:
				sample_allele = self.sequences[sample].sequence[i]
				other_alleles = [self.sequences[sample2].sequence[i] for sample2 in self.samples if not sample2 == sample and not self.sequences[sample2].sequence[i] in ['N','n','-','?']]
				if len(set(other_alleles)) == 0:
					continue
				elif sample_allele in other_alleles or sample_allele in ['N','n','-','?']:
					continue
				else:
					print(sample)
					print(sample_allele)
					print(other_alleles)
					privsubs[sample] += 1
					
		return privsubs					

	def writePartitions(self,outfile,format="raxml"):
		# function to write partition files, right now only tailored to raxml
		with open(outfile, "w") as of:
			for i,part in enumerate(self.partitions):
				of.write('DNA, partition_{} = {}-{}\n'.format(i,part[0],part[1]))
	# function to return a copy of self
	def copy(self):
		msa = MSA()
		msa.samples = [s for s in self.samples]
		msa.sequences = {k: Sequence(sample = k, sequence = v.sequence, meta = v.meta) for k,v in self.sequences.items()}
		msa.called_positions = self.called_positions
		msa.chrom = self.chrom
		msa.start = self.start
		msa.end = self.end
		msa.length = self.length
		msa.missingness = self.missingness
		msa.partitions = self.partitions
		msa.name = self.name
		return msa
	def subsetSamples(self, samples):
		# loop through the msa samples and delete if not in target samples
		seqs = [s for s in self.sequences.keys()]
		for s in seqs:
			if not s in samples:
				del self.sequences[s]
		# and update the sample list too
		self.samples = [s for s in self.sequences.keys() if s in samples]
		# check that all samples were found in the alignment, throw an error otherwise
		for sample in samples:
			if not sample in self.sequences.keys():
				sys.stderr.write("Error: sequence {} was not found in alignment.\n".format({sample}))
				sys.exit(1)
		# update missingness
		self.checkAlnLength()
	def subset(self,start,end,index=0):
		for seq in self.sequences.values():
			seq.subsetSequence(start,end,index)
		self.start = start
		self.end = end
		self.checkAlnLength()
	def upperCase(self):
		# change all sequences to uppercase letters
		for seq in self.sequences.values():
			seq.upperCase()
	
		# # initiate a new msa object
		# msa2 = MSA()
		# msa = MSA()
		# print(msa2.samples)
		# # metadata from msa object we can copy directly
		# msa.called_positions = self.called_positions
		# msa.chrom = self.chrom
		# msa.start = self.start
		# msa.end = self.end
		# msa.length = self.length
		# msa.missingness = self.missingness
		# msa.partitions = self.partitions
		# msa.name = self.name
		# # go through the samples in the current msa, and check if it's in samples
		# for sample in samples:
		# 	msa.addSample(name = sample, seq = self.sequences[sample].sequence, meta = self.sequences[sample].meta)
		# # and that's it
		# return msa
	def nt2aa(self):
		# function to translate all sequences in the msa to amino acids
		for sample,seq in self.sequences.items():
			print("Translating sequence {} to amino acids.\n".format(sample))
			newseq = translate(seq.sequence, break_on_stop=False)
			seq.sequence = newseq
		self.checkAlnLength()
	def removeTrailingStopCodon(self, remove_blindly=True):
		# if remove blindly, simply remove the last three bases of each sequence
		if remove_blindly:
			for seq in self.sequences.values():
				seq.sequence = seq.sequence[0:-3]
		# otherwise check if the last three make up a syop codon, and only remove if so
		else:
			for seq in self.sequences.values():
				if codon2aa(seq.sequence[-3:]) == "*":
					seq.sequence = seq.sequence[0:-3]
		# now if they differ in length fill the shorter ones with N's
		if len(set([len(seq.sequence) for seq in self.sequences.values()])) > 1:
			longest = max([len(seq.sequence) for seq in self.sequences.values()])
			for seq in self.sequences.values():
				if len(seq.sequence) < longest:
					if longest - len(seq.sequence) == 3:
						seq.sequence = seq.sequence + "N" * (longest - len(seq.sequence))
					else:
						sys.stderr.write("Error: sequences differ in length, but not by three, which they should. Please check your alignment.\n")
						sys.exit(1)
		self.checkAlnLength()
	def codonAwareFiltration(self, max_missingness, break_on_stops = True):
		# keep filtered sequences in a dict for now
		filtered = {sample: "" for sample in self.sequences.keys()}
		# keep track of internal stop codons
		internal_stops = {}
		# loop through the alignment in triplets
		for i in range(0,self.length,3):
			missing = 0
			stop_codons = 0
			for seq in self.sequences.values():
				# translate
				aa = codon2aa(seq.sequence[i:i+3])
				if aa == "?":
					missing += 1
				# check also for internal stop codons
				if aa == "*":
					stop_codons += 1
					# if there is an internal stop codon we probably want to skip this one altogether, but what we definitely want to do is set everything from here on in this sample to N, so add this info to the internal stop dict
					if not seq.sample in internal_stops.keys():
						internal_stops[seq.sample] = len(filtered[seq.sample])
			# if missingness is not too high, add this triplet to the filtered sequences
			if missing / len(self.sequences) <= max_missingness:
				for seq in self.sequences.values():
					filtered[seq.sample] += seq.sequence[i:i+3]
		# check stop codons
		if len(internal_stops.keys()) > 0:
			#print(internal_stops)
			# if we break on stop codons
			if break_on_stops:
				sys.stderr.write("Internal stop codon(s) found in sample(s) {}. Please check your alignment.\n".format(",".join(internal_stops.keys())))
				sys.exit(1)
			#otherwis filter out everything after the first stop codons
			else:
				for sample,position in internal_stops.items():
					# get the two pieces
					p1 = filtered[sample][0:position]
					p2 =  "N" * (len(filtered[sample]) - position)
					#print(p1)
					filtered[sample] = p1 + p2
					#print("Pos {}-{} will be converted to N in sample {}.".format(position,self.length,sample))
					#print("N" * (self.length - position))
					# and print a warning about this
					sys.stderr.write("Internal stop codon found in sample {} at position {}. All bases after this one will be converted to N for this sample.\n".format(sample, position))
		# then replace the sequences with the filtered ones
		for sample,seq in self.sequences.items():
			self.sequences[sample].sequence = filtered[sample]
		# and check that the lengths are still ok
		self.checkAlnLength()
	def renameSequences(self, names):
		# takes a dictionary where keys are old names and values are new names
		# loop through the sequences and rename them
		for old,new in names.items():
			self.sequences[new] = self.sequences.pop(old)
			self.sequences[new].sample = new






				




# this one gets chromosomes and their lengths from a vcf file
def getChromLengths(vcf, only_contigs_with_records=True):
	all_contigs = {contig: 0 for contig in list(vcf.header.contigs)}
	# fetch all the chromosomes to begin with from the vcf header
	for contig in all_contigs.keys():
		length = vcf.header.contigs[contig].length
		all_contigs[contig] = length
	contigs = {}
	if only_contigs_with_records:
		# usually we'll only want the chromosome lengths for the chromosome(s) that has records in this vcf file
		for contig,length in all_contigs.items():
			try: 
				# check this by trying to fetch the first record on each chromosome, if successfull, add this contig to the returned dict
				next(vcf.fetch(contig))
				contigs[contig] = (0,length)
			except:
				pass
	else:
		contigs = all_contigs 
	return contigs


# and this one splits the output from previous function into windows
def generateWindows(intervals, window, step_size):
	print(intervals)
	# first, if the intervals are given as a simle chrom: length dictionary, add 0 as start position to each such interval
	if not type(intervals[list(intervals.keys())[0]]) is tuple:
		for c in intervals.keys():
			intervals[c] = (0, intervals[c])
    #create a dictionary with windows for each chromosome represented in intervals
	windows = []
	window_number = 1
	start = intervals[list(intervals.keys())[0]][0]
	end = start + window
	for c in intervals.keys():
		start = intervals[c][0]
		length = intervals[c][1]
		chrom = c
		while start < length:
			windows.append({'window_number': window_number, 'chrom': chrom, 'start':start, 'end':end})
			start = start + step_size
			end = start + window if start + window < length else length
			window_number = window_number + 1
		else:
			break
	return windows

# function to parse gff file into a list of features
def parseGFF(gff):
	features = []
	with open(gff) as f:
		for line in f.readlines():
			if line.startswith("#"):
				continue
			else:
				name,source,feature,start,end,score,strand,frame,attribute = line.strip().split("\t")
				# extract gene name/ID if possible
				splitattr = [re.split(' |=', i.strip().replace('"','')) for i in attribute.split(";") if not i == ""]
				attrDict = {i[0].lower(): i[1] for i in splitattr}
				if 'gene_name' in attrDict.keys():
					gene_name = attrDict['gene_name']
				elif 'name' in attrDict.keys():
					gene_name = attrDict['name']
				elif 'gene_id' in attrDict.keys():
					gene_name = attrDict['gene_id']
				else:
					gene_name = None
				if not feature == 'source':
					splitattr = zip(attribute.split(";"))
					features.append({
						'seq': name,
						'source': source,
						'feature_type': feature,
						'start': int(start) - 1,
						'end': int(end),
						'score': score,
						'strand': strand,
						'frame':frame,
						'attribute':attribute,
						'gene_name': gene_name,
					})
	return features


# just a wrapper to fetch an input file by determining file format
def readSequenceFile(file):
	if file.endswith((".fa",".fasta",".fa.gz",".fasta.gz")):
		return MSA.fromFasta(file)
	elif file.endswith((".phy",".phylip",".phylip.gz",".phy.gz")):
		return MSA.fromPhylip(file)

# and a similar one for writing the output
def writeSequenceFile(msa, file):
	if file.endswith((".fa",".fasta",".fa.gz",".fasta.gz")):
		return msa.writeFasta(file)
	elif file.endswith((".phy",".phylip",".phylip.gz",".phy.gz")):
		return msa.writePhylip(file)
	else:
		sys.stderr.write("Error: unknown file format. Exiting.")
		sys.exit(1)

# function to parse a two column popfile into a dictionary
def parsePopfile(popfile):
	pops = {}
	with open(popfile) as f:
		for line in f.readlines():
			if line.strip().startswith("#"):
				# header lines starting with hash should be ignored
				continue
			if line.strip() == "":
				# we of course also want to skip empty lines (e.g. trailings)
				continue
			try:
				sample,pop = line.strip().split()
				if pop in pops.keys():
					pops[pop].append(sample)
				else:
					pops[pop] = [sample]
			except:
				# it's ok not to have a sample assigned, but print this as warning
				sys.stderr.write("Sample {} was not assigned to any population")
	return pops