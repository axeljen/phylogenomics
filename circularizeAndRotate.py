import functions as fn
import sys
import argparse
import os
import re

#function that will build a numerical list from a sequence and kmercounts, to identify duplicate regions
def kmerProfile(seq, kmercounts):
	#we'll only look at kmers appearing > 1, so starting of with a list of 1s of the same length as the sequence
	seqlen = len(seq)
	counts = [1 for i in range(0, seqlen)]
	#then find the sequences which kmers are non-unique
	for kmer,i in kmercounts.items():
		if i > 1:
			matches = re.finditer(kmer, seq)
			for match in matches:
				start = match.span()[0]
				end = match.span()[1]
				for y in range(start, end):
					if counts[y] < i:
						counts[y] = i
	return counts

#this function will check for duplicated (or more) kmers at the beginning and end of sequence, and cut off the end sequence if matching
def CircularizeGenome(seq, counts):
	subseq = ""
	for i in range(0,len(counts)):
		if counts[i] > 1:
			subseq = subseq + seq[i]
		else:
			break
	#double check so that it matches with the ending string as well
	if not subseq == "":
		match = re.finditer(subseq, seq)
		s = {'start':0,'end':0}
		for m in match:
			s['start'] = m.span()[0]
			s['end'] = m.span()[1]
		if s['end'] == len(seq) and s['start'] < s['end']:
			circularized_sequence = seq[0:s['start']]
			sys.stderr.write("Found evidence of circularization, cutting the following {} bp trail from the end:\n{}\n".format(s['end'] - s['start'], seq[s['start']:s['end']]))
		else:
			circularized_sequence = None
	else:
		circularized_sequence = None
		sys.stderr.write("Didn't find any signs of circularization.\n")
	if not circularized_sequence:
		# if circularization failed, check if there are duplications near the edges, and print a warning about this in that case
		first_dup = 0
		last_dup = 0
		for i,c in enumerate(counts):
			if c > 1:
				first_dup = i
				break
		counts_rev = counts[::-1]
		for i,c in enumerate(counts_rev):
			if c > 1:
				last_dup = i
				break
		if first_dup < 50 or last_dup < 50:
			sys.stderr.write("Warning: there are duplicated sequences close to the beginning or end of sequence, allthough I didn't manage to circularize it. Might be something to check manually!\n")
			sys.stderr.write("Duplicated kmer starting {} bp from beginning, and {} bp from end.\n".format(first_dup,last_dup))
	if not circularized_sequence:
		circularized_sequence = seq
		circularization = False
	else:
		circularization = True
	return circularized_sequence,circularization

# function to find the corresponding starting position from a blast results dict
def findStartPos(blastres):
	positions = [] #will store blast hits as tuples, (querystart,subjectstart)
	for hit in blastres.values():
		positions.append((int(hit['q. start']), int(hit['s. start'])))
	if len(positions) > 1:
		#find the first starting position on the query seq, and use the equivalent on the subjectseq as startpoint
		qpos = None
		for q,s in positions:
			if qpos == None or q < qpos:
				startpos = s
				qpos = q
	else:
		startpos = positions[0][1]
	return startpos

parser = argparse.ArgumentParser()

parser.add_argument('-i', '--input', type=str, help='Input sequence to search for circularization and rotate.', required=True)
parser.add_argument('-o', '--output', type=str, help="Output file.")
parser.add_argument('--force-rotation', action='store_true', help="Rotate even if evidence for circularization could not be found. Can be ok if genome already is circularized, for example, but better be sure about that.")
parser.add_argument('-k', '--kmer-size', type=int, default=13)

parser.add_argument('-ref', '--reference', type=str, help="File after which to align the starting position.")

args = parser.parse_args()

# read sequence input
sys.stderr.write("Reading input file:\n\{}\n".format(args.input))
msa = fn.readSequenceFile(args.input)

# loop through the sequences in msa and count kmers
kmers = []
for n,seq in msa.sequences.items():
	kmers.append(seq.countKmers(13))
	prof = kmerProfile(seq.sequence, kmers[-1])
	circ,is_circularized = CircularizeGenome(seq.sequence, prof)
	if is_circularized == False:
		# if not circularized, default is to quit
		if not args.force_rotation:
			sys.stderr.write('Circularization was not found, enforce --force-rotation to carry on with shifting of starting position anyway.\n\n')
			sys.exit()
		else:
			sys.stderr.write("Warning: Proceeding with rotation of genome in the absence of circularization due to --force-rotation invoked.\n")
	seq.sequence = circ
	msa.checkAlnLength()
	# blast the reference sequence to get the starting position
	blast = seq.blastToMe(args.reference)
	# fetch the starting position of the lowest pos blast hit
	startpos = findStartPos(blast)
	# now, we simply cut the part 1-startpos, and paste this to the end of the sequence
	rotseq = seq.sequence[startpos - 1:] + seq.sequence[0:startpos - 1]
	seq.sequence = rotseq
	msa.checkAlnLength()

# write circularized and rotated sequence file
fn.writeSequenceFile(msa,args.output)

sys.stderr.write("\nShifted starting position to match {}. Rotated genome written here:\n{}\n\n\n".format(args.reference,args.output))