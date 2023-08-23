import os
import re
import argparse
import pysam

def openFile(filepath, mode='r', append=False):
    if filepath.endswith(".gz"):
        if append:
            return bgzf.open(filepath, 'rb')
        else:
            return bgzf.open(filepath, 'wt')
    else:
        if append:
            return open(filepath, 'a')
        else:
            return open(filepath, 'w')

def closeFile(f, filepath):
    if filepath.endswith(".gz"):
        return bgzf.close(f)
    else:
        return f.close()


def WriteSeqDictToFasta2(seqDict, outfile, append=False):
    with openFile(outfile) as f:
        for name,seq in seqDict.items():
            f.write('>' + name + "\n" + seq + "\n")
        closeFile(f, outfile)


#UPDATED SEQDICT FORMAT function to handle subsetting of fasta file with pysam, returning a seqdict
def FastaFetch2(fastafile, sequences=None, start=None, stop=None):
    f = pysam.FastaFile(fastafile)
    seqDict = {}
    if sequences is not None:
        if os.path.isfile(sequences):
            seqList = []
            with open(sequences) as s:
                for line in s.readlines():
                    if not line == "":
                        seqList.append(line.rstrip().lstrip(">")) 
        else:
            seqList = []
            for s in sequences.split(","):
                if not s == "":
                    seqList.append(s.rstrip().lstrip(">"))
        if start and stop:
            for seq in seqList:
                seqDict[seq] = f.fetch(seq, start, stop)
        elif start:
            for seq in seqList:
                seqDict[seq] = f.fetch(seq, start)
        else:
            for seq in seqList:
                seqDict[seq] = f.fetch(seq)
    else:
        if start and stop:
            seqDict = []
            for seq in f.references:
                seqDict[seq] = f.fetch(seq, start, stop)
        elif start:
            for seq in f.references:
                seqDict[seq] = f.fetch(seq, start)
        else:
            for seq in f.references:
                seqDict[seq] = f.fetch(seq)
    return seqDict

def parseGFF(gff):
	features = []
	with open(gff) as f:
		for line in f.readlines():
			feature = [i.rstrip().lstrip() for i in line.split("\t")]
			seqname = feature[0]
			source = feature[1]
			feat = feature[2]
			start = feature[3]
			end = feature[4]
			score = feature[5]
			strand = feature[6]
			frame = feature[7]
			attribute = feature[8]
			features.append({
				'seqname':seqname,
				'source':source,
				'feature':feat,
				'start':int(start),
				'end':int(end),
				'score':score,
				'strand':strand,
				'frame':int(frame) if not frame == "." else frame,
				'attribute':attribute
			})
	return features

def subsetSequence(sequence, start, end, index=1):
	subsequence = sequence[start - index:end]
	return subsequence

def complementery_sequence(seq, reverse=True):
    complements = {
        'A':'T',
        'G':'C',
        'T':'A',
        'C':'G',
    }
    seq_complement = ""
    if reverse:
        while len(seq_complement) < len(seq):
            original_base = seq[len(seq) - (len(seq_complement) + 1)]
            complement_base = complements[original_base]
            seq_complement = seq_complement + complement_base
    else:
        while len(seq_complement) < len(seq):
            original_base = seq[len(seq_complement)]
            complement_base = complements[original_base]
            seq_complement = seq_complement + complement_base
    return seq_complement


parser = argparse.ArgumentParser()

#arguments for in and output:
parser.add_argument("-i", "--input", help="Input fasta-file.", required=True)
parser.add_argument("-o", "--output-directory", help="Directory to put the extracted features in", required=True)
parser.add_argument("-r", "--regions", help="Regions to extract in tab-separated format. Minimum two colums: start, end. Third column can be given to name the resulting outfile.")
parser.add_argument("--gff", help="GFF-annotation file can be given as well, will be parsed and named according to feature name.")
parser.add_argument("--incorporate-strand-info", action="store_true")
args = parser.parse_args()

input_fasta=args.input
gff=args.gff

seqDict = FastaFetch2(input_fasta)

#check so that there's a single sequence in fasta and get the name
if len(seqDict) > 1:
	print("Can only do this on single fasta sequences")
	sys.exit()

#get sample name
sample_name = list(seqDict.keys())[0]

#parse annotations
features = parseGFF(gff)

#will store all extracted features in this featuredict
featureDict = {}

for feature in features:
	name = feature['attribute'].split("=")[1].split(" ")[0]
	feat_type = feature['feature']
	if feat_type in ['tRNA','rRNA','CDS']:
		subseq = subsetSequence(seqDict[sample_name], feature['start'], feature['end'])
		if feature['strand'] == "-":
			subseq = complementery_sequence(subseq)
		#add Ns to get first coding position to be first position
		Ns = "N" * int(feature['frame'])
		featureDict[sample_name + "@" + name] = Ns + subseq

#write each feature to a separate file
for name,sequence in featureDict.items():
	outfile = os.path.join(args.output_directory, name + ".fa")
	WriteSeqDictToFasta2({name: sequence}, outfile)