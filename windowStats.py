#!/usr/bin/python
import sys
import os
import argparse
import random
import functions as fn
from multiprocessing import Process
from multiprocessing import SimpleQueue
from threading import Thread
import time
import datetime
from datetime import timedelta
from time import sleep
import tempfile
from datetime import datetime
import numpy as np
import pysam

# class to keep store window coordinates, metadata and results
class windowStats:
	def __init__(self,samples,chrom,start,end):
		self.chrom = chrom
		self.start = start
		self.end = end
		self.firstpos = None
		self.lastpos = None
		self.samples = samples
		self.sites = 0
		self.results = {sample: {'hetsites':0,'heterozygosity':0,'homsites':0,'calledsites':0,'missingsites':0,'missing_fraction':0,'mean_depth':0,'median_depth':0,'depth':0} for sample in samples}
	def parse_rec(self,rec):
		# function to read a single vcf record and add to window
		if self.firstpos == None:
			self.firstpos = rec.pos
		self.lastpos = rec.pos
		self.sites += 1
		# loop through all samples and fetch stats
		for sample,value in rec.samples.items():
			gt,dp = value['GT'],value['DP']
			if None in gt:
				# if record is missing in sample, continue
				self.results[sample]['missingsites'] += 1
				continue
			if dp:
				self.results[sample]['depth'] += dp
			self.results[sample]['calledsites'] += 1
			if gt[0] != gt[1]:
				# if heterozygous genotype
				self.results[sample]['hetsites'] += 1
				continue
			else:
				self.results[sample]['homsites'] += 1
	def iterate(self,vcf):
		for rec in vcf.fetch(self.chrom,self.start,self.end):
			self.parse_rec(rec)
	# function to summarize window once it's finished
	def sumwindow(self):
		for sample in self.samples:
			if self.results[sample]['calledsites'] == 0:
				self.results[sample] = {k: 'NA' for k in self.results[sample].keys()}
				continue
			self.results[sample]['heterozygosity'] = self.results[sample]['hetsites'] / self.results[sample]['calledsites']
			self.results[sample]['missing_fraction'] = self.results[sample]['missingsites'] / (self.results[sample]['calledsites'] + self.results[sample]['missingsites'])
			self.results[sample]['mean_depth'] = self.results[sample]['depth'] / self.results[sample]['calledsites']

# function to analyze windows in parallell
def analyseWindow(windowQueue,invcf,tmpout):
	while True:
		if _FINISH:
			break
		window = windowQueue.get()
		# open vcf file
		vcf = pysam.VariantFile(invcf)
		window.iterate(vcf)
		window.sumwindow()
		with open(tmpout, 'a') as tf:
			for sample in window.samples:
				tf.write('{chrom}\t{start}\t{end}\t{firstpos}\t{lastpos}\t{sites}\t{sample}\t{het}\t{depth}\t{missfrac}\t{called}\t{missing}\n'.format(
					chrom=window.chrom,
					start=window.start + 1,
					end=window.end,
					firstpos=window.firstpos,
					lastpos=window.lastpos,
					sites=window.sites,
					sample=sample,
					het=window.results[sample]['heterozygosity'],
					depth=window.results[sample]['mean_depth'],
					missfrac=window.results[sample]['missing_fraction'],
					called=window.results[sample]['calledsites'],
					missing=window.results[sample]['missingsites']
				))
		global windows_finished
		windows_finished +=1

# indefinite loop to write some stats every tenth second as the program runs:
def checkStats(tmpfile, multiplier, totwindows,report_interval=5):
	while True:
		finished = int(len([line for line in open(tmpfile).readlines() if not line.strip() == ""]) / multiplier)
		sleep(report_interval)
		global windows_finished
		global analysis_started
		windows_finished = finished
		if analysis_started:
			sys.stderr.write("{time}: Analysed {finished} / {total} windows ({percentage} %).\n".format(time=time.ctime(time.time()),finished=windows_finished,total=totwindows,percentage=(windows_finished / totwindows) * 100))
		if _FINISH:
			return

# function to fetch the tempfile when finished to a dictionary, for clean outputwriting
def fetchResults(tmpout):
	tmp = open(tmpout.name)
	windows_read = []
	results = {}
	for line in tmp.readlines():
		windString = "{chrom}_{start}_{end}".format(chrom=line.split("\t")[0],start=line.split("\t")[1],end=line.split("\t")[2])
		if not windString in windows_read:
			windows_read.append(windString)
			chrom,start,end,firstpos,lastpos,sites,sample,het,depth,missfrac,called,missing = line.strip().split("\t")
			results[windString] = {'chrom': line.split()[0], 'start': line.split()[1], 'end': line.split()[2], 'results':{}}
			results[windString]['results'][sample] = {
				'firstpos':firstpos,
				'lastpos':lastpos,
				'sites':sites,
				'heterozygosity':het,
				'depth':depth,
				'missing_fraction':missfrac,
				'called_sites':called,
				'missing_sites':missing
				}
		else:
			chrom,start,end,firstpos,lastpos,sites,sample,het,depth,missfrac,called,missing = line.strip().split("\t")
			windString = "{chrom}_{start}_{end}".format(chrom=chrom,start=start,end=end)
			results[windString]['results'][sample] = {
				'firstpos':firstpos,
				'lastpos':lastpos,
				'sites':sites,
				'heterozygosity':het,
				'depth':depth,
				'missing_fraction':missfrac,
				'called_sites':called,
				'missing_sites':missing
				}
	return results

# function to write output file
def writeFinalResults(results,windows,outfile):
		with open(args.out, "w") as of:
			# write header
			of.write('\t'.join(['chrom','start','end','firstpos','lastpos','sites','sample','heterozygosity','mean_dp','missing_fraction','called_sites','missing_sites']) + "\n")
			# loop through the windows to get the results in the correct order
			for window in windows:
				windString = "{}_{}_{}".format(window['chrom'],window['start'] + 1,window['end'])
				windres = results[windString]
				chrom = windres['chrom']
				start = windres['start']
				end = windres['end']
				for sample,res in windres['results'].items():
					of.write('\t'.join([chrom,str(start),str(end)]) + '\t')
					of.write('\t'.join([str(res['firstpos']),str(res['lastpos']),str(res['sites']),sample,str(res['heterozygosity']),str(res['depth']),str(res['missing_fraction']),str(res['called_sites']),str(res['missing_sites'])]) + "\n")


# input arguments
parser = argparse.ArgumentParser(description="Collect some individual stats from genomic windows.")

parser.add_argument('-v','--vcf', type=str, help='Input vcf')

parser.add_argument('-o','--out', type=str, help='Path to write output tsv file to.')

parser.add_argument('-w','--window-size',type=int, help='Window size.')
parser.add_argument('-s','--step-size',type=int,help="Step size, defaults to window size (non-overlapping).")

parser.add_argument('--samples', type=str, help="File listing the samples to analyze, one per line.", default=None)

parser.add_argument('-T', '--threads', type=int, help="Specify how many threads to run.", default=1)

args = parser.parse_args()

# get step size manually if not given
if not args.step_size:
	args.step_size = args.window_size

# fetch samples either from argparser or vcf header
	samples = []
	if args.samples:
		for line in open(args.samples).readlines():
			if not line.strip() == "":
				samples.append(line.strip())
	else:
		args.samples = list(pysam.VariantFile(args.vcf).header.samples)
		samples = args.samples

if __name__ == "__main__":
	analysis_started = False
	# write time and parameters
	starttime = time.time()
	sys.stderr.write('{t}: Running windowStats.py with the following input parameters:\n'.format(t=time.ctime(starttime))) 
	sys.stderr.write('\n' + "#" * 100 + '\n\n' + '\n'.join(f'{k}: {v}' for k, v in vars(args).items()) + "\n\n" + "#" * 100 + "\n\n")
	# start some variables for tracking
	_FINISH=False
	windows_queued = 0
	global windows_finished
	windows_finished = 0

	chromlengths = fn.getChromLengths(pysam.VariantFile(args.vcf))
	windows = fn.generateWindows(chromlengths,args.window_size,args.step_size)

	# store all the window results in a list
	results = []

	# traverse through the windows and get the stats
	threads = args.threads
	tmpout = tempfile.NamedTemporaryFile(mode="w", prefix="temp.{}".format(os.path.basename(args.out)), suffix=".txt", dir = ".", delete=False)
	#list to put all processes in:
	workerprocesses = []

	windowQueue = SimpleQueue()

	# start threads, analyses will not start until we feed it data
	for thread in range(threads):
		worker =Process(target=analyseWindow, args=(windowQueue,args.vcf,tmpout.name))
		worker.daemon = True
		worker.start()
	analysisStarted = True

	
	stats = Thread(target=checkStats, args=(tmpout.name,len(samples), len(windows)))
	stats.daemon = True
	stats.start()

	start_time = time.time()
	analysis_started = True
	sys.stderr.write('Starting to analyze {w} windows on {t} threads.\n\n'.format(w=len(windows),t=threads)) 

	# loop through the windows and add them to the queue
	for window in windows:
		windowQueue.put(windowStats(samples, chrom = window['chrom'],start = window['start'], end = window['end']))
		windows_queued += 1

	while windows_finished < len(windows):
		sleep(5)
	_FINISH = True

	sys.stderr.write('{t}: Finished with analysis, sorting results and writing output.\n\n'.format(t = time.ctime(time.time()))) 
	# grab the tempfile output to an (unsorted still) dictionary when finished
	results = fetchResults(tmpout)

	# last, loop through the window dictionary and write the corresponding results to get it in order
	writeFinalResults(results,windows,args.out)

	# remove temfile
	os.remove(tmpout.name)

	sys.stderr.write('{t}: All done. Elapsed time: {t2}.'.format(t=time.ctime(time.time()), t2 = timedelta(seconds=time.time() - starttime)))
