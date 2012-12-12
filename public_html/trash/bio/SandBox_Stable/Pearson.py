#!/usr/bin/env python
import sys
import os.path
import re
import random, math
import commands
from math import *
from decimal import Decimal

DEBUG = True

class basePairProbability:
	#constructor
	def __init__(self, sequence, entry):
		self.entries = []
		self.bpp = []
		self.seq = sequence
		self.sequenceSize = 0	
		self.box = entry

	def parse(self):
		self.sequenceSize = len(self.seq)
		##initialize bpp
		for i in range(self.sequenceSize):
			self.bpp.append(float(0))
		for line in self.box:
			data = re.split('\s+', line)
			entry = []
			entry.append(int(data[0]))
			entry.append(int(data[1]))
			entry.append(float(data[2]))
			self.entries.append(entry)
			self.bpp[entry[0]-1]+=entry[2]
			self.bpp[entry[1]-1]+=entry[2]
class pearson: 
	#onstructor 
	def __init__(self, typeA, typeB):
		self.X = typeA
		self.Y = typeB
		self.r = 0
	def compute(self):
		sum_x = 0
		sum_y = 0
		sum_XX = 0
		sum_YY = 0
		sum_XY = 0
		for i in range(len(self.X.bpp)):
			x = self.X.bpp[i]
			sum_x+=x		
			y = self.Y.bpp[i]
			sum_y+=y
		p = sum_x/(len(self.X.bpp))		#mean of x
		q = sum_y/(len(self.Y.bpp))		#mean of y

		for i in range(len(self.X.bpp)):
			x = self.X.bpp[i] - p
			y = self.Y.bpp[i] - q	
			XX = pow(x,2)
			YY = pow(y,2)
			sum_XX+=XX
			sum_YY+=YY
			sum_XY+=(x*y)
		self.r = sum_XY / (sqrt(sum_XX*sum_YY))  ##pearson's coefficient correlation
##end of class correlation

class Correlation:
	def __init__(self, WildSeq):
		self.Wild = WildSeq
		w_fname = self.RNAfold(self.Wild)
		w_entry = self.getbpp(w_fname)
		self.WildBpp = basePairProbability(self.Wild, w_entry)
		self.WildBpp.parse()

	def getCorr(self, seq):
		self.Seq = seq
		fname = self.RNAfold(self.Seq)
		seq_entry = self.getbpp(fname)
		self.seqBpp = basePairProbability(self.Seq, seq_entry)
		self.seqBpp.parse()
		correlation = pearson(self.WildBpp, self.seqBpp)
		correlation.compute()
		self.corr = correlation.r
		return self.corr	
	def getbpp(self, fname):
		f = open(fname+'_dp.ps','r')
		entry = []
		for line in f:
			line.strip()
			if(re.search('\d+.*ubox',line)):
				line = re.sub('\\n','',line)
				entry.append(line)
		commands.getoutput('rm '+fname+'*.ps')
		return entry
		#computes RNAfold
	def RNAfold(self, seq):
		rant = random.randint(1,100)
		rant_t = random.randint(1,100)
		fname = str(rant)+'-'+str(rant_t)
		if os.path.isfile(fname+'.tmp'):
			return RNAfold(self, seq)
		f = open(fname+'.tmp','w')
		f.write('> '+fname+'\n')
		f.write(seq)
		f.close()
		out = commands.getoutput('cat '+fname+'.tmp @- | RNAfold -p -d0')		
		commands.getoutput('rm '+fname+'.tmp')
		return fname
