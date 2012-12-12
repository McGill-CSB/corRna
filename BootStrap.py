#! /usr/bin/env python
import sys
import re
import commands
import shutil
import getopt
import os, glob, sys
import os.path
import random, math
from multiprocessing import Process, Manager, Array, Value
from decimal import Decimal
import time
from math import *
import Pearson

DEBUG = False
class bootstrap:
	def __init__(self, wildtype, depth):
		self.wild = wildtype
		self.depth = depth
		if DEBUG:
			print "to grab correlation"
		self.corr = Pearson.Correlation(self.wild)
		self.boot = Manager().list()
		if DEBUG:
			print "done grab correlation"
	def start(self):
		if DEBUG:
			print "ready to start"
		self.startSampling()
	def getStanding(self,*args):
		if len(args) == 1:
			seq = args[0]
			depth = self.detectDepth(seq)
			if depth == 0:
				return "---"
			sampling = self.boot[depth]
			corr = self.corr.getCorr(seq)
			sampling = self.boot[int(depth)-1]
			sampling.append(corr)
			sampling.sort()
			return float(float(sampling.index(corr))/float(len(sampling)))
		elif len(args) == 2:
			depth = self.detectDepth(args[0])
			if depth == 0:
				return "---"
			sampling = self.boot[int(depth)-1]
			sampling.append(args[1])
			sampling.sort()
			return float(float(sampling.index(args[1]))/float(len(sampling)))
	def detectDepth(self,seq):
		counter = 0
		for i in range(len(self.wild)):
			if(self.wild[i] != seq[i]):
				counter+=1
		return counter
	def startSampling(self):
		proc = []
		for i in range(1,int(self.depth)+1):
			if DEBUG:
				print i
			p = Process(target=self.generate, args=(i,))
			p.start()
			proc.append(p)
		for p in proc:
			p.join()
	def generate(self, depth):
		kSequence = []
		for i in range(1000):
			kSequence.append(self.computeSeq(depth))
		self.boot.insert(int(depth)-1,kSequence)
	def computeSeq(self,depth):
		seq = self.sequenceRandomizer(self.wild, depth)
		corr = self.corr.getCorr(seq)
		return corr
			
	def sequenceRandomizer(self, wildtype, n):
		random.seed()
		for i in range(n):
			x = random.randint(1,len(wildtype)-1)
			wildtype = self.substring(wildtype,x,self.ranRNA(wildtype[x]))
		return wildtype
	def substring(self, stri, posi, ch):
		return stri[:posi]+ch+stri[posi+1:]
	def ranRNA(self,c):	
		random.seed()	
		x = random.randint(0,3)	
		rna = ["A","C","G","U"]
		if(rna[x] == c): 
			return self.ranRNA(c)
		return rna[x]
		
