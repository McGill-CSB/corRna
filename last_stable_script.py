#! /usr/bin/env python
import sys
import re
import commands
import shutil
import getopt
import os, glob, sys
import os.path
import random, math
import time
from multiprocessing import Process, Manager, Array, Value
from decimal import Decimal
from time import clock, time
from math import *
### note ### we are assumming the squence length of typeA and typeB are the same!

DEBUG = False
RANDOM = random.randint(1,100)
DUMMY = 0
start = time()
manager = Manager()

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
class sampling:
	def __init__(self, sequence, SETR, SETF):
		self.wild = sequence
		if RNAMUTANT:
			self.setR = SETR
		if FIXEDCGSAMPLING:
			self.setF = SETF
		self.allSeq = []
		#collect wild type bpp
		self.wild_entry = []
		fname = RNAfold(self.wild)
		self.wild_entry = getbpp(fname)
		self.wild_bpp = basePairProbability(self.wild, self.wild_entry)
		self.wild_bpp.parse()
		#print self.wild_bpp.bpp
	def compute(self):	
		entry = manager.list()
		if(DEBUG):
			print "start grabbing corr"
		proc = []
		if RNAMUTANT:
			for seq in self.setR:
				p=Process(target=self.split, args=(seq,'R',entry))
				p.start()
				proc.append(p)
		if FIXEDCGSAMPLING:
			for seq in self.setF:
				p=Process(target=self.split, args=(seq,'F',entry))
				p.start()
				proc.append(p)
		for p in proc:
			p.join()
		if(DEBUG):
			print "finish grabbing corr"
		self.allSeq = entry
	def split(self, seq, st, entry):
		if(re.search('F',st)):
			sequence = seq[0]
		else:
			sequence = seq
		tuples = []
		fname = RNAfold(sequence)
		seq_entry = getbpp(fname)	
		seq_bpp = basePairProbability(sequence, seq_entry)
		seq_bpp.parse()
		correlation = pearson(self.wild_bpp, seq_bpp)
		correlation.compute()	
		corr = correlation.r
		tuples.append(st)
		tuples.append(sequence)
		tuples.append(corr)
		if(re.search('F',st)): #append mutation Depth 
			tuples.append(seq[1])
		else: #get mutation Depth
			tuples.append(self.getMDepth(sequence))
		entry.append(tuples)
	def getMDepth(self,seq):
		mDepth = 0
		for i in range(len(self.wild)):
			if(self.wild[i] != seq[i]):
				mDepth+=1
		return mDepth
def getbpp(fname):
	f = open(fname+'_dp.ps','r')
	entry = []
	for line in f:
		line.strip()
		if(re.search('\d+.*ubox',line)):
			line = re.sub('\\n','',line)
			entry.append(line)
	f.close()
	commands.getoutput('rm '+fname+'*.ps')
	return entry
	#computes RNAfold
def RNAfold(seq):
	#print "> " +str(seq)
	rant = random.randint(1,99)
	rant_t = random.randint(1,99)
	fname = str(rant)+'-'+str(rant_t)
	f = open(fname+'.tmp','w')
	f.write('> '+fname+'\n')
	f.write(seq)
	f.close()
	out = commands.getoutput('cat '+fname+'.tmp @- | RNAfold -p -'+DANGLING)		
	try:
		f = open(fname+'_dp.ps', 'r')
		f.close()
	except IOError, e:
		if not re.match('',seq):
			tmp = seq
			tmp = re.sub('A|G|C|U|a|g|c|u', '', tmp)
			if re.match('', tmp):
				time.sleep(1)	
				try: 
					f = open (fname+'_dp.ps', 'r')
					f.close()
				except IOError, e:
					return RNAfold(seq)
	commands.getoutput('rm '+fname+'.tmp')
	return fname
	#computes fixedCGSampling 
#end of class sampling
class bootstrap:
	def __init__(self,wildtype, size):
		self.wild = wildtype
		self.seq = []
		self.n = size
		self.bootstrap = manager.list()
		self.wild_entry = []
		fname = RNAfold(self.wild)
		self.wild_entry = getbpp(fname)
		self.wild_bpp = basePairProbability(self.wild, self.wild_entry)
		self.wild_bpp.parse()
		##start bootstrap
		self.strap()
	def strap(self):
		proc = []
                for i in range(1,self.n+1):
			if(DEBUG):
                        	print i
                        p = Process(target=self.sequenceGenerator, args=(self.wild, i,))
                        p.start()
                        proc.append(p)
                for p in proc:
                        p.join()
	def sequenceGenerator(self, wildtype, n):
		entry = manager.list()
		kSequence = []	
		proc = []
		for i in range(1000):	
			p = Process(target=self.computeSeq, args=(wildtype, n, entry))
			p.start()
			p.join()
		kSequence = entry
		self.bootstrap.insert(n, kSequence)
		#return kSequence
	def computeSeq(self, wildtype, n, entry):
		pair = [] 
		seq = self.sequenceRandomizer(wildtype, n)
		fname = RNAfold(seq)
		seq_entry = getbpp(fname)
		if(len(seq_entry) == 0):
			return self.computeSeq(wildtype, n, entry) 
		else:
			seq_bpp = basePairProbability(seq, seq_entry)
			seq_bpp.parse()
			correlation = pearson(self.wild_bpp, seq_bpp)
			correlation.compute()	
			corr = correlation.r
			entry.append(corr)
		##compute corre,ation
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
class filterSeq:
	def __init__(self, boot, allSequence):
		self.strap = boot
		self.allSeq = allSequence	
		self.start()
	def start(self):
		out = sorted(self.allSeq, self.compare)
		self.result(out)
	def result(self,out):
		for s in out:
			ext = ''
			#if(re.search(s[0],'F')):
			ext = self.getStanding(float(s[2]),int(s[3]))
			#else: 
			#	ext = "---"
			if(DEBUG):
				print str(s) + " bootstrap: " + str(ext)	
			else:
				for i in range(len(s)):
					print str(s[i])+", ",
				print str(ext)+"<br>"	
	def getStanding(self, corr, lvl):
		if(lvl == 0):
			return "---"
		sample = self.strap[lvl-1]
		sample.append(corr)
		sample.sort()
		return float(float(sample.index(corr))/float(len(sample)))
	def compare(self,a, b):
		return cmp(float(b[2]), float(a[2]))
			
class collect_Sequence:
	def __init__(self, sequence):
		self.seq = sequence
		RNA_d = manager.list()
		fixedCG_d = manager.list()
		MA = Value('i',0)
		if(RNAMUTANT):
			rnamut =Process(target=self.RNAmutants, args=(self.seq,RNA_d,MA))
		else:
			self.RNA_entry = 0
		if(FIXEDCGSAMPLING):
			fixedCG = Process(target=self.fixedGC, args=(self.seq,fixedCG_d,MA))
		else:
			self.fixedCG_entry = 0
		if(RNAMUTANT):
			rnamut.start()
		if(FIXEDCGSAMPLING):
			fixedCG.start()
			fixedCG.join()
		if(RNAMUTANT):
			rnamut.join()
		if(DEBUG):
			print str(time()-start)+" Completed: Collecting Seq"
		if(RNAMUTANT):
			self.RNA_entry = RNA_d
		if(FIXEDCGSAMPLING):
			self.fixedCG_entry = fixedCG_d
		self.MAX = MA.value
	def GCcontent(self,sequence):
		gc=0;
		l = 0;
		s = re.compile("C|G|c|g")
		for c in sequence:
			l+=1
			if(s.match(c)):
				gc+=1
		return (Decimal(str(gc)) / Decimal(str(l))).quantize(Decimal("0.0001"))
	## retrives the base pair probability
	def fixedGC(self,seq,fixedCG_d, MA):
		gc = self.GCcontent(seq)
		out = commands.getoutput('./fixedCGSampling.py '+seq+' -n '+str(TOP)+' -g '+str(gc)+' -e '+str(ERROR)+' -m '+str(MUTATION_DEPTH))
		##parse the data out
		flag = 0
		num = ''
		pair = []
		tmp = ''
		for line in out.splitlines():
			if(re.search('Sampled',line)):
				tmp = line
				line = re.sub('^>\s+\w+\s+',"",line)
				line = re.sub('\s+.*$',"",line)
				if(int(line) == 0):
					flag = 0
				else:
					num = re.sub('^>\s+\w+\s+\d+\s+.*(?=\d)',"",tmp)
					num = re.sub('\s.*$',"",num)
					if(int(num) > MA.value):
						MA.value = int(num)
					flag = 1 
					continue
			if(flag == 1):
				if(re.search('\w+\s+',line)):
					line = re.sub('\-?\d+\.?\d+',"",line)
					line = re.sub('^\s+',"",line)
					line = re.sub('\\t',"",line)
					pair.append(line)
					pair.append(num)
					if pair not in fixedCG_d:
						fixedCG_d.append(pair)
					pair = []	
	#comptes RNAmutants
	def RNAmutants(self,seq,RNA_d,MA):
		out = commands.getoutput('./RNAmutants -l lib -m '+str(MUTATION_DEPTH)+' -s '+seq +' -n '+str(TOP))	
		for line in out.splitlines():
			if re.match('\w+',line):
				if re.match("C|G|c|G|A|a|U|u", line):
					line = re.sub('\s+.*$',"",line)
					mDepth = self.getMDepth(seq,line)	
					if(MA.value <  mDepth):
						MA.value = mDepth
					if line not in RNA_d:
						RNA_d.append(line)	
	def getMDepth(self,wild, seq):
		mDepth = 0
		for i in range(len(wild)):
			if(wild[i] != seq[i]):
				mDepth+=1
		return mDepth
###end of class collect_sequence
#get RANDOM variable for session
def setRandom():
	RANDOM = random.randint(1,100)
	return
def start_sampling(sequence):
	start = time()
	seq = collect_Sequence(sequence)
	smp = sampling(sequence, seq.RNA_entry, seq.fixedCG_entry)
	smp.compute()
	print "Definitions: R > RNAmutants, F > fixedCGSampling<br>"
	if(BOOTSTRAP==False):
		if(DEBUG):
			print smp.allSeq
		else:	
			for i in smp.allSeq:
				print str(i) +"<br>"
	if(BOOTSTRAP):
		boot = bootstrap(sequence, seq.MAX)
		fil = filterSeq(boot.bootstrap,smp.allSeq)
	end = time() - start
	print "Time(s): "+str(end)

def usage():
	print "-s [sequence]"
	print "-b BOOTSTRAP OPTION"
	print "-m [mutation depth]"
	print "-n [top x results]"
	print "-d [RNAfold dangling energyi]"
	print "-e [error]"
	print "-r ENABLES RNAmutant"
	print "-g ENABLES FixedCGSampling"
	print "-v Verbose option"
	print "sample"
	print "./script.py -b -r -g -n 10 -m 2 -s AGCGGGGGAGACAUAUAUCAUAGCCUGUCUCGUGCCCGACCCCGC"
	print "Wuff --- Wuff"
if __name__=='__main__':
	## collect input
	try: 
		optlist, args = getopt.getopt(sys.argv[1:], 'hbs:m:n:e:d:f:rgv', ["help"])
	except getopt.GetoptError, err:
		print str(err)
		print usage()
		sys.ext(2)
	seq = ''
	BOOTSTRAP = False 
	MUTATION_DEPTH = 5 
	TOP = 10
	DANGLING = "d0"
	ERROR = 0.05 
	FILE = ""
	RNAMUTANT = False
	FIXEDCGSAMPLING = False
	for opt, query in optlist:
		if opt == "-s":
			seq = query
		elif opt in ("-h", "--help"):
			usage()
			sys.exit
		elif opt == "-b":
			BOOTSTRAP = True
		elif opt == "-m":
			MUTATION_DEPTH = int(query)
		elif opt == "-n":
			TOP = query
		elif opt == "-e":
			ERROR = query
		elif opt == "-d":
			DANGLING = query
		elif opt == "-f":
			FILE = query
		elif opt == "-r":
			RNAMUTANT = True
		elif opt == "-g":
			FIXEDCGSAMPLING = True
		elif opt == "-v":
			DEBUG = True
		else:
			print "ERROR: incorrect usage"
			usage()
			sys.exit
	##start of program

	if(re.match('\w+',seq)):
		start_sampling(seq)
	else: 
		print "ERROR: Buggy sequence"
		sys.exit
##end of main
