#! /usr/bin/env python
#import numeric
from math import *
import sys
import re
import getopt

### note ### we are assumming the squence length of typeA and typeB are the same!

class basePairProbability:
	#constructor
	def __init__(self, *args):
		self.filename = ""
		self.entries = []
		self.bpp = []
		self.sequence = ''
		self.sequenceSize = 0	
		self.sumColumn = []
		if len(args) != 0:
			self.filename = args[0]	

	def parse(self):
		if self.filename != "":
			infile = open(self.filename,"r")
			line = infile.readline()		
			#	line = sys.stdin.readline()
			self.sequence = line
			self.sequenceSize = len(line) -1
			#print self.sequenceSize
			for i in range(self.sequenceSize):
				self.bpp.append(float(0))
			while infile:
				#for line in sys.stdin.readlines():
				line = infile.readline()
				line.strip()
				data = re.split('\s+', line)
				entry = []
	#			print data
				if data[0] == '':
					break 
				entry.append(int(data[0]))
				entry.append(int(data[1]))
				entry.append(float(data[2]))
				self.entries.append(entry)
				self.bpp[entry[0]-1]+=entry[2]
				self.bpp[entry[1]-1]+=entry[2]
				#delete the stuff below when you see it
				#self.bpp[self.sequenceSize-entry[0]-1]+=entry[2]
			print "Base Pair Probability"
			print  self.bpp	

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

		print "sum_XX =", sum_XX
        	print "sum_YY =", sum_YY
                print "sum_XY =", sum_XY

		self.r = sum_XY / (sqrt(sum_XX*sum_YY))  ##pearson's coefficient correlation

def usage():
	print "-x [file one]"
	print "-y [file two]"
## end of class 

#line = sys.stdin.readline()
#args = line.split()
try: 
	optlist, args = getopt.getopt(sys.argv[1:], 'x:y:h', ["help"])
except getopt.GetoptError, err:
	print str(err)
	sys.ext(2)
for opt, query in optlist:
	if opt == "-x":
		typeA = basePairProbability(query)
		typeA.parse()
	elif opt == "-y":
		typeB = basePairProbability(query)
		typeB.parse()
	elif opt in ("-h", "--help"):
		usage()
		sys.exit()	
	else:
		assert False, "unhandled option"
		 	
correlation = pearson(typeA, typeB)
correlation.compute()
#print "The Pearson's coefficient of correlation is %s against %s " %(correlation.r,typeB.filename)
print "wildtype size: " + str(typeA.sequenceSize)
print "the other size: " + str(typeB.sequenceSize)
print str(typeB.filename)
print "%s" %(correlation.r)
"""
f = open('bpp-results.txt','a')
#f.write(str(typeA.filename)+ '\n')
#for i in typeA.bpp:
#	f.write(str(i)+'\n')
f.write(str(typeB.filename)+ '\n')
for i in typeB.bpp:
	f.write(str(i)+'\n')
f.write('correlation:'+str(correlation.r)+'\n')
f.close()
"""	
