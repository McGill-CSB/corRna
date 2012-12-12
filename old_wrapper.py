#!/usr/bin/env python
import sys
import re
import commands
import shutil
import getopt
import os, glob, sys
import os.path
import random, math
from decimal import Decimal
from time import clock, time
from math import *
import Pearson

DEBUG = True
MUTATION_DEPTH = 2
NSIZE = 10
SEQUENCE = ''
LOOP_SIZE = 10
BUCKETS = []
NODE = ''
Question_Marks = ''
AVG = float(0)
nAVG = float(0)
totalCORRELATION = 0
totalSEQUENCE = 0
totalTIME = 0
SIGNIFICANT = []

def usage():
	print "Displaying help:"
	print "-s [sequence]"
	print "-h display help file"
	print "-m [x] mutation depth"
	print "-n [x] number of samples"
	print "-l [x] loop size"
	print "-t [x] correlation threshold"
	print "sample usage:\n ./wrapper.py -m 2 -n 10 -l 10 -t 0.4 -s AGCGGGGGAGACAUAUAUCAUAGCCUGUCUCGUGCCCGACCCCGC"
def SCOUT(seq):
#!/bin/bash
	if(re.match('AGCGGGGGAGACAUCUAUCACAGCCUGUCUCGUGCCCGACCCCGC', seq)):
		print "A68C"
	if re.match('AGCGGGGGAGACAUAUACCACAGCCUGUCUCGUGCCCGACCCCGC', seq):
		print "U71C"
	if re.match('AGCGGGGGAGACAUAUAUCAUAGCCUGUCUCGUGCCCGACCCCGC', seq):
		print "C74U"
	if re.match('AGCGGGGGAGACAUAUAUCACUCUCUGUCUCGUGCCCGACCCCGC', seq):
		print "A75U_G76C_C77U"
	if re.match('AGCGGGGGAGACAUAUAUCACAGUCUGUCUCGUGCCCGACCCCGC', seq):
		print "C77U"
	if re.match('AGCGGGGGAGACAUAUAUCACAGCCUGUCUAGGGCCCGACCCCGC', seq):
		print  "C84A_U86G"
	if re.match('AGCGGGGGAGACAUAUAUCACAGCCUGUCUAGAGCCCGACCCCGC', seq):
		print  "C84A_U86A"
	if re.match('AGCGGGGGAGACAUAUAUCACAGCCUGUCUCGGGCCCGACCCCGC', seq):
		print "U86G"
	if re.match('AGCGGGGGAGACAUAUAUCACAGCCUGUCUCGUGCCAGGCCCCGC', seq):
		print "C90A_A92G"
	if re.match('AGCGGGGGAGACAUAUAUCACAGCCUGUCUCGUGCCCGGCCCCGC', seq):
		print "C92G"
	if re.match('AGCGGGGGAGACAUAUAUCACAGCCUGUCUCGUGCCAGACCCCGC', seq):
		print  "C90A"
	if re.match('AGCGGGGGAGACAUAUAUCACAGCCUGUCUCGUGCCCGUCCCCGC', seq):
		print "A92U"
	if re.match('AGCGGGGGAGACAUAUAUCACAGCCUGUCUCGUGCCCGCCCCCGC', seq):
		print "A92C"
def getAvg():
	x = float(0.0)
	x += P.getCorr('AGCGGGGGAGACAUCUAUCACAGCCUGUCUCGUGCCCGACCCCGC');
	x += P.getCorr('AGCGGGGGAGACAUAUACCACAGCCUGUCUCGUGCCCGACCCCGC');
	x += P.getCorr('AGCGGGGGAGACAUAUAUCAUAGCCUGUCUCGUGCCCGACCCCGC');
	x += P.getCorr('AGCGGGGGAGACAUAUAUCACUCUCUGUCUCGUGCCCGACCCCGC');
	x += P.getCorr('AGCGGGGGAGACAUAUAUCACAGUCUGUCUCGUGCCCGACCCCGC');
	x += P.getCorr('AGCGGGGGAGACAUAUAUCACAGCCUGUCUAGGGCCCGACCCCGC');	
	x += P.getCorr('AGCGGGGGAGACAUAUAUCACAGCCUGUCUAGAGCCCGACCCCGC');
	x += P.getCorr('AGCGGGGGAGACAUAUAUCACAGCCUGUCUCGGGCCCGACCCCGC');
	x += P.getCorr('AGCGGGGGAGACAUAUAUCACAGCCUGUCUCGUGCCAGGCCCCGC');
	x += P.getCorr('AGCGGGGGAGACAUAUAUCACAGCCUGUCUCGUGCCCGGCCCCGC');
	x += P.getCorr('AGCGGGGGAGACAUAUAUCACAGCCUGUCUCGUGCCAGACCCCGC');
	x += P.getCorr('AGCGGGGGAGACAUAUAUCACAGCCUGUCUCGUGCCCGUCCCCGC');
	x += P.getCorr('AGCGGGGGAGACAUAUAUCACAGCCUGUCUCGUGCCCGCCCCCGC');
	y = float(x)/float(13) 
	global AVG
	AVG = y
	global nAVG
	#nAVG = P.getCorr('AGCGGGGGAGACAUAUAUCACAGCCUGUCUAGGGCCCGACCCCGC');
	#modified to use threshold from command line
	print "<Under Avg> :: under "+str(AVG)+" across the 13 given confirmed sequence"
	print "<Under Threshold> :: under "+str(nAVG)+" on A84_U86G or user input threshold"

def initial():
	node = ''
	getAvg()
	global Question_Marks
	global NODE
	global totalTIME
	global totalCORRELATION
	global totalSEQUENCE
	sequenceCOUNTER = 0
	sequenceCORRELATION = 0
	for i in range(len(SEQUENCE)):
		BUCKETS.append(int(0))
		Question_Marks+='?'
		node+="o"
	#Question_Marks = '?((((((????????????????????????????????))))))'
	NODE = node
	start = time()
	out = commands.getoutput('./RNAmutants -l lib -m '+str(MUTATION_DEPTH)+' -s '+ SEQUENCE +'-n '+str(NSIZE))
	end = time() - start
	totalTIME+=end
	for line in out.splitlines():
		if re.match('\w+', line):
			if re.match("C|G|c|G|A|a|U|u", line):
				sequenceCOUNTER+=1
				line = re.sub('\s+.*$',"",line)
				x = P.getCorr(line)
				sequenceCORRELATION+=x
				note = ''
				if(AVG > x):
					note = " <Under Avg>"
				if(nAVG > x):
					SIGNIFICANT.append(line)
					note+=" <Under Threshold>"
				#print line +" "+str(x) + note
				print '{0} {1:.4f} {2}'.format(line,x,note)
				SCOUT(line)
				node = locate(line, node)				
	totalCORRELATION+=sequenceCORRELATION
	totalSEQUENCE+=sequenceCOUNTER
	print "avg correlation for this round "+str(sequenceCORRELATION/sequenceCOUNTER)
	print 'time took:   '+str(end)
	loop(node)
def loop(node):	
	#print LOOP_SIZE
	global totalTIME
	global totalSEQUENCE
	global totalCORRELATION
	#lock file name for this session
	rant = random.randint(1, 100);
	wrapper = str(rant)+'-wrapper.in'	
	while(os.path.isfile(wrapper)):
		rant = random.randint(1, 100);
		wrapper = str(rant)+'-wrapper.in'	
	for i in range(int(LOOP_SIZE)):
		sequenceCOUNTER =0
		sequenceCORRELATION = 0
		f = open(wrapper, 'w')
		f.write(SEQUENCE+'\n')
		f.write(node+'\n')
		f.write(Question_Marks+'\n')
		f.close()
		print node
		node = NODE
		start = time()
		out = commands.getoutput('./RNAmutants -l lib -m '+str(MUTATION_DEPTH) +' -n '+str(NSIZE) +' -f '+ wrapper +' -C')
		end = time() - start
		totalTIME+=end
	#	node = lockS(node)
		for line in out.splitlines():
			if re.match('\w+', line):
				if re.match("C|G|c|G|A|a|U|u", line):
					sequenceCOUNTER+=1
					line = re.sub('\s+.*$',"",line)
					x = P.getCorr(line)
					sequenceCORRELATION+=x
					note = ''
					if(AVG > x):
						note = " <Under Avg>"
					if(nAVG > x):
						SIGNIFICANT.append(line)
						note+=" <Under Threshold>"
					#print line +" "+str(x) + note
					print '{0} {1:.4f} {2}'.format(line,x,note)
					SCOUT(line)
					node = locate(line, node)				
		print "avg correlation for this round "+str(sequenceCORRELATION/sequenceCOUNTER)
		print 'time took:   '+str(end)
	out = commands.getoutput('rm '+wrapper)
	print "total avg correlation "+str(totalCORRELATION/totalSEQUENCE)
	print "total TIME took:: "+str(totalTIME) + " avg time took:: "+str(totalTIME/(LOOP_SIZE+1))
	print BUCKETS	
	print 'No of SIGS Under Threshold threshold {0} = {1}'.format(nAVG,len(SIGNIFICANT))
	print "# [sequence] [correlation]"
	#print SIGNIFICANT
	for i, seq in enumerate(SIGNIFICANT):
		print i, seq, "%.4f" % P.getCorr(seq) 
def lockS(node):
	for i in range(len(BUCKETS)):
		if(BUCKETS[i] != 0 and BUCKETS[i] % 10 == 0):
			node = node[:i]+"X"+node[i+1:]
			BUCKETS[i] = 0
	return node
def locate(seq, node):
	i=0
	for x in seq:
		if re.match("c|g|a|u", x):
			if not re.match("X", node[i]):
				node = node[:i]+"X"+node[i+1:]
				BUCKETS[i]+=1;
		i+=1;
	return node
			
if __name__ == "__main__":
	try:
		optlist, args = getopt.getopt(sys.argv[1:], 'n:l:s:m:t:h', ["help"])
	except getopt.GetoptError, err:
		print str(err)
		print usage()
		sys.exit(2)
	for opt, query in optlist:
		if opt == "-s":
			SEQUENCE = query
			P = Pearson.Correlation(SEQUENCE)
		elif opt in ("-h", "--help"):
			usage()
			sys.exit(2)
		elif opt == "-m":
			MUTATION_DEPTH = int(query)
		elif opt == "-n":
			NSIZE = int(query)
		elif opt == "-l":
			LOOP_SIZE = int(query)
		elif opt == "-t":
			print 'Setting threshold to {0}'.format(query)
			nAVG = float(query)
			
		else:
			print "ERROR: incorrect usage"
			usage()
			sys.exit
	if(re.match('\w+',SEQUENCE)):
			initial()
	else:
		print "ERROR: Buggy or Null Sequence"
		sys.exit
