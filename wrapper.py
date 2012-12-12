#!/usr/bin/env python
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
from time import clock, time
from math import *
import Pearson
import smail
import BootStrap

DEBUG = False
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
OUTPUT_DATA = ""
EMAILMSG = ""

def usage():
	print "Displaying help:"
	print "-s [sequence]"
	print "-h display help file"
	print "-m [x] mutation depth"
	print "-n [x] number of samples"
	print "-l [x] loop size"
	print "-t [x] correlation threshold"
	print "-u [url location]"
	print "-d Debugging Mode"
	print "-c CSV *note*, requires debugging mode to be enabled"
	print "-a EMAIL"
	print "sample usage:\n ./wrapper.py -m 2 -n 10 -l 10 -t 0.4 -s AGCGGGGGAGACAUAUAUCAUAGCCUGUCUCGUGCCCGACCCCGC"

def SCOUT(seq):
#!/bin/bash

	if(DEBUG) :
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
	"""
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
	"""
	global AVG
	AVG = 0.814673348838
	global nAVG
	#nAVG = P.getCorr('AGCGGGGGAGACAUAUAUCACAGCCUGUCUAGGGCCCGACCCCGC');
	#modified to use threshold from command line
	if EMAIL:
		global EMAILMSG
		EMAILMSG = "<p>Under Avg: "+str(AVG)+" Threshold : "

		if(nAVG < 0):
			EMAILMSG+="none"
		else:
			EMAILMSG+=str(nAVG)
		EMAILMSG+=" </p>"
	if(DEBUG):
		if not CSV:
			print "<Under Avg> :: under "+str(AVG)+" across the 13 given confirmed sequence"
			print "<Under Threshold> :: under "+str(nAVG)+" on A84_U86G or user input threshold"
	else:
		global OUTPUT_DATA
		OUTPUT_DATA+='[{"method":"mutation","setAvg":"'+str(AVG)+'","setThres":"'+str(nAVG)+'","wildtype":"'+SEQUENCE+'","data":['

def getRNAfoldOptStruct(seq):
	rant = random.randint(1,100)
	rant_t = random.randint(1,100)
	fname = str(rant)+'-'+str(rant_t)
	f = open(fname+'.tmp','w')
	f.write('> '+fname+'\n')
	f.write(seq)
	f.close()
	out = commands.getoutput('cat '+fname+'.tmp @- | RNAfold -d0')		
	commands.getoutput('rm '+fname+'.tmp')
	commands.getoutput('rm '+fname+'*.ps')
	
	for line in out.splitlines():	
		if re.match("^[.|?]",line):
			line = re.sub("\s+.*$", "", line)
			return line

def initial():
	node = ''
	getAvg()
	global Question_Marks
	global NODE
	global totalTIME
	global totalCORRELATION
	global totalSEQUENCE
	global OUTPUT_DATA
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
	OUTPUT_DATA+='{"seq":['
	bool_first = 0
	underAvg = 0
	underThres = 0
	if EMAIL: 
		global EMAILMSG
		EMAILMSG+="<table border = 0><thead><td>Sequence</td><td>Structure</td><td>Correlation</td><td>MFE</td><td>Under Threshold</td></thead><tbody>"
	for line in out.splitlines():
		if re.match('\w+', line):
			if re.match("C|G|c|G|A|a|U|u", line):
				mfe = line
				mfe = re.sub('^\w+\s+\(',"",mfe)
				mfe = re.sub('\).*$',"",mfe)
				sequenceCOUNTER+=1
				line = re.sub('\s+.*$',"",line)
				struct = getRNAfoldOptStruct(line)
				x = P.getCorr(line)
				sequenceCORRELATION+=x
				note = ''
				if(AVG > x):
					underAvg = 1 
					note = " <Under Avg>"
				if(nAVG > x):
					underThres = 1
					SIGNIFICANT.append(line)
					note+=" <Under Threshold>"
				#print line +" "+str(x) + note
				if EMAIL:
					EMAILMSG+="<tr><td>"+str(line)+"</td><td>"+str(struct)+"</td><td>"+str(x)+"</td><td>"+str(mfe)+"</td><td>"	
					#if(AVG > x):	
					#	EMAILMSG+="Y"
					#EMAILMSG+="</td><td>"
					if(nAVG > x):
						EMAILMSG+="Y"
					EMAILMSG+="</td></tr>"
				if(DEBUG):
					if(CSV):
						print str(line) +","+struct+","+str(x)+","+str(mfe)
					else:
						print str(line) +","+struct+","+str(x)+","+str(mfe)+","+str(note)
					##print '{0}, {0}, {1:.4f}, {1:.4f}, {2}'.format(line,struct,x,mfe,note)
						SCOUT(line)
				else: 
					if(bool_first == 0):
						bool_first = 1
					elif bool_first == 1:
						OUTPUT_DATA+=','
					OUTPUT_DATA+='{"seq":"'+str(line)+'", "struct":"'+str(struct)+'", "corr":"'+str(x)+'","mfe":"'+str(mfe)+'","avg":"'
					if(underAvg == 1):
						OUTPUT_DATA+='t'
						underAvg = 0
					else:
						OUTPUT_DATA+='f'
					OUTPUT_DATA+='","navg":"'
					if underThres == 1:
						OUTPUT_DATA+='t'
						underThres = 0
					else:
						OUTPUT_DATA+='f'
					OUTPUT_DATA+='"}'
				node = locate(line, node)				
	totalCORRELATION+=sequenceCORRELATION
	totalSEQUENCE+=sequenceCOUNTER
	if EMAIL:
		EMAILMSG+="</tbody></table><p>Mutations: "+str(node)+" Avg: "+str(sequenceCORRELATION/sequenceCOUNTER)+" Time: "+str(end)+"</p>"
	if(DEBUG):
		if not CSV:
			print line
			print "avg correlation for this round "+str(sequenceCORRELATION/sequenceCOUNTER)
			print 'time took:   '+str(end)
	else:
		OUTPUT_DATA+='],"node":"'+str(node)+'","avg":"'+str(sequenceCORRELATION/sequenceCOUNTER)+'","time":"'+str(end)+'"}'
	loop(node)
def loop(node):	
	#print LOOP_SIZE
	global totalTIME
	global totalSEQUENCE
	global totalCORRELATION
	global OUTPUT_DATA
	if EMAIL:
		global EMAILMSG
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
		if(DEBUG):
			if not CSV:
				print node
		node = NODE
		start = time()
		out = commands.getoutput('./RNAmutants -l lib -m '+str(MUTATION_DEPTH) +' -n '+str(NSIZE) +' -f '+ wrapper +' -C')
		end = time() - start
		totalTIME+=end
		underAvg = 0
		underThres = 0
		bool_first = 0
		if EMAIL:
			EMAILMSG+="<table border = 0><thead><td>Sequence</td><td>Structure</td><td>Correlation</td><td>MFE</td><td>Under Threshold</td></thead><tbody>"
		OUTPUT_DATA += ',{"seq":['
	#	node = lockS(node)
		for line in out.splitlines():
			if re.match('\w+', line):
				if re.match("C|G|c|G|A|a|U|u", line):
					mfe = line
					mfe = re.sub('^\w+\s+\(',"",mfe)
					mfe = re.sub('\).*$',"",mfe)
					sequenceCOUNTER+=1
					line = re.sub('\s+.*$',"",line)
					struct = getRNAfoldOptStruct(line)
					x = P.getCorr(line)
					sequenceCORRELATION+=x
					note = ''
					if(AVG > x):
						note = " <Under Avg>"
						underAvg = 1
					if(nAVG > x):
						underThres = 1
						SIGNIFICANT.append(line)
						note+=" <Under Threshold>"
					#print line +" "+str(x) + note
					if EMAIL:
						EMAILMSG+="<tr><td>"+str(line)+"</td><td>"+str(struct)+"</td><td>"+str(x)+"</td><td>"+str(mfe)+"</td><td>"	
						#if(AVG > x):	
						#	EMAILMSG+="Y"
						#EMAILMSG+="</td><td>"
						if(nAVG > x):
							EMAILMSG+="Y"
						EMAILMSG+="</td></tr>"
					if(DEBUG):
						if(CSV):
							print str(line) +","+str(struct)+","+str(x)+","+str(mfe)
						else:
							print str(line) +","+str(struct)+","+str(x)+","+str(mfe)+","+str(note)
						##print '{0}, {0}, {1:.4f}, {1:.4f}, {2}'.format(line,struct,x,mfe,note)
							SCOUT(line)
					else:
						if(bool_first == 0):
							bool_first = 1
						elif bool_first == 1:
							OUTPUT_DATA+=','
						OUTPUT_DATA+='{"seq":"'+str(line)+'", "struct":"'+str(struct)+'","corr":"'+str(x)+'", "mfe":"'+str(mfe)+'","avg":"'
						if(underAvg == 1):
							OUTPUT_DATA+='t'
							underAvg = 0
						else:
							OUTPUT_DATA+='f'
						OUTPUT_DATA+='","navg":"'
						if underThres == 1:	
							OUTPUT_DATA+='t'
							underThres = 0
						else:		
							OUTPUT_DATA+='f'
						OUTPUT_DATA+='"}'
					node = locate(line, node)				
		if EMAIL:
			EMAILMSG+="</tbody></table><p>Mutations: "+str(node)+" Avg: "+str(sequenceCORRELATION/sequenceCOUNTER)+" Time: "+str(end)+"</p>"
		if(DEBUG):
			if not CSV:
				print "avg correlation for this round "+str(sequenceCORRELATION/sequenceCOUNTER)
				print 'time took:   '+str(end)
		else:
			OUTPUT_DATA+='],"node":"'+str(node)+'","avg":"'+str(sequenceCORRELATION/sequenceCOUNTER)+'","time":"'+str(end)+'"}'
	out = commands.getoutput('rm '+wrapper)
	if EMAIL:
		EMAILMSG += "<p> Avg Correlation: "+str(totalCORRELATION/totalSEQUENCE)+" Total Sequence: "+str(totalSEQUENCE)+"<br>Total Time: "+str(totalTIME)+" Total Avg Time: "+str(totalTIME/(LOOP_SIZE+1))+"</p>"	
		EMAILMSG += "Total Sequence Under Threshold: "+str(len(SIGNIFICANT))+"<br>"
		EMAILMSG += "<table border = 0><thead><td>Sequence</td><td>Correlation</td></thead><tbody>"
		for i in SIGNIFICANT:
			EMAILMSG +="<tr><td>"+str(i)+"</td><td>"+str(P.getCorr(i))+"</td></tr>"
		EMAILMSG+="</tbody></table>"
		global email
		email.sent(EMAILMSG)
		
	if(DEBUG):
		print "total avg correlation "+str(totalCORRELATION/totalSEQUENCE)
		print "total TIME took:: "+str(totalTIME) + " avg time took:: "+str(totalTIME/(LOOP_SIZE+1))
		print BUCKETS	
		print 'No of SIGS Under Threshold threshold {0} = {1}'.format(nAVG,len(SIGNIFICANT))
		print "# [sequence] [correlation]"
	#print SIGNIFICANT
		for i, seq in enumerate(SIGNIFICANT):
			print i, seq, "%.4f" % P.getCorr(seq) 
	else:
		OUTPUT_DATA+='], "tAvgCorr":"'+str(totalCORRELATION/totalSEQUENCE)+'", "tSeq":"'+str(totalSEQUENCE)+'", "tTime":"'+str(totalTIME)+'", "tAvgTime":"'+str(totalTIME/(LOOP_SIZE+1))+'","tThres":"'+str(len(SIGNIFICANT))+'"}]'
		print OUTPUT_DATA
	print "Time(s): "+str(totalTIME)
	
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
	EMAIL = False
	ADDRESS = ""
	CSV = False
	try:
		optlist, args = getopt.getopt(sys.argv[1:], 'cn:l:s:m:t:hda:u:', ["help"])
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
		elif opt =="-d":
			DEBUG = True
		elif opt == "-t":
			if(DEBUG):
				print 'Setting threshold to {0}'.format(query)
			nAVG = float(query)
		elif opt == "-u":
			location = query
		elif opt == "-a":
			EMAIL = True
			ADDRESS=query
		elif opt == "-c":
			CSV = True
		else:
			print "ERROR: incorrect usage"
			usage()
			sys.exit
	if EMAIL:
		email = smail.mail(ADDRESS, SEQUENCE,location)	
	if(re.match('\w+',SEQUENCE)):
		sampling = BootStrap.bootstrap(seq)
		initial()
	else:
		print "ERROR: Buggy or Null Sequence"
		sys.exit
