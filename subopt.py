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
from time import clock, time
from math import *
import Pearson
import BootStrap
import smail

DEBUG = False
BOOTSTRAP = False
CSV = False
DUMMY = 0
start = time()
manager = Manager()
QUESTION_MARK = False
ALL = False
SEQUENCE = ""
MUTATION_DEPTH = 2
TOP_BREAK = 10
NSIZE = 10
EMAIL = False
JSON = False
location = ""

class subopt:
	def __init__(self):	
		return
def usage():
	print "./subopt.py -k -m 2 -n 10 -s AGCGGGGGAGACAUAUAUCAUAGCCUGUCUCGUGCCCGACCCCGC"
	print "-s [sequence]"
	print "-x question marks mode"
	print "-d debugging mode"
	print "-h help mode"
	print "-k prints the single brackets as well"
	print "-m [mutation depth]"
	print "-n [top seq]"
	print "-u [location]"
	print "-a [email]"
	print "-b [top N breakNumber]"
	print "-j json mode"
	print "-c CSV mode"
	print "-t bootstrap"
	print "Wuff --- Wuff"
def getbpp(fname):
	f = open(fname+'_dp.ps','r')
	ubox = []
	lbox = []
	pair = []
	index = -1
	for line in f:
		line.strip()
		if(re.search('\d+.*ubox',line)):
			line = re.sub('\\n','',line)
			pair = re.split('\s+',line)	
			ubox.append(pair)
		elif(re.search('\d+.*lbox',line)):
			line = re.sub('\\n','',line)
			pair = re.split('\s+',line)	
			lbox.append(pair)	
	for l in lbox:  ##optimize this search?
		for u in range(len(ubox)):
			if((l[0]==ubox[u][0])and(l[1]==ubox[u][1])):
				index = u	
		if(index >= 0 ):
			#print "pop"
		#	ubox.pop(u)
			index = -1
		else:
			print "ERROR: " + str(l) + " has no ubox match"
			sys.ext(2)
	out = sorted(ubox, compare)
	commands.getoutput('rm '+fname+'*.ps')
	return out
	#computes RNAfold
def compare(a,b):
	return cmp(b[1],a[1])
def RNAfold(seq):
	rant = random.randint(1,100)
	rant_t = random.randint(1,100)
	fname = str(rant)+'-'+str(rant_t)
	f = open(fname+'.tmp','w')
	f.write('> '+fname+'\n')
	f.write(seq)
	f.close()
	out = commands.getoutput('cat '+fname+'.tmp @- | RNAfold -p -d0')		
	commands.getoutput('rm '+fname+'.tmp')
	return fname
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
		if re.match("^[.|\(|?]",line):
			line = re.sub("\s+.*$", "", line)
			return line

def pickCandidates(bpp):
	G = sorted(bpp, compare)	
	visited = []
	optimal = []
	##initialize marker
	for vertex in G:
		if vertex not in visited:
		#	print ">>>"+str(vertex)
			visited.append(vertex)
			stack = []	
			stack.append(vertex)
			state = G.index(vertex)
			for v in range(state+1,len(G)):
				parent = stack.pop()
				try:
					x = int(G[v][0]) - int(parent[0])
					y = int(parent[1]) - int(G[v][1])	
					if(x == 1 and y ==1 ): #found a match
						stack.append(parent)
						stack.append(G[v])
						visited.append(G[v])
					elif( x < 0 and y < 0):
						break
					else:
						stack.append(parent)
				except IndexError as e:
					stack.append(parent)
			if(DEBUG):
				print stack
			if not ALL:
				if(len(stack)>1):
					optimal.append(stack)
			else:
				optimal.append(stack)
	opt = sorted(optimal, comp)
	return opt
def comp(a, b):
	return cmp(len(b),len(a))
def breakNumber(opt, struct):
##visualize break number
	Marker = ""
	number =""
	if DEBUG:
		for i in range(len(struct)):
			number+=" " +str(i)

		print number
		for i in range(len(struct)):
			print struct[i],
		print "	"

	L = list()
	T = list()
	R = list()
	W = list()
	buk = list()
	t_buk = list()
	

	for i in range(len(struct)):
		if(struct[i]=="("):	
			L.append(i)
		if(struct[i]==")"):
			##print str(L.pop()) + " " + str(i)
			T.append(L.pop())
			T.append(i)
			W.append(T)
			T = list()
	if(DEBUG):
		for i in W:
			print i


	for p in opt:
		v = []
		tmp = ""
		totalCount = 0
		leftCount = 0
		rightCount = 0
		detect = False
		detect2 = False
		node = ""
		for x in range(len(seq)):
				v.append(".")
		while p:
			n = p.pop()
			v[int(n[0])-1] = "("
			v[int(n[1])-1] = ")"
		for i in range(len(v)):
			tmp+=v[i]
			if(re.match('\.',v[i])):
				node+="?"
			else:
				node+=v[i]
			
		for i in range(len(tmp)):
			if(tmp[i]=="("):	
				L.append(i)
			if(tmp[i]==")"):
				##print str(L.pop()) + " " + str(i)
				T.append(L.pop())
				T.append(i)
				R.append(T)
				T = list()
		if(DEBUG):
			for i in R:
				print i
		breakNumber = 0
		for i in W:
			for j in R:
				if i[0] == j[0] and  i[1] == j[1]:
					continue
				elif i[0] == j[0] and i[1] != j[1]:
					breakNumber+=1
					break
				elif i[1] == j[1] and i[0] != j[0]:
					breakNumber+=1
					break
				elif i[0] == j[1] or j[0] == i[1]:
					breakNumber+=1
					break
				elif i[0]>j[0] and i[0]<j[1] and i[1] > j[1]:
					breakNumber+=1
					break
				elif i[1]>j[0] and i[1]<j[1] and i[0] < j[0]:
					breakNumber+=1 
					break
		if DEBUG:
			print struct
			print tmp + " b# " + str(breakNumber) ##+ " " + str(leftCount) + " " + str(rightCount)
		t_buk.append(node)
		t_buk.append(tmp)
		t_buk.append(breakNumber)
		buk.append(t_buk)
		t_buk = list()
		R = list()
	out = sorted(buk, compare)
	if DEBUG:
		print out
	return out
def compare(a, b):
	return cmp(float(b[2]),float(a[2]))

def RNAmutant(struct,breakSeq):
	if int(TOP_BREAK) > len(breakSeq):
		loop = len(breakSeq)
	else:
		loop = int(TOP_BREAK)
	Marker = ""
	proc = []
	bucket = manager.list()
	for i in range(len(SEQUENCE)):
		Marker+="o"
	for i in range(int(loop)):
		p = Process(target=_RNAmutant, args=(breakSeq[i][2],breakSeq[i][1],breakSeq[i][0],Marker,bucket))
		p.start()
		proc.append(p)
	for p in proc:
		p.join()
	return bucket

def printSeq(struct, boot):
	EMAILMSG = "<table><thead><td>Sequence</td><td>Structure</td><td>Correlation</td><td>MFE</td><td>BreakNumber Set</td><td>Structure Set</td><td>Significant</td></thead><tbody>"
	json = '[{"method":"structural","wild":"'+SEQUENCE+'","struct":"'+struct+'","seq":['
	for i in range(len(boot)):
		if(DEBUG):
			print boot[i],
			print "boot: "+str(sampling.getStanding(boot[i][0],boot[i][2])) 
		if JSON:
			json+= '{"seq":"'+str(boot[i][0])+'", "struct":"'+str(boot[i][1])+'", "corr":"'+str(boot[i][2])+'", "mfe":"'+str(boot[i][3])+'", "breakNumber":"'+str(boot[i][4])+'", "breakStruct":"'+str(boot[i][5])+'","boot":"'+str(sampling.getStanding(boot[i][0],boot[i][2]))+'"}'
			if i != len(boot)-1:
				json+=','
			
		if EMAIL:
			EMAILMSG += "<tr><td>"+str(boot[i][0])+"</td><td>"+str(boot[i][1])+"</td><td>"+str(boot[i][2])+"</td><td>"+str(boot[i][3])+"</td><td>"+str(boot[i][4])+"</td><td>"+str(boot[i][5])+"</td><td>"+str(sampling.getStanding(boot[i][0],boot[i][2]))+"</td></tr>"
		if CSV:
			print str(boot[i][0])+","+boot[i][1]+","+str(boot[i][2])+","+str(boot[i][3])+","+str(boot[i][4])+","+str(boot[i][5])+","+str(sampling.getStanding(boot[i][0],boot[i][2]))
	if JSON:
		json+=']}]'
		print json
	if EMAIL:
		EMAILMSG += "</tbody>"
		global email
		email.sent(EMAILMSG)


def _RNAmutant(bNumber,bStruct, breakStruct, Marker, bucket):
		tmp = list()
		rant = random.randint(1, 100);
		wrapper = str(rant)+'-subopt.in'	
		while(os.path.isfile(wrapper)):
			rant = random.randint(1, 100);
			wrapper = str(rant)+'-subopt.in'	
		f = open(wrapper, 'w')
		f.write(SEQUENCE+'\n')
		f.write(breakStruct+'\n')
		f.write(Marker+'\n')
		f.close()
		out = commands.getoutput('./RNAmutants -l lib -m '+str(MUTATION_DEPTH)+' -n '+str(NSIZE) +' -f '+wrapper+ ' -C')
		commands.getoutput('rm '+wrapper)	
		#begin parsing RNAmutant output
		for line in out.splitlines():
			if re.match("C|G|c|G|A|a|U|u", line):
				s = re.sub('\s+.*$',"",line)
				corr = P.getCorr(s)
				mfe = re.sub('^\w+\s+\(',"",line)
				mfe = re.sub('\).*$',"",mfe)
				struct = getRNAfoldOptStruct(line)
				if DEBUG:
					print str(s) + " " + str(struct)+" "+str(corr) + " " + str(mfe)
				tmp.append(s)
				tmp.append(struct)
				tmp.append(corr)
				tmp.append(mfe)
				tmp.append(bNumber)
				tmp.append(bStruct)
				bucket.append(tmp)
				tmp = list()
def echo(arr):
	for i in arr:
		print i
def start_sampling(seq):
	start = time()
	if BOOTSTRAP:
		if DEBUG:
			print "bootstrap started"
		P = Process(target=start_bootstrapAlgo, args=())
		P.start()
	fname = RNAfold(seq)
	bpp = getbpp(fname)
	if(DEBUG):
		echo(bpp)
	optimal = pickCandidates(bpp)
	struct = getRNAfoldOptStruct(seq)
	out = breakNumber(optimal, struct)  ##returns list with cmp b > a, returns struct1(?), struct2(.), breaknumber
	out = RNAmutant(struct, out)
	if BOOTSTRAP:
		P.join()
	printSeq(struct,out)
	end = time()-start
	print "Time(s) "+str(end)
def start_bootstrapAlgo():
	sampling.start()
	##visualize(optimal, seq)	
##get input
if __name__ == "__main__":
	try: 
		optlist, args = getopt.getopt(sys.argv[1:], 'khs:xda:m:n:b:ju:ct', ["help"])
	except getopt.GetoptError, err:
		print str(err)
		print usage()
		sys.ext(2)
	seq = ''
	for opt, query in optlist:
		if opt == "-s":
			seq = query
			SEQUENCE = query
			P = Pearson.Correlation(SEQUENCE)
		elif opt == "-x":
			QUESTION_MARK = True	
		elif opt == "-d":
			DEBUG = True
		elif opt == "-n":
			NSIZE = query
		elif opt == "-m":
			MUTATION_DEPTH = query
		elif opt == "-k":
			ALL = True
		elif opt in ("-h", "--help"):
			usage()
		elif opt == "-b":
			TOP_BREAK = query
		elif opt == "-j":
			JSON = True
		elif opt == "-a":
			EMAIL = True
			ADDRESS = query
		elif opt == "-u":
			location = query
		elif opt == "-c":
			CSV = True
		elif opt == "-t":
			BOOTSTRAP = True
		else:
			print "ERROR: incorrect usage"
			print usage()
			sys.exit
	if EMAIL:
		email = smail.mail(ADDRESS, SEQUENCE, location)
	if BOOTSTRAP:
		sampling = BootStrap.bootstrap(seq,MUTATION_DEPTH)
	if(re.match('\w',seq)):
		start_sampling(seq)
		
	else:
		print "ERROR: Buggy Sequence"
		sys.exit
