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

DEBUG = False
DUMMY = 0
start = time()
manager = Manager()
QUESTION_MARK = False

class subopt:
	def __init__(self):	
		return
def usage():
	print "./subopt.py -s AGCGGGGGAGACAUAUAUCAUAGCCUGUCUCGUGCCCGACCCCGC"
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
			#"""
			if(len(stack)>1):
				optimal.append(stack)
			#"""
			#optimal.append(stack)
	opt = sorted(optimal, comp)
	return opt
def comp(a, b):
	return cmp(len(b),len(a))
def visualize(opt, seq):
	for x in seq:
		print x,
	print ""
	for p in opt:
		v = []
		for x in range(len(seq)):
			if(QUESTION_MARK):
				v.append("?")
			else:
				v.append(".")
		while p:
			n = p.pop()
			v[int(n[0])-1] = "("
			v[int(n[1])-1] = ")"
		for i in range(len(v)):
			print v[i],
		print ""
def echo(arr):
	for i in arr:
		print i
def start_sampling(seq):
	fname = RNAfold(seq)
	bpp = getbpp(fname)
	if(DEBUG):
		echo(bpp)
	optimal = pickCandidates(bpp)
	visualize(optimal, seq)	
##get input
if __name__ == "__main__":
	try: 
		optlist, args = getopt.getopt(sys.argv[1:], 'hs:xd', ["help"])
	except getopt.GetoptError, err:
		print str(err)
		print usage()
		sys.ext(2)
	seq = ''
	for opt, query in optlist:
		if opt == "-s":
			seq = query
		elif opt == "-x":
			QUESTION_MARK = True	
		elif opt == "-d":
			DEBUG = True
		elif opt in ("-h", "--help"):
			usage()
		else:
			print "ERROR: incorrect usage";
			usage()
			sys.exit
	if(re.match('\w',seq)):
		start_sampling(seq)
	else:
		print "ERROR: Buggy Sequence"
		sys.exit
