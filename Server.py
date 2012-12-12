#!/usr/bin/env python
from multiprocessing import Process, Queue
import socket
import sys
import commands
import script
import re
maxLOAD = 300
LOCATION = "/var/www/bio/results/"
#LOCATION = "../public_html/results/"
def process(client_socket, address, q):
	while 1:
		data = client_socket.recv(1024)
		if not data: break
		else: q.put(data)
def MP():
	while 1:
		while not q.empty():
			load = commands.getoutput("ps aux|awk 'NR > 0 { s +=$3 }; END {print s}'")
			if(float(maxLOAD) > float(load)):
				p = Process(target=execute, args=(q.get()))
				p.start()
def execute(*args):
	loc = ""
	param = ""
	flag = True
	for i in args:
		if(flag):
			if(re.match(" ",i)):
				flag = False
				next
			loc+=i
		else:
			param+=i
	loc = loc.strip()
	print "l: "+str(loc) + " p: " + str(param)
	print str(LOCATION+loc+"/index.php")
	f = open(LOCATION+loc+"/index.php",'w')
	f.write('<head> <META HTTP-EQUIV="Refresh" CONTENT="5"> </head>')
	f.write('Processing<br>')
	f.write('Page will auto refresh every 5 seconds<br>') 
	f.close()

	output = commands.getoutput(param)
	#output = re.sub("^.*[", "", output)
	out = output.split("\n")
	f = open(LOCATION+loc+"/index.php",'w')
	f.write('<!DOCTYPE HTML PUBLIC "-//W#C??DTD HTML 4.01 Transitional//EN">')
	f.write('<html xmlns="http://www.w3.oeg/1999/xhtml" xml:land="en"><head></head><body>')
	f.write('<script type="text/javascript">var output='+out[0]+';</script><div id="body" class="body"></div><br>'+out[1]+'<br>')
	f.write('<script src="http://code.jquery.com/jquery-1.4.3.min.js" type="text/javascript"></script><script src="../../app.js" type="text/javascript"></script>')
	f.write('<div style="display:none;"')
	f.write(output)
	f.write('</div>')
	f.write('</body></html>')
	f.close()
	print "done"

if __name__ == "__main__":
	q = Queue()
	server_socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
	server_socket.bind(("localhost", 5000))
	server_socket.listen(5)
	p = Process(target=MP, args=())
	p.start()

	print "TCP Socket Server Started"

	while 1:
		client_socket, address = server_socket.accept()
		p = Process(target=process, args=(client_socket, address,q))
		p.start()
		
			
	"""
	while 1:
		client_socket, address = server_socket.accept()
		print " i got a connection from ", address
		while 1:
			data = raw_input ( "SEND( TYPE q or Q to Quit):" )
			if (data == 'Q' or data == 'q'):
				client_socket.send (data)
				client_socket.close()
				break;
			else:
				client_socket.send(data)
	 
			data = client_socket.recv(512)
			if ( data == 'q' or data == 'Q'):
				client_socket.close()
				break;
			else:
				print "RECIEVED:" , data
	"""
