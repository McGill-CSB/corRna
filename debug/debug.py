#!/usr/bin/env python
import commands
import re

if __name__=="__main__":
	out = commands.getoutput("cat test.txt")	
	for line in out.splitlines():
		if(re.match("\[\{",line)):	
			print line

	
	
