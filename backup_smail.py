#! /usr/bin/env python
import smtplib
from email.MIMEText import MIMEText
from email.mime.multipart import MIMEMultipart

class mail:
	def __init__(self,address, wild,location):
		self.email = address
		self.wildtype = wild
		self.loc = location
	def sent(self, emsg):
		msg = MIMEMultipart("alternative")
		me = "no-reply@corrna.cs.mcgill.ca"
		msg['Subject'] = "Results " + str(self.wildtype)
		msg["From"]  = me
		msg["To"] = str(self.email)
		html = "<html>link: <a href='"+self.loc+"'>"+self.loc+"</a><br>"+emsg+"</html>"

		##part1 = MIMEText(text, 'plain')
		part2 = MIMEText(html, 'html')
		
		##msg.attach(part1)
		msg.attach(part2)
		
		try:
			s = smtplib.SMTP('localhost')
			s.sendmail(str(me), str(self.email), msg.as_string())
			s.quit
		except:
			print "Error: unable to send email"
			
