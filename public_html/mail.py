import smtplib
from email.mime.text import MIMEText
from email.mime.multipart import MIMEMultipart

class mail:
	def __init__(self,address, wild,location):
		self.email = address
		self.wildtype = wild
		self.loc = location
	def sent(self, emsg,wildtype):
		me = "no-reply@corrna.cs.mcgill.ca"
		msg['Subject'] = "Results " + str(self.wildtype)
		msg["From"]  = me
		msg["To"] = address
		html = "<html>link: <a href='"+self.loc+"'>"+self.loc+"</a><br>"+emsg+"</html>"
		msg.attach(MIMEText(html, 'html'))
		s = smtplib.SMTP('localhost')
		s.sendemail(me, address, msg.as_string())
		s.quit
			
