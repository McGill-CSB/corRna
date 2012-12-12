#!/usr/bin/env python
import smtplib
from email.mime.text import MIMEText

msg = "test"

s = smtplib.SMTP()
s.sendmail("kam.alfred@gmail.com","kam.alfred@gmail.com", msg)
