# script to send email when job finishes

# before this need to set up a gmail app password by going to account security page
#   enable 2 step verification
#   then under "how you sign in to google" --> click "app passwords" & generate a new password (i named it my app "BioHPC")
#   it will generate a 16 character password that is needed below 

import smtplib
from email.mime.text import MIMEText

#! EDIT EMAIL AND PASSWORD
sender_email = "your_email@gmail.com"
receiver_email = "your_email@gmail.com"
app_password = "your_16_char_app_password"
#! EDIT EMAIL AND PASSWORD ^

subject = "ur script is complete :)"
body = "howdy partner ur script is complete.\n\nFinished at: {}".format(__import__('datetime').datetime.now())

# compose message
msg = MIMEText(body)
msg['Subject'] = subject
msg['From'] = sender_email
msg['To'] = receiver_email

# send
try:
    with smtplib.SMTP_SSL("smtp.gmail.com", 465) as server:
        server.login(sender_email, app_password)
        server.sendmail(sender_email, receiver_email, msg.as_string())
    print("email sent successfully!")
except Exception as e:
    print("failed to send email:", e)
