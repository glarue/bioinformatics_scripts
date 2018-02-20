#!/usr/bin/env python3

"""
usage: cli_mail.py [-h] -sadd SERVER_ADDRESS -fadd FROM_ADDRESS -pw PASSWORD
                   -tadd TO_ADDRESS [-p PORT] [-s [SUBJECT]] [-b [BODY]]
                   [--ID ID] [--suffix_tag]

Send an email over SSL from the command line

optional arguments:
  -h, --help            show this help message and exit
  -p PORT, --port PORT  Port to use on email server (default: 587)
  -s [SUBJECT], --subject [SUBJECT]
                        A subject line for the email, wrapped in single quotes
                        ('') (default: None)
  -b [BODY], --body [BODY]
                        The body of the email, wrapped in single quotes ('')
                        (default: None)
  --ID ID               Arbitrary string to include in email header (default:
                        )
  --suffix_tag          place time/server info tag at end of subject rather
                        than beginning (default: False)

required named arguments:
  -sadd SERVER_ADDRESS, --server_address SERVER_ADDRESS
                        Server address for outgoing email using the SSL
                        protocol (default: None)
  -fadd FROM_ADDRESS, --from_address FROM_ADDRESS
                        Full email address from which to send email (default:
                        None)
  -pw PASSWORD, --password PASSWORD
                        Account password for the 'from' address. Wrap in
                        single quotes ('') (default: None)
  -tadd TO_ADDRESS, --to_address TO_ADDRESS
                        Destination email address (default: None)

"""

import sys
import smtplib
import subprocess
import argparse
import time
from time import strftime
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText


def send_ssl_mail(
        from_address,
        to_address,
        server_address,
        port,
        password,
        subject=None,
        body=None):
    """
    Sends an email via specified server, with optional subject-line
    and body arguments.

    Adapted from http://naelshiab.com/tutorial-send-email-python/

    """
    msg = MIMEMultipart()
    msg['From'] = from_address
    msg['To'] = to_address
    if subject:
        msg['Subject'] = subject
    if body:
        msg.attach(MIMEText(body, 'plain'))
    try:
        server = smtplib.SMTP_SSL(server_address, port)
    except smtplib.SMTPConnectError:
        print('[#] Server connection error - retrying', file=sys.stderr)
        time.sleep(10)
        server = smtplib.SMTP_SSL(server_address, port)
    retries = 2
    success = False
    while retries > 0:  # in case server rejects attempt
        # server.starttls()  # not used for SMTP_SSL class
        try:
            server.login(from_address, password)
            success = True
        except:
            server.quit()
            time.sleep(30)  # sleep for 30 seconds
            retries -= 1
            continue
        break
    if not success:
        sys.exit("{} error: connection to server could not be established".
                 format(sys.argv[0]))
    text = msg.as_string()
    server.sendmail(from_address, to_address, text)
    server.quit()


parser = argparse.ArgumentParser(
    description="Send an email over SSL from the command line",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# create an argument group for all of the required "optional" arguments,
# because Python by default labels all arguments with names as optional
# in the usage info
required_named = parser.add_argument_group('required named arguments')
required_named.add_argument(
    "-sadd",
    '--server_address',
    help='Server address for outgoing email using the SSL protocol',
    required=True
)
required_named.add_argument(
    "-fadd",
    "--from_address", 
    help="Full email address from which to send email",
    required=True)
required_named.add_argument(
    "-pw",
    "--password",
    help="Account password for the 'from' address. "
    "Wrap in single quotes ('')",
    required=True)
required_named.add_argument(
    "-tadd",
    "--to_address", 
    help="Destination email address",
    required=True)
parser.add_argument(
    "-p",
    "--port",
    help='Port to use on email server',
    type=int,
    default=587
)
parser.add_argument(
    "-s",
    "--subject",
    nargs="?",
    help="A subject line for the email, "
    "wrapped in single quotes ('')")
parser.add_argument(
    "-b",
    "--body",
    nargs="?",
    help="The body of the email, wrapped in "
    "single quotes ('')")
parser.add_argument(
    '--ID',
    type=str,
    help='Arbitrary string to include in email header',
    default=''
)
parser.add_argument(
    "--suffix_tag",
    action="store_true",
    help="place time/server info tag at end of subject rather than beginning"
)

if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)

# Get machine hostname and time
host = subprocess.check_output(["hostname"])

args = parser.parse_args()

# get system time as string
tstring = "%m-%d-%y@%H:%M"
sys_time = '[{}]'.format(time.strftime(tstring))

if not host:
    host = ""
else:
    host = "[{}]".format(host.decode("utf-8").strip())

host_prefix = '{}{}{}'.format(sys_time, host, args.ID)

if args.subject:
    subject_bits = [host_prefix, args.subject]
    if args.suffix_tag:
        subject_bits.reverse()
    args.subject = "{} {}".format(*subject_bits)
else:
    args.subject = host_prefix

arg_dict = vars(args)  # make a dictionary of args to use in function call
del arg_dict['ID']
del arg_dict['suffix_tag']
send_ssl_mail(**arg_dict)  # unpack dictionary to kwargs

sys.exit(0)
