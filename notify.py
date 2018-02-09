#!/usr/bin/env python3

"""
usage: notify.py [-h] [-e EMAIL] [--add_email] [--view_config] [--ID ID]
                 [external commands [external commands ...]]

Automatically sends an email to the specified address upon completion of the
specified command. Useful primarily for very long-running processes. In many
cases, the command being run (including the name of the other program) will
need to be placed in quotes. If no email address is provided with the -e flag,
a prompt will be displayed based upon the configuration file.

positional arguments:
  external commands     External commands to run, including external program
                        call (default: None)

optional arguments:
  -h, --help            show this help message and exit
  -e EMAIL, --email EMAIL
                        the email address to notify (default: None)
  --add_email           add or change an email address in the config file
                        (default: False)
  --view_config         view the contents of the configuration file (default:
                        False)
  --ID ID               additional string to include in email subject
                        (default: None)
"""

import sys
import subprocess
import time
import os
import argparse
from collections import defaultdict
import smtplib
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText


def get_runtime(start_time, p=3):
    """
    Takes a start time and optional decimal precision p,
    returns a string of the total run-time until current
    time with appropriate units.

    """
    total = time.time() - start_time  # starts with seconds
    divided = total/60.0
    if divided < 2:
        run_time = total
        units = "seconds"
    elif divided < 60:
        run_time = divided
        units = "minutes"
    else:
        run_time = divided/60.0
        units = "hours"
    rounded = round(run_time, p)

    return "{} {}".format(rounded, units)


def names_from_config(config):
    """
    Read name/email pairs out of config file, and
    return a dictionary of names and emails.

    Config format for users is: user    name    email
    
    """
    user_info = defaultdict(dict)
    index = 1
    try:
        with open(config) as infile:
            for line in infile:
                if not line.startswith('user'):
                    continue
                _, name, email = line.strip().split('\t')
                user_info[name]['email'] = email
                user_info[name]['index'] = index
                index += 1
    except FileNotFoundError:
        pass

    return user_info


def get_config_info(config):
    """
    Checks config for necessary information and prompts for
    additional info as needed.
    
    """
    config_info = {}
    with open(config) as conf:
        for l in conf:
            if l.startswith('#'):
                continue
            if l.startswith('user'):
                continue
            try:
                key, value = l.strip().split()
                config_info[key] = value
            except ValueError:  # blank line
                continue
    info_prompts = {
        (0, 'server'): 'Server address for outgoing mail over SSL',
        (1, 'from_address'): 'Email address to send mail from',
        (2, 'password'): 'Server password for email (stored in plaintext)',
        (3, 'port'): 'Server port for outgoing mail over SSL (usually 587)'
    }
    missing_info = {
        k: v for k, v in info_prompts.items() 
        if k[1] not in config_info.keys()}
    if missing_info:
        print('Please provide the following server config info: ')
        for k, v in sorted(missing_info.items()):
            target = k[1]
            config_info[target] = input('{}: '.format(v))
        should_write = input(
            'Write the provided information to config (y/n): '
        )
        if should_write.lower() == 'y':
            write_config_info(config, config_info)
            print('Information written to \'{}\''.format(config))
    
    return config_info


def write_config_info(config, info):
    """
    Updates info in config with dict >info<, maintaining 
    existing information unless overridden by >info<.

    """
    # keep any existing user info separate
    config_info = {}
    user_info = []
    with open(config) as conf:
        for l in conf:
            if l.startswith('#'):
                continue
            l = l.strip()
            if l.startswith('user'):
                user_info.append(l)
            else:
                try:
                    key, value = l.split()
                    config_info[key] = value
                except ValueError:  # blank line
                    continue
    # add new config info to existing, replacing
    # as needed
    config_info.update(info)
    with open(config, 'w') as conf:
        conf.write('#' * 80)
        conf.write('\n')
        for k, v in sorted(config_info.items()):
            conf.write('\t'.join([k, v]) + '\n')
        conf.write('#' * 80)
        conf.write('\n\n')
        conf.write('\n'.join(user_info) + '\n')


def email_from_config(config):
    """
    Retrieves user/email pairs from the config file and
    prompts the user to pick. Allows user to add new name to config.

    Returns an email address.
    
    """
    user_dict = names_from_config(config)
    if not user_dict:
        print('No users/emails found in config. Please add one now.')
        name, email = get_user_info()
        add_replace_email(name, email, config)
        return email
    selection_dict = {}
    for name, info in sorted(user_dict.items(), key=lambda v: v[1]['index']):
        index = str(info['index'])
        email = info['email']
        selection_dict[index] = name
        print('{}. {}'.format(index, name))
    user_choice = input('Select user number (\'a\' to add, \'q\' to quit): ')
    if user_choice == 'a':
        new_name, new_email = get_user_info()
        add_replace_email(new_name, new_email, config)
        return new_email
    elif user_choice == 'q':
        sys.exit('Exiting at user request')
    user_name = selection_dict[user_choice]

    return user_dict[user_name]['email']


def add_replace_email(name, email, config):
    """
    Replaces the >email< for >name< in >config< in place.
    
    """
    new_config = '{}.temp'.format(config)
    replaced = False
    new_entry = '\t'.join(['user', name, email])
    with open(config) as oldfile, open(new_config, 'w') as newfile:
        for line in oldfile:
            line = line.strip()
            if name in line:
                line = new_entry
                replaced = True
            newfile.write(line + '\n')
        if not replaced:  # is a new name entry
            newfile.write(new_entry + '\n')
    os.rename(new_config, config)
    print('\'{}\' added to \'{}\'.'.format(name, config))


def get_user_info():
    user_name = input('Name to add to config: ')
    user_email = input('Email address for user \'{}\': '.format(user_name))

    return (user_name, user_email) 
                

def add_name(name, email, config):
    """
    Adds a new name to the config file.
    
    """
    current_names = names_from_config(config)
    if name in current_names:
        current_email = current_names[name]['email']
        new_email = None
        replace_name = None
        while replace_name not in ('y', 'n'):
            replace_name = input(
                'Name already found in config file with email address '
                '\'{}\'\nReplace email (y/n): '.format(current_email))
            if replace_name.lower() == 'y':
                # new_email = input('Enter new email address: ')
                add_replace_email(name, email, config)
            elif replace_name.lower() == 'n':
                print('Leaving email for user \'{}\' intact'.format(name))
            else:
                print('Input not understood (must be \'y\' or \'n\')')
    else:
        add_replace_email(name, email, config)


def view_config(config):
    with open(config) as f:
        print(f.read(), end='')


def send_ssl_mail(
        from_address,
        to_address,
        server_address,
        port,
        password,
        subject=None,
        body=None):
    """
    Sends an email via notify@roylab.science, with optional subject-line
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
    description='Automatically sends an email to the specified address upon '
    'completion of the specified command. Useful primarily for very long-'
    'running processes. In many cases, the command being run (including '
    'the name of the other program) will need to be placed in quotes. '
    'If no email address is provided with the -e flag, a prompt will be '
    'displayed based upon the configuration file.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument(
    'commands',
    metavar='external commands',
    nargs='*',
    help='External commands to run, including external program call'
)
parser.add_argument(
    '-e',
    '--email',
    help='the email address to notify'
)
parser.add_argument(
    '--add_email',
    help='add or change an email address in the config file',
    action='store_true'
)
parser.add_argument(
    '--view_config',
    help='view the contents of the configuration file',
    action='store_true'
)
parser.add_argument(
    '--ID',
    help='additional string to include in email subject',
    type=str
)

if len(sys.argv) == 1:
    sys.exit(parser.print_help())

SCRIPT_HOME = os.path.dirname(os.path.realpath(sys.argv[0]))
HOME_DIR = os.path.expanduser('~')
CONFIG = os.path.join(HOME_DIR, '.notify.py.config')
if not os.path.isfile(CONFIG):
    open(CONFIG, 'w').close()  # hacky!

args = parser.parse_args()

if args.view_config:
    view_config(CONFIG)
    sys.exit(0)

if args.add_email:
    name, TARGET_EMAIL = get_user_info()
    add_name(name, TARGET_EMAIL, CONFIG)
    sys.exit(0)

if args.email:
    TARGET_EMAIL = args.email
else:  # use config
    TARGET_EMAIL = email_from_config(CONFIG)

# Make sure we have a correct email address to send to
if "@" not in TARGET_EMAIL:
    sys.exit("Email address missing '@' symbol. Exiting.")

# check other config information
config_info = get_config_info(CONFIG)

CMDS = args.commands

# assume first argument is external program name
REF_NAME = CMDS[0]

CMD_STRING = ' '.join(CMDS)

# Run external script
start_time = time.time()
run_dir = os.getcwd()
return_code = subprocess.call(CMD_STRING, shell=True)
run_time = get_runtime(start_time)

# For list-style format of email
result = return_code

# Get machine hostname and time
host = subprocess.check_output(["hostname"])
tstring = "%Y.%m.%d-%H.%M"
sys_time = '[{}]'.format(time.strftime(tstring))

# format the email subject line depending on provided info 
if not host:
    host = ""
else:
    host = "[{}]".format(host.decode("utf-8").strip())

if args.ID:
    id_string = ' {}'.format(args.ID)
else:
    id_string = ''

host_prefix = '{}{}{}'.format(sys_time, host, id_string)

# Message subject line
msg_subject = "{}: '{}' completed".format(host_prefix, REF_NAME)

# Completion message
msg_body = "\n".join([
    "Command-line arguments: {}".format(CMD_STRING),
    "Total runtime: {}".format(run_time),
    "Return value: {}".format(result),
    "Location: {}".format(run_dir)
    ])

send_args = {
    'from_address': config_info['from_address'],
    'to_address': TARGET_EMAIL,
    'server_address': config_info['server'],
    'port': config_info['port'],
    'password': config_info['password'],
    'subject': msg_subject,
    'body': msg_body}

send_ssl_mail(**send_args)

print('[#] Command completed in {}'.format(run_time), file=sys.stderr)

sys.exit(0)
