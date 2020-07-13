import sys
import os
import time
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--dir', type=str, default=None, help='Top-level project directory, otherwise will use current working directory.')
parser.add_argument('--env', type=str, default='env/lib', help='Path to your virtual env site-packages location.')

args = parser.parse_args()

#-- Create a system link in the users local python site-packages.

LINK_NAME = 'chp'

python_version_str = 'python{}.{}'.format(sys.version_info.major, sys.version_info.minor)
HOME_DIR = os.getenv("HOME")
if args.dir is None:
    WORKING_DIR = os.getcwd()
else:
    WORKING_DIR = args.dir
INSTALL_PATH = os.path.join(WORKING_DIR, args.env, python_version_str, 'site-packages')

LINK_PATH = os.path.join(INSTALL_PATH, LINK_NAME)


#-- Remove the system link if it exists
if os.path.exists(LINK_PATH):
    print('CHP system link already exists, so it is being removed and I am creating a new one.')
    os.system("rm {}".format(LINK_PATH))

if args.dir is None:
    os.system("ln -s {} {}".format(os.path.join(WORKING_DIR, 'chp'), LINK_PATH))
else:
    os.system("ln -s {} {}".format(os.path.join(WORKING_DIR, 'submodules/bkb-pathway-provider/chp'), LINK_PATH))

print('Now running PyBKB setup...')
#-- Run setup for PyBKB
os.system("cd submodules/PyBKB; python3 setup_for_development.py --env {} --no_copy".format(os.path.join(WORKING_DIR, args.env)))

print('Completed setup.')
