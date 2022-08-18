'''
This script takes files from the current directory and submits them to
the online Tapenade interface for automatic differentiation.
Then this script copies the downloaded files back to the current directory and
extract them to a folder named TapenadeResults_d or TapenadeResults_b for
forward or reverse differentiation, respectively.

Usage:
1- Check if the variable within the 'USER INPUTS' field are correct

2- Run the code with one of the following options:

   For forward AD:
   $ python3 send_to_tapenade.py f

   For reverse AD:
   $ python3 send_to_tapenade.py r

Ney Rafael Secco
2019-03-05
'''

import os
import shutil
import selenium
from selenium import webdriver
import time
import argparse
import fileinput

# GET DIFF MODE FROM COMMAND LINE ARGUMENT

# Create parser
parser = argparse.ArgumentParser(description='Send files to the Tapenade web interface for automatic differentiation.')
parser.add_argument('mode', action='store',
                    help='differentiation mode: f (forward)) or r (reverse)')
args = parser.parse_args()

if args.mode == 'f':
    mode = 'forward'
elif args.mode == 'r':
    mode = 'reverse'
else:
    raise ValueError('Differentiation mode not recognized. Use either f or r.')

# USER INPUTS:

# Default download directory of the browser
download_directory = os.path.expanduser('~/Downloads')

# Standard filename of the Tapenade Download File
download_filename = 'TapenadeResults.zip'

# List of files that will be differentiated
files = ['fem_module.f90', 'llt_module.f90', 'asa_module.f90']

# Name of top routine
top_routine = 'asa_analysis'

# Output variables (separate with whitespace)
output_vars = 'resllt resfem liftExcess margins KSmargin FB Weight'

# Inputs variables (separate with whitespace)
input_vars = 'Gama alpha0 d t'

# Log for the Tapenade website (can be any text)
log = 'tapenade test'

# Close session after download?
close_session = False

# Replace push/pops from real4 to real8?
real4toreal8 = True

# Replace diff by the differention mode in the module names?
diff2mode = True

#===================================================

# Append the current directory to the file names
cwd = os.getcwd()
files = list(map(lambda f: os.path.join(cwd,f), files))

# BROWSING SESSION

# Using Firefox to access web
driver = webdriver.Firefox()

# Set implicit wait so that the driver waits until page is loaded
#driver.implicitly_wait(10) # seconds

# Open the website
driver.get('http://www-tapenade.inria.fr:8080/tapenade/form.jsp')

# Uploading files
for f in files:
    browse = driver.find_element_by_name('uploaded_file_or_include')
    browse.send_keys(f)
    button = driver.find_element_by_name('buttonfile')
    button.click()

# Filling boxes
box_names = ['head',           # top routine box
             'depoutputvars',  # output variables box
             'indepinputvars', # input variables box
             'Tapenade-Log']   # application box

box_data = [top_routine,           # top routine box
            output_vars,  # output variables box
            input_vars, # input variables box
            log]   # application box

for name,data in zip(box_names,box_data):

    # Select the box
    box = driver.find_element_by_name(name)

    # Send the data
    box.send_keys(data)

# Differentiate
if mode == 'forward':
    button = driver.find_element_by_name('buttondiffforward')
elif mode == 'reverse':
    button = driver.find_element_by_name('buttondiffbackward')
else:
    raise ValueError('Could not recognize differentiation mode')
button.click()

# Remove previous download file
full_download_filename = os.path.join(download_directory,download_filename)
if os.path.exists(full_download_filename):
    os.remove(full_download_filename)

# Download files
driver.switch_to_frame('diffCG')
form = driver.find_element_by_name('Download')
form.submit()
driver.switch_to_default_content()

# Wait for the download
while (not os.path.exists(full_download_filename)):
    time.sleep(2)
    print('Waiting for ' + full_download_filename)

# Close everything
if close_session:
    driver.quit()
os.remove('geckodriver.log')

# HANDLING DOWNLOADED FILES

# Copy the downloaded file to the current directory
des = os.path.join('.',download_filename)
shutil.move(full_download_filename, des)

# Generate directory name to extract differentiated files
dirname = os.path.splitext(download_filename)[0]
if mode == 'forward':
    dirname = dirname + '_d'
elif mode == 'reverse':
    dirname = dirname + '_b'

# Unpack the downloaded files
os.system('unzip -d %s %s'%(dirname,download_filename))

# Check if we need to replace real4 by real8 in reverse differentiated files
if (mode == 'reverse') and real4toreal8:

    # Print log
    print('Replacing REAL4 by REAL8')

    # Get list of files in the extracted directory
    files = os.listdir(dirname)

    # Loop over the differentiated file names
    for ff in files:

        # Check if this is a differentiated file
        if '_b.f90' in ff:

            # Get full filename
            ff_full = os.path.join(dirname,ff)

            # Replace real4 by real8
            for line in fileinput.input([ff_full], inplace=True):
                print(line.replace('REAL4', 'REAL8'), end='')

# Check if we need to replace diff by the differention mode in module names
if diff2mode:

    # Print log
    print('Replacing DIFF in module names')

    # Get list of files in the extracted directory
    files = os.listdir(dirname)

    # Set file differentiation tag
    if mode == 'forward':
        filetag = '_d.f90'
        modetag = '_D'
    elif mode == 'reverse':
        filetag = '_b.f90'
        modetag = '_B'

    # Loop over the differentiated file names
    for ff in files:

        # Check if this is a differentiated file
        if filetag in ff:

            # Get full filename
            ff_full = os.path.join(dirname,ff)

            # Replace real4 by real8
            for line in fileinput.input([ff_full], inplace=True):
                print(line.replace('_DIFF', modetag), end='')
