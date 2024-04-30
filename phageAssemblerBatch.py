#!/usr/bin/env python3

"""
Arguments:
	fastq_dir      	# Required. Directory with fastq files.
    num_reads		# Optional. Number of reads to try assembling.  Default: 100000
    cutoff_percent	# Optional. Contigs with more than this percentage of assembled reads will be blasted and AceUtiled. Default: 5%

"""

PATH_TO_PHAGEASSEMBLER = "/Users/DAR78/Desktop/Software/GitHub/phageAssemblerUnicycler/phageAssembler.py"

#from datetime import datetime
import subprocess
import argparse
import sys
import os
import re
from shutil import copy, move

# Make parser for options listed above
parser = argparse.ArgumentParser(description='A script to run phageAssembler on multiple phage genomes.',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--fastq_dir', help="Directory with fastq or fastq.gz files to be assembled. One file per genome.")
parser.add_argument('-n','--num_reads', help="The number of reads to try assembling from each file.", type=int, default=100000)
parser.add_argument('-c','--reads_percent_cutoff', help="Contigs with more than this percentage of assembled reads will be blasted and AceUtiled.", type=int, default=5)

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def print_cwd():
    print(os.path.abspath(os.getcwd()))

# Parse command-line options and save values
args = parser.parse_args()
num_reads = args.num_reads
reads_percent_cutoff = args.reads_percent_cutoff

#DEBUG
#print_cwd()

print(f"{bcolors.BOLD}Command: " + " ".join(sys.argv)+f"{bcolors.ENDC}")

fastq_dir = args.fastq_dir
if not os.path.isdir(fastq_dir):
    sys.exit("ERROR: Couldn't find the directory %s, or %s is not a directory." % (fastq_dir,fastq_dir))

print(f"{bcolors.OKBLUE}\tDirectory for assembly: %s{bcolors.ENDC}" % fastq_dir)

# Get list of files with a particular extension
def get_files(dir, extension):
    files = []
    for file in os.listdir(dir):
        if file.endswith(extension):
            files.append(os.path.join(dir, file))
            print("\t\tFound file %s" % file)
    return files

# Unzip fastq.gz files
def gunzip_fastqs(file_list):
    for file in file_list:
        subprocess.run(["gunzip","-f","%s" % file])
        print("\t\tUnzipped file %s" % file)

# Run phageAssembler
def phageAssemble(fastq,nreads,rcutoff):
    print("\t\tCommand: phageAssembler.py -n %s -c %s --quiet %s" % (nreads,rcutoff,fastq))
    result = subprocess.run([PATH_TO_PHAGEASSEMBLER,"-n",str(nreads),"-c",str(rcutoff),"--quiet",fastq])
    if result.returncode == 0:
        print(f"{bcolors.OKGREEN}\t\t...finished running phageAssembler on %s.{bcolors.ENDC}" % fastq)
    else:
        sys.exit(f"{bcolors.FAIL}ERROR: Something went wrong running phageAssembler: %s{bcolors.ENDC}" % result.returncode)

print(f"{bcolors.OKBLUE}\tLooking for *.fastq.gz files in %s...{bcolors.ENDC}" % fastq_dir)
gz_files = get_files(fastq_dir,".fastq.gz")

if gz_files:
    print("\t\tFound %s *.fastq.gz files." % len(gz_files))
    print(f"{bcolors.OKBLUE}\tExpanding *.fastq.gz files...{bcolors.ENDC}")
    gunzip_fastqs(gz_files)
    print("\t\tUnzipped %s *.fastq.gz files." % len(gz_files))
else:
    print("\t\tNo *.fastq.gz files in selected directory.")

print(f"{bcolors.OKBLUE}\tLooking for *.fastq files in %s...{bcolors.ENDC}" % fastq_dir)
fastq_files = get_files(fastq_dir,".fastq")

if not fastq_files:
    sys.exit(f"{bcolors.FAIL}ERROR: No fastq files available for assembly in the directory %s.{bcolors.ENDC}" % fastq_dir)

print("\t\tFound %s *.fastq files." % len(fastq_files))

print(f"{bcolors.OKBLUE}\tWill now run phageAssembler on %s files...{bcolors.ENDC}" % len(fastq_files))
print("\t\tWill use %s reads from each fastq." % num_reads)
print("\t\tWill BLAST and run AceUtil on contigs with over %s%% of reads." % reads_percent_cutoff)

os.chdir(fastq_dir)

g = 0
for fastq in fastq_files:
    print(f"{bcolors.OKBLUE}\tRunning phageAssembler on %s...{bcolors.ENDC}" % fastq)
    phageAssemble(fastq,num_reads,reads_percent_cutoff)
    g = g+1

print(f"{bcolors.OKBLUE}\tFinished running phageAssembler on %s files.{bcolors.ENDC}" % str(g))

print(f"{bcolors.BOLD}All done!{bcolors.ENDC}")

