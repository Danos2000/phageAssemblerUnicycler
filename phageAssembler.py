#!/usr/bin/env python3

"""
Arguments:
	fastq      		# Required. The file of reads in fastq format.
	genome_name		# Optional. Specify the name of this genome. Default will be the input fastq file name before any special characters.
    num_reads		# Optional. Number of reads to try assembling.  Default: 100000
    adapter_list	# Optional. Specify a file of adapters to use.  Default: Adapters.fasta included in package
    cutoff_percent	# Optional. Contigs with more than this percentage of assembled reads will be blasted and AceUtiled. Default: 5%
    no-cutadapt		# Flag. Skips NextSeq-specific cutadapt trimming.
    no-skewer		# Flag. Skips skewer trimming.
    no-clean		# Flag. Skips cleanup steps, resulting in more output files.

"""

PATH_TO_UNICYCLER = "/usr/local/genome/bin/"
PATH_TO_ACEUTIL = "~/Desktop/GitHub/phageAssemblerUnicycler/AceUtil"
PATH_TO_ALIGN2REF = "/usr/local/genome/bin/"
DEFAULT_NUM_READS = 100000
DEFAULT_ADAPTER_LIST = "~/Desktop/GitHub/phageAssemblerUnicycler/Adapters/Adapters.fasta"
DEFAULT_BLAST_DATABASE = "~/Desktop/GitHub/phageAssemblerUnicycler/BLASTdbs/Actino_DB"

#from datetime import datetime
import subprocess
import argparse
import sys
import os
import re
from shutil import copy, move
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from Bio.Sequencing import Ace
from Bio.Blast import NCBIXML

# Make parser for options listed above
parser = argparse.ArgumentParser(description='A script to assemble phage genomes.',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('fastq', help="A file of reads in fastq format.")
parser.add_argument('-g','--genome_name', help="The name of the genome you are assembling.  Default will be name of the fastq before any special characters.")
parser.add_argument('-n','--num_reads', help="The number of reads to try assembling.", type=int, default=DEFAULT_NUM_READS)
parser.add_argument('-a','--adapter_list', help="A fasta-formatted file of adapter sequences to be trimmed. Default file contains common adapters from Illumina, NEB, and more.", default=DEFAULT_ADAPTER_LIST)
parser.add_argument('-c','--reads_percent_cutoff', help="Contigs with more than this percentage of assembled reads will be blasted and AceUtiled.", type=int, default=5)
parser.add_argument('--no-cutadapt', dest='cutadapt', action='store_false')
parser.add_argument('--no-skewer', dest='skewer', action='store_false')
parser.add_argument('--no-cleanup', dest='cleanup', action='store_false')
parser.set_defaults(cutadapt=True, skewer=True, cleanup=True)

# Parse command-line options and save values
args = parser.parse_args()

fastq = args.fastq
if not os.path.isfile(fastq):
    sys.exit("ERROR: Couldn't find the file %s" % fastq)

cleanup = args.cleanup
num_reads = args.num_reads
reads_percent_cutoff = args.reads_percent_cutoff
cutadapt = args.cutadapt
skewer = args.skewer

if args.genome_name:
    genome_name = args.genome_name
else:
    genome_name = re.split('\W+|_', os.path.basename(args.fastq))[0]

if args.adapter_list:
    adapter_list = args.adapter_list
else:
    adapter_list = DEFAULT_ADAPTER_LIST

#Log file
log_file_name = '%s_phageAssembler.log' % genome_name
log_file = open(log_file_name,'w')
def printlog(message):
    print(message)
    log_file.write(message + "\n")
printlog("Command: " + " ".join(sys.argv))

#Input
printlog("\n***INPUT AND SETTINGS***")
printlog("\tGenome: %s" % genome_name)
printlog("\tWill begin with %s reads from the file %s" % (str(num_reads), fastq))
if skewer:
    printlog("\tAdapters file in fasta format for skewer trimming: %s" % adapter_list)
printlog("\tWill BLAST and run aceUtil on contigs with over %s%% of reads" % str(reads_percent_cutoff))
if cutadapt:
    printlog("\tWill perform NextSeq-specific trimming with cutadapt.")
else:
    printlog("\tWill NOT perform NextSeq-specific trimming with cutadapt since --no-cutadapt was invoked.")
if skewer:
    printlog("\tWill perform adapter and quality trimming with skewer.")
else:
    printlog("\tWill NOT perform adapter and quality trimming with skewer since --no-skewer was invoked.")
if cleanup:
    printlog("\tWill delete some intermediate files since --no-cleanup wasn't invoked.")
else:
    printlog("\tWill NOT delete intermediate files since --no-cleanup was invoked.")

def wc(filename):
    from subprocess import check_output
    return int(check_output(["wc", "-l", filename]).split()[0])

# Total reads in initial fastq file
total_reads = int(wc(fastq)/4)

#Subset of reads function
def subsample_fastq_file(filename, number_of_reads, head_tail="head", new_file_if_all_reads=True):
    reads = number_of_reads
    if total_reads < number_of_reads:
        reads = total_reads
        printlog("\tFewer reads in file than total requested.")
        if new_file_if_all_reads: 
            new_filename = '%s_AllReads.fastq' % (genome_name)
            copy(filename,new_filename)
            printlog("\tCreated new file, %s, with all reads for assembly." % new_filename)
            return new_filename, reads
        else:
            printlog("\tWill use entire original file, %s, for assembly." % filename)
            return filename, reads
    new_filename = '%s_%skReads.fastq' % (genome_name, str(int(number_of_reads/1000)))
    subprocess.call('head -n %s %s > %s' % (str(number_of_reads*4),filename,new_filename), shell=True)
    printlog("\tCreated new file %s with %s reads." % (new_filename, str(number_of_reads)))
    return new_filename, reads

#Function to run cutadapt
def run_cutadapt(filename):
    new_filename = os.path.splitext(filename)[0]+'_cutadapt.fastq'
    command = 'cutadapt --nextseq-trim 30 -o %s %s > cutadapt.log' % (new_filename,filename)
    subprocess.call(command, shell=True)
    printlog("\t\t...cutadapt trimming complete, output in %s" % new_filename)
    if cleanup:
        printlog("\tRemoving %s and cutadapt.log to tidy up..." % filename)
        subprocess.call(["rm","%s" % filename])
        subprocess.call(["rm","cutadapt.log"])
        printlog("\t\t...files removed.")
    return new_filename

#Function to run skewer
def run_skewer(filename):
    subprocess.call('skewer -x %s -q 20 -Q 30 -n -l 50 -o %s -t 8 %s' % (adapter_list,genome_name,filename), shell=True)
    new_filename = genome_name+'-trimmed.fastq'
    printlog("\t\t...skewer trimming complete, output in %s" % new_filename)
    os.rename(genome_name+'-trimmed.log','skewer.log')
    if cleanup:
        printlog("\tRemoving %s and skewer.log to tidy up..." % filename)
        subprocess.call(["rm","%s" % filename])
        subprocess.call(["rm","skewer.log"])
        printlog("\t\t...files removed.")
    return new_filename

#Prepare assembly fastq by downsampling and trimming
printlog("\n***PREPARING FASTQ FOR ASSEMBLY***")
assembly_fastq,num_reads = subsample_fastq_file(fastq, num_reads)
if cutadapt:
    printlog("\tRemoving NextSeq artifacts with cutadapt...")
    assembly_fastq = run_cutadapt(assembly_fastq)
if skewer:
    printlog("\tTrimming reads with skewer...")
    assembly_fastq = run_skewer(assembly_fastq)
printlog("\tFinal fastq file for assembly: %s" % assembly_fastq)

#Assemble with unicycler
def run_unicycler(fastq):
    unicycler_command = "%s/unicycler -s %s -o unicycler_assembly --no_rotate" % (PATH_TO_UNICYCLER, fastq)
    printlog("\tRunning: %s" % unicycler_command)
    subprocess.call(unicycler_command, shell=True)

printlog("\n***ASSEMBLING READS WITH UNICYCLER***")
run_unicycler(assembly_fastq)
printlog("\tUnicycler assembly complete. Details in /unicycler_assembly/unicycler.log")

def align_reads(fasta,fastq):
    command = "%s/align2Reference.py --fasta %s --fastq %s" % (PATH_TO_ALIGN2REF, fasta, fastq)
#    printlog("\tWill run: %s" % newbler_command)
    subprocess.call(command, shell=True)

printlog("\n***CREATING CONSED DIRECTORY STRUCTURE AND ALIGNING READS***")
subprocess.call(["mkdir", "consed"])
printlog("\tCreated consed directory.")
move(assembly_fastq,"./consed/")
printlog("\tMoved fastq to consed directory.")
copy("./unicycler_assembly/assembly.fasta","./consed/UNIout.fasta")
printlog("\tMoved unicycler output fasta to consed directory.")
os.chdir("./consed/")
printlog("\tAligning reads to consensus...")
align_reads("UNIout.fasta","%s" % assembly_fastq)
printlog("\t\t...reads aligned and consed files created.")

printlog("\n***CHECKING ASSEMBLY RESULTS***")
ace = Ace.read(open('edit_dir/UNIout.ace.1'))
printlog("\t%s total contigs were assembled. They are:" % str(ace.ncontigs))
contigs = ace.contigs
printlog("\t\tName\tLength (bp)\t\t# Reads\t\t~% of Initial Reads")
contigs_to_blast=[]
c = 1
for contig in contigs:
    if c > 20:
        printlog("\tOnly showing the largest 20 contigs.")
        break
    percent_reads = round(100*(contig.nreads/num_reads),1)
    out_string = "\t\t%s\t\t%s\t\t%s\t\t%s" % (contig.name, contig.nbases, contig.nreads, str(percent_reads))
    if percent_reads >= reads_percent_cutoff:
        out_string += "*"
        contigs_to_blast.append(contig)
    printlog(out_string)
    c = c+1
printlog("\t* These contigs have > %s%% of the initial reads and will thus be blasted and run through AceUtil." % reads_percent_cutoff)

def blast_contigs(seqfile,database,outfile="blast_results.txt"):
    blast_command = "blastn -db %s -query %s -out %s" % (database,seqfile,outfile)
    subprocess.call(blast_command,shell=True)
    return outfile

def biopy_blast(queryseq,database,outfile="blast_output.xml"):
    blast_command = "blastn -db %s -query %s -out %s -outfmt 5" % (database,queryseq,outfile)
    subprocess.call(blast_command,shell=True)
    result_handle = open(outfile,'r')
    return NCBIXML.read(result_handle)

def display_blast_results(record, reblast=False):
    potential_clusters = []
    base_one = ()
    printlog("\n\n\tQuery: %s" % record.query)
    if len(record.alignments) > 10:
        printlog("\t\tFirst 10 hits:")
    else:
        printlog("\t\tAll hits:")
    if len(record.alignments) == 0:
        printlog("\t\t\tNo hits found.")
    else:
        for alignment in record.alignments[:10]:
            title_split=alignment.title.split()
            outline='\t\t\t' + ' '.join([title_split[1],title_split[2],title_split[3]]) + ' (%s bp)' % str(alignment.length)
            if title_split[-2] == 'Cluster':
                outline += ", Cluster %s" % title_split[-1]
                potential_clusters.append(title_split[-1])
            printlog(outline)
    printlog("\t\tBest hit details:")
    try:
        best = record.alignments[0]
        title_split=best.title.split()
        best_title = ' '.join([title_split[1],title_split[2],title_split[3]]) + ' (%s bp)' % str(best.length)
        printlog("\t\t\tSubject: %s" % best_title)
        i = 1
        for hsp in best.hsps[:10]:
            printlog("\t\t\t\tHit #%s" % str(i))
            printlog("\t\t\t\tLength: %s bp Score: %s E-value: %s" % (hsp.align_length,hsp.score,hsp.expect))
            printlog("\t\t\t\tIdentities: %s/%s Gaps: %s/%s." % (hsp.identities,hsp.align_length,hsp.gaps,hsp.align_length))
            if hsp.frame[1] == -1:
                strand = "Minus"
            else:
                strand = "Plus"
            printlog("\t\t\t\tStrand: Plus/%s" % strand)
            printlog("\t\t\t\tQuery: %s" % hsp.query_start)
            printlog("\t\t\t\tSbjct: %s" % hsp.sbjct_start)
            if hsp.sbjct_start == 1:
                base_one = (best_title,hsp.query_start,hsp.sbjct_start,reblast,record.query,)
            i += 1
    except:
        pass
    return (potential_clusters, base_one,)

def parse_blast(blastresults):
    best_hits = []
    with open(blastresults,'r') as in_file:
        for line in in_file:
            if (len(best_hits) <= 10) and " phage " in line:
                data = line.split()
                best_hits.append((data[2],data[-2],data[-1],))
    return best_hits

def rc_reblast(seq_file):
    out_file = '%s_rc.fasta' % seq_file.split('.')[0]
    rc_out = open(out_file,'w')
    seq = SeqIO.read(open(seq_file,'r'),'fasta')
    rc_out.write(">%s_reverse_complement\n" % seq.name)
    seq = seq.reverse_complement()
    rc_out.write(str(seq.seq))
    rc_out.close()
    return biopy_blast(out_file, DEFAULT_BLAST_DATABASE, outfile='%s_blast.xml' % out_file.split('.')[0])

os.chdir("..")

#BLAST
printlog("\n***BLAST***")
printlog("\tRunning local blast of %s contig(s) against %s database..." % (str(len(contigs_to_blast)),DEFAULT_BLAST_DATABASE))

all_contig_objects=SeqIO.parse(open('./unicycler_assembly/assembly.fasta','r'),'fasta')
blasted_contigs = []
contig_names = [x.name for x in contigs_to_blast]
for contig in all_contig_objects:
    if contig.id in contig_names:
        SeqIO.write(contig, '%s.fasta' % contig.id, 'fasta')
        blasted_contigs.append(biopy_blast('%s.fasta' % contig.id, DEFAULT_BLAST_DATABASE, outfile='%s_blast.xml' % contig.id))

reblasted_contigs = []
cluster_guesses = []
base_ones = []
for result in blasted_contigs:
    cg = display_blast_results(result)
    cluster_guesses.append((result.query.split()[0],cg[0],))
    if cg[1]:
        base_ones.append(cg[1])
    try:
        if result.alignments[0].hsps[0].frame[1] == -1:
            reblasted_contigs.append(rc_reblast('%s.fasta' % result.query.split()[0]))
    except:
        pass

if reblasted_contigs:
    printlog("\n\tRe-blasting %s contig(s) in reverse orientation." % str(len(reblasted_contigs)))
    for result in reblasted_contigs:
        cg = display_blast_results(result, reblast=True)
        if cg[1]:
            base_ones.append(cg[1])

#AceUtil
def run_AceUtil(acefile,contig=None):
    try:
        outfile = acefile.rsplit('.',1)[0] + "." + str(int(acefile.rsplit('.',1)[1])+1)
    except:
        outfile = acefile + ".aceUtil"
    command = "java -jar %s/AceUtil.jar %s %s " % (PATH_TO_ACEUTIL, acefile, outfile)
    if contig:
        command += contig
    command += " >> %s" % (log_file_name)
    print("\t\tCommand: %s" % command)
    subprocess.call(command, shell=True)
    return outfile

printlog("\n***ACE UTIL***")
aceutil_infile = "./consed/edit_dir/UNIout.ace.1"
for contig in contigs_to_blast:
    printlog("\tRunning AceUtil on %s..." % contig.name)
    log_file.close()
    aceutil_outfile = run_AceUtil(aceutil_infile,contig=contig.name)
    aceutil_infile = aceutil_outfile
    log_file = open(log_file_name,'a')
printlog("\tAceUtil analysis complete.")

#Report
printlog("\n***REPORT***")
printlog("\tCluster Guess")
if cluster_guesses:
    for contig in cluster_guesses:
        gs = ', '.join(contig[1][:5])
        printlog("\t\t%s\tCluster of top hits: %s" % ("Contig_"+contig[0], gs))
    if len(set(contig[1][:5])) == 1:
        printlog("\t\tProbable cluster: %s" %  contig[1][0])
    else:
        printlog("\t\tUnable to make single cluster guess from blast results.")
else:
    printlog("\t\tUnable to determine a likely cluster.")

printlog("\n\tBase One Guess")
if base_ones:
    for base_one in base_ones:
        out = "\t\tIn the blast hit to %s, query position %s matches subject position %s." % (base_one[0], str(base_one[1]), str(base_one[2]))
        out2 = "\t\tLikely Base 1 position: %s in %s" % (base_one[1], base_one[4]) 
        if base_one[3]:
            out += "  (After contig was reverse-complemented.)"
        printlog(out)
        printlog(out2)
else:
    printlog("\t\tUnable to find Base 1.")

printlog("\n\tGC Info")
i=0
all_contig_objects=SeqIO.parse(open('./unicycler_assembly/assembly.fasta','r'),'fasta')
for contig in all_contig_objects:
    if i>=10:
        printlog("\t\tOnly showing first 10 contigs.")
        break
    printlog("\t\t%s\t%s %%" % ("Contig_"+contig.id, round(100*(gc_fraction(contig.seq)),1)))
    i += 1

printlog("\n\tCoverage Info")
i=0
for contig in contigs:
    if i>=10:
        printlog("\t\tOnly showing first 10 contigs.")
        break
    total_bases = 0
    for read in contig.reads:
        total_bases=total_bases+len(read.rd.sequence)
    assembled_coverage = int(total_bases/contig.nbases)
    total_coverage = int((assembled_coverage/num_reads)*total_reads)
    printlog("\t\t%s\t%s (assembled)\t\t%s (estimated for entire initial fastq)" % ("Contig_"+contig.name,str(assembled_coverage),str(total_coverage))) 
    i += 1

log_file.close()

