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
    no-rotate		# Flag. Skips trying to re-orient fasta to similar Base 1s
    quiet			# Flag. Suppresses all output

"""

DEFAULT_NUM_READS = 100000
PATH_TO_ACEUTIL = "~/phageAssemblerUnicycler/AceUtil"
DEFAULT_ADAPTER_LIST = "~/phageAssemblerUnicycler/Adapters/Adapters.fasta"
DEFAULT_BLAST_DATABASE = "~/phageAssemblerUnicycler/BLASTdbs/Actino_DB"

#from datetime import datetime
import subprocess
import argparse
import sys
import os
import re
from shutil import copy, move, which
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
parser.add_argument('--no-rotate', dest='rotate', action='store_false')
parser.add_argument('--quiet', dest='quiet', action='store_true', help="Suppresses output to the terminal but keeps output in log file.")
parser.set_defaults(cutadapt=True, skewer=True, cleanup=True, rotate=True, quiet=False)

cwd = os.path.abspath(os.getcwd())

# Parse command-line options and save values
args = parser.parse_args()

fastq = args.fastq
if not os.path.isfile(fastq):
    sys.exit("ERROR: Couldn't find the file %s" % fastq)
fastq = os.path.abspath(fastq)

cleanup = args.cleanup
num_reads = args.num_reads
reads_percent_cutoff = args.reads_percent_cutoff
cutadapt = args.cutadapt
skewer = args.skewer
rotate = args.rotate
quiet = args.quiet

if args.genome_name:
    genome_name = args.genome_name
else:
    genome_name = re.split('\W+|_', os.path.basename(args.fastq))[0]

if args.adapter_list:
    adapter_list = args.adapter_list
else:
    adapter_list = DEFAULT_ADAPTER_LIST

#Log file
log_file_name = cwd+'/%s_phageAssembler.log' % genome_name
log_file = open(log_file_name,'w')
def printlog(message):
    if not quiet:
        print(message)
    log_file.write(message + "\n")
printlog("Command: " + " ".join(sys.argv))


# Check for required programs
printlog("\n***CHECKING FOR REQUIRED PROGRAMS***")
def check_dependency(name):
    if which(name) is None:
        printlog("\tERROR: Couldn't find %s on your path" % name)
        sys.exit("ERROR: Couldn't find %s on your path" % name)
    else:
        printlog("\t%s found at %s..." % (name,which(name)))

check_dependency("unicycler")
check_dependency("align2Reference.py")
check_dependency("blastn")
check_dependency("skewer")
check_dependency("cutadapt")

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

#Create and move to project directory
subprocess.call(["mkdir", "%s/%s" % (cwd, genome_name)])
project_dir = cwd+"/"+genome_name
os.chdir(project_dir)

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
    if quiet:
        subprocess.call(command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    else:
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
    if quiet:
        subprocess.call('skewer -x %s -q 20 -Q 30 -n -l 50 -o %s -t 8 %s' % (adapter_list,genome_name,filename), shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    else:
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

#move(assembly_fastq,".")

#Assemble with unicycler
def run_unicycler(fastq):
    unicycler_command = "unicycler -s %s -o unicycler_assembly --no_rotate --kmers 87,115,121,127" % fastq
    printlog("\tRunning: %s" % unicycler_command)
    if quiet:
        subprocess.call(unicycler_command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    else:
        subprocess.call(unicycler_command, shell=True)

printlog("\n***ASSEMBLING READS WITH UNICYCLER***")
run_unicycler(assembly_fastq)
printlog("\tUnicycler assembly complete. Details in /unicycler_assembly/unicycler.log")

def align_reads(fasta,fastq):
    command = "align2Reference.py --fasta %s --fastq %s" % (fasta, fastq)
#    printlog("\tWill run: %s" % newbler_command)
    if quiet:
        subprocess.call(command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    else:
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
    if quiet:
        subprocess.call(blast_command,shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    else:
        subprocess.call(blast_command,shell=True)
    return outfile

def biopy_blast(queryseq,database,outfile="blast_output.xml"):
    blast_command = "blastn -db %s -query %s -out %s -outfmt 5" % (database,queryseq,outfile)
    if quiet:
        subprocess.call(blast_command,shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    else:
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

def blast_cleanup():
    subprocess.call("rm *.fasta", shell=True)
    subprocess.call("rm *blast.xml", shell=True)
#    subprocess.call(["rm","*.fasta"])
#    subprocess.call(["rm","*blast.xml"])

os.chdir(project_dir)

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

#blast_cleanup()
if cleanup:
    printlog("\tRemoving blast files...")
    blast_cleanup()
    printlog("\t\t...blast files removed.")

#Reorient to Base 1 if possible
def reorient_fasta(fasta,complement,position):
    seq = SeqIO.parse(open(fasta,'r'),'fasta')
    s = list(seq)[0]
    if complement:
        s2 = s.reverse_complement()
        s2.id = s.id
        s2.name = s.name
        s2.description = s.description
        s = s2
    return s[position-1:] + s[:position-1]

reorient = False
if rotate:
    printlog("\n***ATTEMPTING TO ORIENT TO BASE 1***")
    with open('%s/%s/unicycler_assembly/assembly.gfa' % (cwd,genome_name)) as f:
        for line in f:
            pass
        last_line = line.rstrip().split("\t")
    if len(contigs) != 1:
        printlog("\tMore than one contig was assembled, so skipping Base 1 orientation.")
    elif (len(last_line) > 3 and last_line[0] != "L" and last_line[1] != "1" and last_line[5] != "0M"):
        printlog("\tOnly one contig was assembled, but ends don't cleanly link, so skipping Base 1 orientation.")
    elif (not base_ones):
        printlog("\tNo BLAST matches to known Base 1s, so skipping Base 1 orientation.")
    else:
        reorient = True
        printlog("\tThe assembly produced a single contig, and its ends link together. Great!")
        printlog("\tLet's try to reorient the genome to match Base 1s of similar genomes...")
        newseq = reorient_fasta('%s/unicycler_assembly/assembly.fasta' % project_dir,base_ones[0][3],base_ones[0][1])
        rc_name = '%s_reoriented.fasta' % genome_name
        SeqIO.write(newseq, rc_name, 'fasta')
        printlog("\t\t...grabbing fastq from old consed/solexa_dir.")
        move("%s/consed/solexa_dir/%s" % (project_dir,assembly_fastq), project_dir)
        printlog("\t\t...removing previous consed directory.")
        subprocess.call(("rm -r %s/consed" % project_dir), shell=True)
        printlog("\t\t...previous consed dir removed.")
        printlog("\t\t...creating new consed directory.")
        subprocess.call(["mkdir", "%s/consed" % project_dir])
        printlog("\t\t...created new consed directory.")
        move("%s/%s" %(project_dir,assembly_fastq),"%s/consed/" % project_dir)
        printlog("\t\t...moved fastq to new consed directory.")
        move(rc_name,"%s/consed" % project_dir)
        printlog("\t\t...moved reoriented fasta to new consed directory.")
        os.chdir("./consed/")
        printlog("\t\t...aligning reads to reoriented consensus.")
        align_reads(rc_name,"%s" % assembly_fastq)
        printlog("\t\t...reads aligned and consed files created.")
        
        # Remove old consed dir, make new consed dir with newly aligned deal

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
    if quiet:
        subprocess.call(command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    else:
        print("\t\tCommand: %s" % command)
        subprocess.call(command, shell=True)
    return outfile

printlog("\n***ACE UTIL***")
if reorient:
    aceutil_infile = "%s/consed/edit_dir/%s_reoriented.ace.1" % (project_dir,genome_name)
    printlog("\tRunning AceUtil on lone contig...")
    log_file.close()
    run_AceUtil(aceutil_infile,contig="1")
    log_file = open(log_file_name,'a')  
else:
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
if reorient:
    printlog("\t\tSince there was a single contig with a blastn match to Base 1 of another genome, the contig was reoriented to match that genome.")
    printlog("\t\tThat means hopefully the correct Base 1 is at the first base of the current contig in the consed folder.")
elif base_ones:
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
if not reorient:
    all_contig_objects=SeqIO.parse(open('./unicycler_assembly/assembly.fasta','r'),'fasta')
else:
    all_contig_objects=SeqIO.parse(open('%s/consed/edit_dir/%s_reoriented.fasta' % (project_dir,genome_name),'r'),'fasta')
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
move(log_file_name,project_dir)
