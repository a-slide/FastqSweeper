#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-

"""
@package    FastqSweeper
@brief      Main file of the program
@copyright  [GNU General Public License v2](http://www.gnu.org/licenses/gpl-2.0.html)
@author     Adrien Leger - 2015
* <adrien.leger@gmail.com>
* <adrien.leger@inserm.fr>
* <adrien.leger@univ-nantes.fr>
* [Github](https://github.com/a-slide)
* [Atlantic Gene Therapies - INSERM 1089] (http://www.atlantic-gene-therapies.fr/)
"""

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# Standard library packages
from os import path, remove
import optparse
import sys
from time import time
from subprocess import Popen, PIPE

# Third party package
import pysam

# Local packages


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class FastqSweeper (object):
    """
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~CLASS FIELDS~~~~~~~#

    VERSION = "FastqSweeper 0.1"
    USAGE = "Usage: %prog -i INDEX -1 FASTQ_R1 [-r -2 FASTQ_R2 -t THREAD -b BWA_OPTIONS -c CUTADAPT_OPTIONS -a ADAPTER]"

    #~~~~~~~CLASS METHODS~~~~~~~#

    @classmethod
    def class_init (self):
        """
        init class method for instantiation from command line. Parse arguments parse CL arguments
        """

        print("\nParse commande line arguments")

        ### Define parser usage, options
        optparser = optparse.OptionParser(usage = self.USAGE, version = self.VERSION)

        optparser.add_option('-i', '--index', dest="index", help= "bwa index path (required)")
        optparser.add_option('-1', '--fastq_R1', dest="fastq_R1", help= "Path to the fastq file (required)")
        optparser.add_option('-2', '--fastq_R2', dest="fastq_R2", help= "Path to the pair fastq file if paired end (facultative)")
        optparser.add_option('-t', '--thread', dest="thread", default=1, help= "Number of thread to use (default: 1)")
        optparser.add_option('-b', '--bwa_opt', dest="bwa_opt", help= "bwa options for the mapping step (facultative and quoted)")
        optparser.add_option('-c', '--cutadapt_opt', dest="cutadapt_opt", default= "-m 25 -q 30,30 --trim-n",
            help= "cutadapt options for the qc step (facultative and quoted) (default: -m 25 -q 30,30 --trim-n)")
        optparser.add_option('-a', '--adapter', dest="adapter", help= "Path to a fasta file containing adapters to be 3' trimmed (facultative)")
        optparser.add_option('-r', '--run', dest="run", action="store_true", default=False, help= "Run command lines (Default print command lines without running)")

        ### Parse arguments
        opt, args = optparser.parse_args()

        ### Init a RefMasker object
        return FastqSweeper (opt.index, opt.fastq_R1, opt.fastq_R2, int(opt.thread), opt.bwa_opt, opt.cutadapt_opt, opt.adapter, opt.run)

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __init__(self, index=None, fastq_R1=None, fastq_R2=None, thread=1, bwa_opt=None, cutadapt_opt=None, adapter=None, run=False):
        """
        General initialization function for import and command line
        """

        print("\nInitialize FastqSweeper")

        ### Verifications
        assert index, "A path to the bwa index is mandatory"
        assert fastq_R1, "A path to the fastq file is mandatory"

        ### Storing Variables
        self.index = index
        self.fastq_R1 = fastq_R1
        self.fastq_R2 = fastq_R2
        self.thread = int(thread)
        self.bwa_opt = bwa_opt
        self.cutadapt_opt = cutadapt_opt
        self.adapter = adapter
        self.run = run
        self.single = False if fastq_R2 else True
        self.basename = fastq_R1.rpartition('/')[2].partition('.')[0]

        ## ADDITIONAL VARIABLES
        self.process_mapped = True
        self.process_unmapped = True
        self.min_mapq = 30
        self.min_size = 25


    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def __call__(self):
        """
        """

        start_time = time()

        if self.run:
            print ("\nRunning in real mode")
        else:
            print ("\nRunning in safe test mode")

        if self.single:
            print ("\nProcessing fastq in single-end mode ")
            self.process_single_end()
        else:
            print ("\nProcessing fastq in paired-end mode ")
            self.process_paired_end()

        # Generate Reports

        # Finalize
        print ("\nDone in {}s".format(round(time()-start_time, 3)))
        return(0)


    def process_single_end (self):

        ### Cutadapt
        print ("\nStarting trimming with CUTADAPT")

        trimmed_fastq = self.basename+"_trim.fastq.gz"
        cutadapt_report = self.basename+"_trim_report.txt"

        cmd = "cutadapt {} {} {} -o {} > {}".format (self.cutadapt_opt, \
        "-a file:"+self.adapter if self.adapter else "", self.fastq_R1, trimmed_fastq, cutadapt_report)

        if self.run:
            proc = Popen(cmd, shell=True)
            proc.communicate()[0]
        else:
            print (cmd)

        ### BWA alignment, compression and sorting
        print ("\nStart aligning with BWA MEM, sort and compress")

        #sam = self.basename+"_all.sam"
        bam = self.basename+"_all.bam"
        cmd1 = "bwa mem -M -t {} {} {}".format (self.thread, self.index, trimmed_fastq)
        cmd2 = "samtools view - -hb "
        cmd3 = "samtools sort - -o a > {}".format ( bam)

        if self.run:
            p1 = Popen(cmd1, stdout=PIPE, shell=True)
            p2 = Popen(cmd2, stdin=p1.stdout, stdout=PIPE, shell=True)
            p3 = Popen(cmd3, stdin=p2.stdout, stdout=PIPE, shell=True)
            p1.stdout.close()  # Allow p1 to receive a SIGPIPE if p2 exits.
            p2.stdout.close()  # Allow p2 to receive a SIGPIPE if p3 exits.
            print p3.communicate()[0]
        else:
            print ("{} | {} | {}".format (cmd1, cmd2, cmd3))

        ## Post processing of reads
        print ("\nStart sorting reads")

        if not self.run:
            print ("Done")
            return

        with pysam.AlignmentFile(bam, "rb") as bamfile:

            # If mapped reads are to be processed init a bam file
            if self.process_mapped:
                pass
                #self.bam_header = bamfile.header

            # If unmapped reads are to be processed init a fastq file
            if self.process_unmapped:
                pass

            # Init a
            count = {"secondary":0, "unmapped":0, "lowmapq":0, "mapped":0, "short_mapped":0 }

            # Parse reads
            for read in bamfile:

                # Always remove secondary alignments
                if read.is_secondary:
                    count["secondary"] +=1

                # Extract mapped read
                elif read.mapq >= self.min_mapq:
                    count["mapped"] += 1
                    if self.process_mapped:
                        pass
                        # Create bam and bedgraph

                # Regenerate fastq from unmapped reads
                # Consider short match, low mapq score and unmapped reads as unmapped
                else:
                    if read.tid == -1:
                        count["unmapped"] += 1

                    elif len(read.query_alignment_sequence) < self.min_size:
                        count["short_mapped"] += 1

                    else: # not unmapped but mapq < min_mapq
                        count["lowmapq"] +=1

                    if self.process_unmapped:
                        pass
                        # Regenerate fastq

        print (count)

        #Removing the original bam file which is no longer needed
        #remove(bam)


    def process_paired_end (self):
        pass
        # cutadapt

        # bwa alignment

        # bwa post processing of mapped reads

        # bwa post processing of unmapped reads

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#   TOP LEVEL INSTRUCTIONS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

if __name__ == '__main__':

    fastq_sweeper = FastqSweeper.class_init()
    fastq_sweeper()
