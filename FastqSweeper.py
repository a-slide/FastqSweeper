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
from gzip import open as gopen
from collections import OrderedDict
from datetime import datetime

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
    USAGE = "Usage: %prog -i INDEX -1 FASTQ_R1 [-r -2 FASTQ_R2 -t INT -b BWA_OPT -c CUTADAPT_OPT -a ADAPTER -m INT -s INT]"

    #~~~~~~~CLASS METHODS~~~~~~~#

    @classmethod
    def class_init (self):
        """
        init class method for instantiation from command line. Parse arguments parse CL arguments
        """

        print("\nParse commande line arguments")

        ### Define parser usage, options
        optparser = optparse.OptionParser(usage = self.USAGE, version = self.VERSION)

        optparser.add_option('-i', '--index', dest="index",
            help= "bwa index path (required)")
        optparser.add_option('-1', '--fastq_R1', dest="fastq_R1",
            help= "Path to the fastq file (required)")
        optparser.add_option('-2', '--fastq_R2', dest="fastq_R2",
            help= "Path to the pair fastq file if paired end (facultative)")
        optparser.add_option('-t', '--thread', dest="thread", default=1,
            help= "Number of thread to use (default: 1)")
        optparser.add_option('-b', '--bwa_opt', dest="bwa_opt", default= "-M",
            help= "bwa options for the mapping step (facultative and quoted) (default: -M)")
        optparser.add_option('-c', '--cutadapt_opt', dest="cutadapt_opt", default= "-m 25 -q 30,30 --trim-n",
            help= "cutadapt options for the qc step (facultative and quoted) (default: -m 25 -q 30,30 --trim-n)")
        optparser.add_option('-a', '--adapter', dest="adapter",
            help= "Path to a fasta file containing adapters to be 3' trimmed (facultative)")
        optparser.add_option('-r', '--run', dest="run", action="store_true", default=False,
            help= "Run command lines (default: False)")
        optparser.add_option('-m', '--min_mapq', dest="min_mapq", default= 0,
            help= "Minimal mapq quality to be considered mapped (default: 0)")
        optparser.add_option('-s', '--min_match_size', dest="min_match_size", default= 0,
            help= "Minimal size of match to be considered mapped (default: 0)")

        ### Parse arguments
        opt, args = optparser.parse_args()

        ### Init a RefMasker object
        return FastqSweeper (
            index=opt.index,
            fastq_R1=opt.fastq_R1,
            fastq_R2=opt.fastq_R2,
            thread=int(opt.thread),
            bwa_opt=opt.bwa_opt,
            cutadapt_opt=opt.cutadapt_opt,
            adapter=opt.adapter,
            run=opt.run,
            min_mapq=int(opt.min_mapq),
            min_match_size=int(opt.min_match_size))

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __init__(self,
        index=None,
        fastq_R1=None,
        fastq_R2=None,
        thread=1,
        bwa_opt="-M",
        cutadapt_opt="-m 25 -q 30,30 --trim-n",
        adapter=None,
        run=False,
        min_mapq=0,
        min_match_size=0):
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
        self.min_mapq = int(min_mapq)
        self.min_size = int(min_match_size)
        self.single = False if fastq_R2 else True
        self.basename = fastq_R1.rpartition('/')[2].partition('.')[0]


    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def __call__(self):
        """
        General function launching either the single end or paired end fonction
        """

        start_time = time()

        if self.run:
            print ("\nRunning in real mode")
        else:
            print ("\nRunning in safe test mode")

        if self.single:
            print ("\nProcessing fastq in single-end mode ")
            count = self.process_single_end()
        else:
            print ("\nProcessing fastq in paired-end mode ")
            count = self.process_paired_end()

        # Generate Reports

        if self.run:
            print ("\nGenerate a Report")
            with open (self.basename+"_FastqSweeper_report.csv", "w") as report:
                report.write ("Program {}\tDate {}\n".format(self.VERSION,str(datetime.today())))
                report.write ("\nRUN PARAMETERS\n")
                report.write("  Index basename\t{}\n".format(self.index))
                report.write("  Fastq R1 path\t{}\n".format(self.fastq_R1))
                report.write("  Fastq R2 path\t{}\n".format(self.fastq_R2 if self.fastq_R2 else "None"))
                report.write("  Number of thread\t{}\n".format(self.thread))
                report.write("  Bwa options\t{}\n".format(self.bwa_opt))
                report.write("  Cutadapt options\t{}\n".format(self.cutadapt_opt))
                report.write("  Adapter file\t{}\n".format(self.adapter if self.adapter else "None"))
                report.write("  Cutoff MapQ score\t{}\n".format(self.min_mapq))
                report.write("  Cutoff match length\t{}\n".format(self.min_size))
                report.write("  Mode used\t{}\n".format("single-end" if self.single else "paired-end"))
                report.write("  Output files basename\t{}\n".format(self.basename))
                report.write("\nREAD COUNT PER CATEGORY\n")
                report.write("  Mapped\t{}\n".format(count["mapped"]))
                report.write("  Unmapped\t{}\n".format(count["unmapped"]+count["short_mapped"]+count["lowmapq"]))
                report.write("  Total\t{}\n\n".format(count["mapped"]+count["unmapped"]+count["short_mapped"]+count["lowmapq"]))
                report.write("  Unaligned\t{}\n".format(count["unmapped"]))
                report.write("  Short match\t{}\n".format(count["short_mapped"]))
                report.write("  Low MapQ\t{}\n".format(count["lowmapq"]))
                report.write("  Secondary\t{}\n".format(count["secondary"]))

        # Finalize
        print ("\nDone in {}s".format(round(time()-start_time, 3)))
        return(0)


    def process_single_end (self):

        # Generates names for output files
        trimmed_fastq = self.basename+"_trim.fastq.gz"
        cutadapt_report = self.basename+"_trim_report.txt"
        all_bam = self.basename+"_all.bam"
        mapped_bam = self.basename+"_mapped.bam"
        unmapped_fastq = self.basename+"_clean.fastq.gz"

        ### Cutadapt
        print ("\nStarting trimming with CUTADAPT")

        cmd = "cutadapt {} {} {} -o {} > {}".format (self.cutadapt_opt, \
        "-a file:"+self.adapter if self.adapter else "", self.fastq_R1, trimmed_fastq, cutadapt_report)

        if self.run:
            proc = Popen(cmd, shell=True)
            proc.communicate()[0]
        else:
            print (cmd)

        ### BWA alignment, compression and sorting
        print ("\nStart aligning with BWA MEM, sort and compress")

        cmd1 = "bwa mem {} -t {} {} {}".format (self.bwa_opt, self.thread, self.index, trimmed_fastq)
        cmd2 = "samtools view - -hb "
        cmd3 = "samtools sort - -o a > {}".format (all_bam)

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

        # Open input and output files within the same context manager block
        with \
            pysam.AlignmentFile(all_bam, "rb") as all_bam_h, \
            pysam.AlignmentFile(mapped_bam, "wb", header=all_bam_h.header) as mapped_bam_h, \
            gopen (unmapped_fastq, "w") as unmapped_fastq_h:

            # Init a dict of counters
            count = {"secondary":0, "unmapped":0, "lowmapq":0, "mapped":0, "short_mapped":0 }

            # Parse reads
            for read in all_bam_h:

                # Always remove secondary alignments
                if read.is_secondary:
                    count["secondary"] +=1

                # Extract mapped read
                elif read.tid != -1 and read.mapq >= self.min_mapq:
                    count["mapped"] += 1
                    mapped_bam_h.write(read)

                # Regenerate fastq from unmapped reads
                # Consider short match, low mapq score and unmapped reads as unmapped
                else:

                    # Update counters according to the category of the read
                    if read.tid == -1:
                        count["unmapped"] += 1
                    elif len(read.query_alignment_sequence) < self.min_size:
                        count["short_mapped"] += 1
                    else: # not unmapped but mapq < min_mapq
                        count["lowmapq"] +=1

                    # Regenerate fastq
                    unmapped_fastq_h.write("{}\n{}\n+\n{}\n".format(read.qname, read.seq, read.qual))

        return count

        #Removing the original bam file which is no longer needed
        #remove(bam)


    def process_paired_end (self):
        pass
        # cutadapt

        # bwa alignment

        # bwa post processing of mapped reads

        # bwa post processing of unmapped reads

    #~~~~~~~PRIVATE METHODS~~~~~~~#

    def _dict_to_report(self, d, tab=""):
        """
        Recursive function to return a text report from nested dict or OrderedDict objects
        """
        report = ""
        for name, value in d.items():
            if type(value) == OrderedDict or type(value) == dict:
                report += "{}{}\n".format(tab, name)
                report += self._dict_to_report(value, tab=tab+"\t")
            else:
                report += "{}{}\t{}\n".format(tab, name, value)
        return report

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#   TOP LEVEL INSTRUCTIONS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

if __name__ == '__main__':

    fastq_sweeper = FastqSweeper.class_init()
    fastq_sweeper()
