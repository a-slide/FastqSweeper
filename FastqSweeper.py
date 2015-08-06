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
import shlex
from gzip import open as gopen
from collections import OrderedDict
from datetime import datetime

# Third party package
import pysam

# Local packages
from BAMHeader import BAMHeader
from BAMSequenceParser import BAMSequenceParser

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class FastqSweeper (object):
    """
    FastqSweeper use Cutadapt to trim reads based on quality, presence of N, and adapters and then
    perform a step of subtractive mapping with bwa Mem. BWA stream is parsed directly to separate
    the properly mapped reads that are written in a bam file, from the non-properly mapped reads
    that are written in a fastq file ready to be realign for further analysis.
    A report containing the options used and reads counts depending of their category is also
    generated.
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~CLASS FIELDS~~~~~~~#

    VERSION = "FastqSweeper 0.1"
    USAGE = "Usage: %prog -i INDEX -1 FASTQ_R1 [-2 FASTQ_R2] [...]"

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
            help= "Path to the bwa index basename")
        optparser.add_option('-1', '--fastq_R1', dest="fastq_R1",
            help= "Path to the fastq file")
        optparser.add_option('-2', '--fastq_R2', dest="fastq_R2",
            help= "[facultative] Path to the pair fastq file if paired end")
        optparser.add_option('-r', '--run', dest="run", action="store_true", default=False,
            help= "[facultative] Run command lines else the sotware run in demo mode (default: False)")
        optparser.add_option('-t', '--thread', dest="thread", default=1,
            help= "[facultative] Number of thread to use (default: 1)")
        optparser.add_option('-b', '--bwa_opt', dest="bwa_opt", default= "-M",
            help= "[facultative] bwa options for the mapping step (facultative and quoted) (default: -M)")
        optparser.add_option('-c', '--cutadapt_opt', dest="cutadapt_opt", default= "-m 25 -q 30,30 --trim-n",
            help= "[facultative] cutadapt options for the qc step (facultative and quoted) (default: -m 25 -q 30,30 --trim-n)")
        optparser.add_option('-a', '--adapter', dest="adapter",
            help= "[facultative] Path to a fasta file containing adapters to be 3' trimmed (facultative)")
        optparser.add_option('-m', '--min_mapq', dest="min_mapq", default= 0,
            help= "[facultative] Minimal mapq quality to be considered mapped (default: 0)")
        optparser.add_option('-s', '--min_match_size', dest="min_match_size", default= 0,
            help= "[facultative] Minimal size of match to be considered mapped (default: 0)")
        optparser.add_option('--skip_cutadapt', dest="skip_cutadapt", action="store_true", default=False,
            help= "[facultative] Skip the cutadapt step (default: False)")

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
            min_match_size=int(opt.min_match_size),
            skip_cutadapt=opt.skip_cutadapt)

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
        min_match_size=0,
        skip_cutadapt=False):
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
        self.min_match_size = int(min_match_size)
        self.skip_cutadapt=skip_cutadapt
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
            print ("\nRunning in demo mode")

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
                report.write("  Index basename:\t{}\n".format(self.index))
                report.write("  Fastq R1 path:\t{}\n".format(self.fastq_R1))
                report.write("  Fastq R2 path:\t{}\n".format(self.fastq_R2 if self.fastq_R2 else "None"))
                report.write("  Number of thread:\t{}\n".format(self.thread))
                report.write("  Bwa options:\t{}\n".format(self.bwa_opt))
                report.write("  Cutadapt options:\t{}\n".format(self.cutadapt_opt))
                report.write("  Adapter file:\t{}\n".format(self.adapter if self.adapter else "None"))
                report.write("  Cutoff MapQ score:\t{}\n".format(self.min_mapq))
                report.write("  Cutoff match length:\t{}\n".format(self.min_match_size))
                report.write("  Mode used:\t{}\n".format("single-end" if self.single else "paired-end"))
                report.write("  Output files basename:\t{}\n".format(self.basename))
                report.write("\nREAD COUNT PER CATEGORY\n")
                for key, value in count.items():
                    report.write("  {}:\t{}\n".format(key, value))

        # Finalize
        print ("\nDone in {}s".format(round(time()-start_time, 3)))
        return(0)

    def process_single_end (self):

        count = OrderedDict()

        ##### CUTADAPT #####
        if self.skip_cutadapt:
            print ("\nSkiping cutadapt step")
            trimmed_fastq=self.fastq_R1

        else:
            print ("\nStarting trimming with CUTADAPT")
            trimmed_fastq = self.basename+"_trim.fastq.gz"
            cutadapt_report = self.basename+"_trim_report.txt"

            cmd = "cutadapt {} {} {} -o {}".format (self.cutadapt_opt, \
            "-a file:"+self.adapter if self.adapter else "", self.fastq_R1, trimmed_fastq)
            print (cmd)

            if self.run:
                with open (cutadapt_report, "w") as fout:
                    for line in self.yield_cmd(cmd):
                        fout.write(line)

                # Extract values from cutadapt_report
                with open (cutadapt_report, "r") as fin:
                    for line in fin:
                        if line.startswith("Total reads processed:"):
                            count["Total reads before trimming"] = int(line.split()[-1].replace(",",""))
                        if line.startswith("Reads with adapters:"):
                            count["Reads with adapters"] = int(line.split()[-2].replace(",",""))
                        if line.startswith("Reads that were too short:"):
                            count["Reads that were too short"] = int(line.split()[-2].replace(",",""))
                        if line.startswith("Reads written (passing filters):"):
                            count["Reads after trimming"] = int(line.split()[-2].replace(",",""))

        ##### BWA #####
        print ("\nStart aligning with BWA MEM and sort reads")
        mapped_bam = self.basename+"_mapped.bam"
        unmapped_fastq = self.basename+"_clean.fastq.gz"

        # Prepare the command line
        cmd = "bwa mem {0} -t {1} {2} {3}".format(self.bwa_opt, self.thread, self.index, trimmed_fastq)
        print (cmd)

        if self.run:
            #counters
            total = mapped = unmapped = 0

            # Initialize the stream line per line generator
            sam = self.yield_cmd(cmd)

            # Initialize and parse the bam header
            h = BAMHeader ()
            for line in sam:

                # Add the line to BAMheader object until the first non header line is found
                if line.startswith("@"):
                    h.add_header_line(line)
                else:
                    break

            # Initialize a bam read parser
            bam_parser = BAMSequenceParser (header=h, skip_secondary=False)

            # Create an output bam file for mapped reads and an ouput fastq file for unmapped reads
            with \
                pysam.AlignmentFile (mapped_bam , "wb", header=h.header) as bam_out,\
                gopen (unmapped_fastq, "w") as fastq_out:

                # Process the first sequence found when parsing header
                read = bam_parser.parse_line(line)
                total += 1
                if read:
                    if read.is_properly_mapped(self.min_mapq, self.min_match_size):
                        bam_out.write(read.to_bam())
                        mapped += 1
                    else:
                        fastq_out.write(read.to_fastq())
                        unmapped += 1

                # Process the remaining sequences
                for line in sam:
                    read = bam_parser.parse_line(line)
                    total += 1
                    if total % 100000 == 0:
                        print ("{} sequence processed".format(total))
                    if read:
                        if read.is_properly_mapped(self.min_mapq, self.min_match_size):
                            bam_out.write(read.to_bam())
                            mapped += 1
                        else:
                            fastq_out.write(read.to_fastq())
                            unmapped += 1

            # Retrieve Count values from the BAMSequenceParser object
            count["Total reads processed by BWA"] = total
            count["Reads Mapped"] = mapped
            count["Reads Unmapped"] = unmapped
            count["Primary read"] = bam_parser.count["primary"]
            count["Secondary read"] = bam_parser.count["secondary"]
            count["Invalid reads"] = bam_parser.count["invalid"]

        return count

    def process_paired_end (self):
        pass

    #~~~~~~~PRIVATE METHODS~~~~~~~#

    def yield_cmd (self, cmd):
        """
        Decompose shell command in list of elementary elements and parse the output line by
        line using a yield statement
        """
        # Split the commands
        split_cmd = shlex.split(cmd)

        # Prepare the popen object
        proc = Popen(split_cmd, stdout=PIPE)

        # yield results line by line
        for line in proc.stdout:
            yield line
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#   TOP LEVEL INSTRUCTIONS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

if __name__ == '__main__':

    fastq_sweeper = FastqSweeper.class_init()
    fastq_sweeper()
