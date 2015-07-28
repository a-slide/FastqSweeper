#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-

"""
@package    RefMasker
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
import optparse
import sys
from time import time

# Third party package
import pysam

# Local packages


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class Fastq_cleaner (object):
    """
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~CLASS FIELDS~~~~~~~#

    VERSION = "Fastq_cleaner 0.1"
    USAGE = "Usage: %prog -i INDEX -r FASTQ_R1 [-s FASTQ_R2 -t THREAD -b BWA_OPTIONS -c CUTADAPT_OPTIONS -a ADAPTER]"

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
        optparser.add_option('-r', '--fastq_R1', dest="fastq_R1", help= "Path to the fastq file (required)")
        optparser.add_option('-s', '--fastq_R2', dest="fastq_R2", help= "Path to the pair fastq file if paired end (facultative)")
        optparser.add_option('-t', '--thread', dest="thread", default=1, help= "Number of thread to use (default 1)")
        optparser.add_option('-b', '--bwa_opt', dest="bwa_opt", help= "bwa options for the mapping step (facultative and quoted)")
        optparser.add_option('-c', '--cutadapt_opt', dest="cutadapt_opt", default= "-m 25 -q 30,30 --trim-n",
            help= "cutadapt options for the qc step (facultative and quoted)")
        optparser.add_option('-a', '--adapter', dest="adapter", help= "Path to a fasta file containing adapters to be 3' trimmed (facultative)")
        ### Parse arguments
        opt, args = optparser.parse_args()

        ### Init a RefMasker object
        return Fastq_cleaner (opt.index, opt.fastq_R1, opt.fastq_R2, int(opt.thread), opt.bwa_opt, opt.cutadapt_opt, opt.adapter)

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __init__(self, index=None, fastq_R1=None, fastq_R2=None, thread=1, bwa_opt=None, cutadapt_opt=None, adapter=None):
        """
        General initialization function for import and command line
        """

        print("\nInitialize Fastq_cleaner")

        ## Verifications
        assert index, "A path to the bwa index is mandatory"
        assert fastq_R1, "A path to the fastq file is mandatory"

        self.index = index
        self.fastq_R1 = fastq_R1
        self.fastq_R2 = fastq_R2
        self.thread = int(thread)
        self.bwa_opt = bwa_opt
        self.cutadapt_opt = cutadapt_opt
        self.adapter = adapter
        self.mode = "PE" if fastq_R2 else "SE"

    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def __call__(self):
        """
        """

        start_time = time()

        if self.mode == "SE":
            self.process_single_end()
        else:
            self.process_paired_end()

        # Generate Reports

        # Finalize
        print ("\nDone in {}s".format(round(time()-start_time, 3)))
        return(0)


    def process_single_end (self):

        print ("\nProcessing fastq in single end mode")

        cutadapt -m 25 -q 30,30 --trim-n -a file:adapter.fa $i -o `basename ${i%.fastq.gz}_trim.fastq.gz` > `basename ${i%.fastq.gz}_report.txt`

        # bwa alignment

        # bwa post processing of mapped reads

        # bwa post processing of unmapped reads


    def process_paired_end (self):

        print ("\nProcessing fastq in single end mode")

        # cutadapt

        # bwa alignment

        # bwa post processing of mapped reads

        # bwa post processing of unmapped reads


#~~~~~~~COMMAND LINE UTILITIES~~~~~~~#


def run_command(cmd, stdin=None, ret_stderr=False, ret_stdout=False):
    """
    Run a command line in the default shell and return the standard output
    @param  cmd A command line string formated as a string
    @param  stdinput    Facultative parameters to redirect an object to the standard input
    @param  ret_stderr  If True the standard error output will be returned
    @param  ret_stdout  If True the standard output will be returned
    @note If ret_stderr and ret_stdout are True a tuple will be returned and if both are False
    None will be returned
    @return If no standard error return the standard output as a string
    @exception  OSError Raise if a message is return on the standard error output
    @exception  (ValueError,OSError) May be raise by Popen
    """
    # Function specific imports
    from subprocess import Popen, PIPE

    # Execute the command line in the default shell
    if stdin:
        proc = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=PIPE)
        stdout, stderr = proc.communicate(input=stdin)
    else:
        proc = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
        stdout, stderr = proc.communicate()

    if proc.returncode == 1:
        msg = "An error occured during execution of following command :\n"
        msg += "COMMAND : {}\n".format(cmd)
        msg += "STDERR : {}\n".format(stderr)
        raise Exception (msg)

    # Else return data according to user choices is returned
    if ret_stdout and ret_stderr:
        return stdout, stderr
    elif ret_stdout:
        return stdout
    elif ret_stderr:
        return stderr
    else:
        return None

def make_cmd_str(prog_name, opt_dict={}, opt_list=[]):
    """
    Create a Unix like command line string from a
    @param prog_name Name (if added to the system path) or path of the programm
    @param opt_dict Dictionnary of option arguments such as "-t 5". The option flag have to
    be the key (without "-") and the the option value in the dictionnary value. If no value is
    requested after the option flag "None" had to be asigned to the value field.
    @param opt_list List of simple command line arguments
    @exemple make_cmd_str("bwa", {"b":None, t":6, "i":"../idx/seq.fa"}, ["../read1", "../read2"])
    """

    # Start the string by the name of the program
    cmd = "{} ".format(prog_name)

    # Add options arguments from opt_dict
    if opt_dict:
        for key, value in opt_dict.items():
            if value:
                cmd += "-{} {} ".format(key, value)
            else:
                cmd += "-{} ".format(key)

    # Add arguments from opt_list
    if opt_list:
        for value in opt_list:
            cmd += "{} ".format(value)

    return cmd


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#   TOP LEVEL INSTRUCTIONS
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

if __name__ == '__main__':

    fastq_cleaner = Fastq_cleaner.class_init()
    fastq_cleaner()
