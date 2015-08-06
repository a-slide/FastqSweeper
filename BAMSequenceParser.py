#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-

"""
@package    FastqSweeper
@brief      Helper class parsing bam aligned sequence line in pysam compatible objects
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
from collections import OrderedDict
from re import compile as recompile

# Third party package
import pysam

# Local packages

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class BAMSequenceParser (object):
    """Create a pysam compatible AlignedSegment from a raw bam line
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __init__ (self, header, skip_secondary=True):

        print("\nParsing sam sequences")

        # Store defined variables
        self.header = header
        self.skip_secondary = skip_secondary

        # counters
        self.count = {"total":0, "invalid":0, "primary":0, "secondary":0}

        # For cigar to tuple code conversion
        self.cigar_to_code = {'M':0, 'I':1, 'D':2, 'N':3, 'S':4, 'H':5, 'P':6, '=':7, 'X':8}
        self.cigar_regex = recompile("(\d+)([MIDNSHP=X])")

    def __str__ (self):
        return "".join(["{}\t{}\n".format(i,j) for i, j in self.count.items()])

    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def parse_line (self, line):

        self.count["total"] += 1
        # Try bloc to skip to the next line in case a non standard line is found
        try:
            # Split the line and verify the number of fields
            split_line = line.strip().split("\t")
            assert len (split_line) >=  11, "Invalid number of field in bam aligned sequence"

            # Init an AlignedSegment
            read = Read(
                qname = split_line[0],
                flag = int(split_line[1]),
                rname = self.header.rname_to_refid(split_line[2]),
                pos = int(split_line[3])-1,
                mapq = int(split_line[4]),
                cigar = self._cigarstring_to_tuple(split_line[5]),
                rnext = self.header.rname_to_refid(split_line[6]),
                pnext = int(split_line[7])-1,
                seq = split_line[9],
                qual = pysam.fromQualityString(split_line[10]),
                tags = self._parse_tags(split_line [11:]) if len (split_line) >= 12 else tuple())

            # skip the read if secondary and required
            if read.is_secondary:
                self.count["secondary"] += 1
                if self.skip_secondary:
                    return None

            # finally return the read
            self.count["primary"] += 1
            return read

        except Exception as E:
            print E
            print "Invalid sequence line in bam: {}".format(line)
            self.count["invalid"] += 1
            return None

    def _cigarstring_to_tuple (self, cigarstring):
        """Transform the cigar string in a pysam compatible tuple"""
        if cigarstring == "*":
            return None
        else:
            return tuple([(self.cigar_to_code[y], int(x)) for x, y in self.cigar_regex.findall(cigarstring)])

    def _parse_tags (self, tags_split):
        """Parse the optional tag fields"""
        tags = []

        for tag in tags_split:
            tag_split = tag.split(":")
            key = tag_split[0]
            value = int(tag_split[2]) if tag_split[1] == 'i' else tag_split[2]
            tags.append((key, value,))

        return tuple(tags)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class Read (object):
    """Wrapping class over a pysam AlignedSegment to surdefine, transform the init and
    add dedicated methods"""
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __init__ (self, qname, flag, rname, pos, mapq, cigar, rnext, pnext, seq, qual, tags):
        """Create a pysam.AlignedSegment object and define fields"""

        self.read = pysam.AlignedSegment()
        self.read.query_name = qname
        self.read.flag = flag
        self.read.reference_id = rname
        self.read.reference_start = pos
        self.read.mapping_quality = mapq
        self.read.cigar = cigar
        self.read.next_reference_id = rnext
        self.read.next_reference_start = pnext
        self.read.template_length = self.read.infer_query_length() if self.read.cigar else 0
        self.read.query_sequence = seq
        self.read.query_qualities = qual
        self.read.tags = tags
        self.is_secondary = self.read.is_secondary

    def __str__ (self):
        return str(self.read)

    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def is_properly_mapped(self, min_mapq, min_match_size):
        """Is properly mapped if mapped on e sequence with a mapq score higher the min_mapq
        and a template lenght higher than min_match_size"""
        return (self.read.tid != -1 and self.read.mapq >= min_mapq and self.read.template_length >= min_match_size)

    def to_bam (self):
        return self.read

    def to_fastq (self):
        return "@{}\n{}\n+\n{}\n".format(self.read.qname, self.read.seq, self.read.qual)
