#!/usr/bin/env python2.7
# -*- coding: utf-8 -*-

"""
@package    FastqSweeper
@brief      Helper class parsing bam header line by line
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

# Third party package

# Local packages

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
class BAMHeader (object):
    """Create a pysam compatible sam header by parsing sam header lines
    """
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

    #~~~~~~~FONDAMENTAL METHODS~~~~~~~#

    def __init__ (self):

        print("\nParsing sam header")

        # Init the header dict with a minimal header line
        self.header = {'HD':{'VN': '1.5', 'SO':"unknown"}, 'SQ':[]}

        # Init enpty dict fo conversion between rname and refid
        self.name_to_id = OrderedDict()
        self.id_to_name = OrderedDict()
        self.current_id = 0

    #~~~~~~~PUBLIC METHODS~~~~~~~#

    def add_header_line (self, line):

        # Try except bloc to skip to the next line in case a non standard line is found
        try:
            split_line = line.strip().split("\t")
            tag = split_line[0].lstrip("@")

            # In case a HD tag is found will replace the minimal header
            if tag == "HD":
                for field in split_line[1:]:
                    key, value = field.split(":")
                    self.header["HD"][key] = value

            # In case a SQ tag is found = build a dict to have correspondance between id and refname
            elif tag == "SQ":
                d={}
                for field in split_line[1:]:
                    key, value = field.split(":")
                    d[key] = int(value) if value.isdigit() else value

                # append to the list of SQ lines
                self.header[tag].append(d)

                # append to the dict of correspondance
                self.name_to_id [d['SN']] = self.current_id
                self.id_to_name [self.current_id] = d['SN']
                self.current_id += 1

            # For all the other tags, create an entry with a list of values found
            else:
                d ={}
                for field in split_line[1:]:
                    key, value = field.split(":")
                    d[key] = int(value) if value.isdigit() else value
                if tag not in self.header:
                    self.header[tag]=[d]
                else:
                    self.header[tag].append(d)

        except ValueError as E:
            print "Invalid line in bam header : {}".format(line)


    def rname_to_refid (self, name):
        """return the id corresponding to the reference sequence name (pysam compatible)"""
        if name == "*":
            return -1
        else:
            return self.name_to_id[name]

    def refid_to_refname (self, id):
        """return the name corresponding to the reference sequence id (pysam compatible)"""
        if id == -1:
            return "*"
        else:
            return self.id_to_name[id]
