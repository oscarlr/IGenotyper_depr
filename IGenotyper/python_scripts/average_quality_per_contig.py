#!/bin/env python
import sys
import pysam

with pysam.FastxFile(sys.argv[1]) as fh:
    for entry in fh:
        qual_array = entry.get_quality_array()
        qual_sum = sum(qual_array)
        qual_avg = qual_sum/len(qual_array)
        length = len(qual_array)
        print "%s\t%s\t%s" % (entry.name,qual_avg,length)
