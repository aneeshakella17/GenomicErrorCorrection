import pysam;
import sys
import numpy as np;
from pysam import TabixFile


class Edge(object):

    def __init__(self, start_contig, end_contig):
        # self.name = ;
        self.source = start_contig;
        self.sink = end_contig;
        self.orientation = {"++": 0, "--":0, "+-":0, "-+":0};
        self.count  = 0;
        self.reads = [];

    def getOrientation(self):
        return self.orientation;

    def getEndContig(self):
        return self.sink;

    def getStartContig(self):
        return self.source;

    def getStartContig(self):
        return self.source;


    def getCount(self):
        return self.count;

    def getStartContig(self):
        return self.source;

    def addRead(self, read):
        self.reads.append(read);
        self.increment(read.getOrientation());

    def increment(self, direction):
        self.orientation[direction] += 1;
        self.count += 1;