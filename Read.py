import pysam;
import sys
import numpy as np;
from pysam import TabixFile

class Read(object):

    def __init__(self, first_contig, second_contig, str):
        self.start_contig = first_contig;
        self.end_contig = second_contig;
        self.startone = int(str[1]);
        self.endone = int(str[2]);
        self.starttwo = int(str[4]);
        self.endtwo = int(str[5]);
        self.read_orientation = str[8] + str[9]
        self.str = str;

    def setStartContig(self, contig):
        self.start_contig = contig;

    def setEndContig(self, contig):
        self.end_contig = contig;

    def getStartContig(self):
        return self.start_contig;

    def getEndContig(self):
        return self.end_contig;

    def getStart1(self):
        return self.startone;

    def getEnd1(self):
        return self.endone;

    def getStart2(self):
        return self.starttwo;

    def getEnd2(self):
        return self.endtwo;

    def getOrientation(self):
        return self.read_orientation;

    def getStr(self):
        return self.str

    def changeOrientation(self, first_orientation, second_orientation):
        self.str[8], self.str[9] = first_orientation, second_orientation;
        self.read_orientation = first_orientation + second_orientation;

    def changeStart1Seq(self, str):
        self.str[1] = str;
        self.startone = int(str);

    def changeEnd1Seq(self, str):
        self.str[2] = str;
        self.endone = int(str);

    def changeStart2Seq(self, str):
        self.str[4] = str;
        self.starttwo = int(str);

    def changeEnd2Seq(self, str):
        self.str[5] = str;
        self.endtwo = int(str);

