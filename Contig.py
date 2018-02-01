import pysam;
import sys
import numpy as np;
from pysam import TabixFile
head = None;

class Contig(object):

    def __init__(self, name, start, end, left):
        self.name = name;
        self.start = start;
        self.end = end;
        self.reads = [];
        self.orientation = {"++": 0, "--":0, "+-":0, "-+":0};
        self.prev = left;
        self.next = None;
        self.mean = 0 ;
        self.std = 0;
        self.edges = [];

    def getName(self):
        return self.name;

    def add_to_read(self, read):
        self.reads.append(read);

    def add_edge(self, edge):
        self.edges.append(edge);

    def getStart(self):
        return self.start;

    def getEnd(self):
        return self.end;

    def setStart(self, start):
        self.start = start;

    def setEnd(self, end):
         self.end = end;

    def setMean(self, mean):
        self.mean = mean;

    def setSTD(self, std):
        self.std = std;

    def getMean(self):
            return self.mean;

    def getSTD(self):
            return self.std;

    def changeEverything(self, contig):
        self.name = contig.name;
        self.start = contig.start;
        self.end = contig.end;
        self.reads = contig.reads;
        self.orientation = contig.orientation;
        self.mean = contig.mean;
        self.std = contig.std;
        self.edges = contig.edges;




    def getDominantOrientation(self):
        best_key = '';
        max_num = float('-inf');
        for key in self.orientation:
            if(self.orientation[key] > max_num):
                max_num = self.orientation[key];
                best_key = key;
        return best_key;



def create_contigs(file):
    contigs = [];
    prev = None;
    head = None;
    head_set = False;
    for row in file.fetch():
        new_row = row.split('\t');
        first_orientation = int(new_row[1]);
        second_orientation = int(new_row[2]);
        name = new_row[3];
        new_contig = Contig(name, first_orientation, second_orientation, prev);
        new_contig.prev = prev;
        if(prev != None):
            prev.next = new_contig;
        if not head_set:
            head_set = True;
            head = new_contig
        prev = new_contig;

    return head;


def getHead():
    return head;