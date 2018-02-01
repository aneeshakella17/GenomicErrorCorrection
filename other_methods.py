import glob
import pysam;
import Contig
import Read
from Read import Read
from Edge import Edge
import copy;
import os
import numpy as np;
from pysam import TabixFile

def compression_check(head):
    contig = head;

    while contig.next is not None:
        lst = np.array([]);
        total_lst  = np.array([])
        for read in contig.reads:
            if(read.getStartContig().getName() == contig.getName() and read.getEndContig().getName() == contig.next.getName()):
                lst = np.append(lst, np.abs(read.getStart1() - read.getEnd2()));


        for read in contig.reads:
            if(read.getStartContig().getName() == contig.getName() and read.getEndContig().getName() == contig.getName() and read.getOrientation() == "+-"):
                total_lst = np.append(total_lst, np.abs(read.getStart1() - read.getEnd2()));


        for read in contig.next.reads:
            if(read.getStartContig().getName() == contig.next.getName() and read.getEndContig().getName() == contig.next.getName() and read.getOrientation() == "+-"):
                total_lst = np.append(total_lst, np.abs(read.getStart1() - read.getEnd2()));

        if(len(total_lst) != 0):
            mean = np.mean(total_lst)
            std = np.std(total_lst)
            new_mean = 0;
            if(len(lst) != 0):
                new_mean = np.mean(lst);
            print(len(lst))
            zscore = (new_mean -  mean)/(std);
            print("ZSCORE is", zscore)


        contig = contig.next;


def crossover_detection(head):
    contig = head;
    potential_flips = [];
    greatest_edge = None;
    greatest_connect = float('-inf');

    while contig is not None:
        greatest_connect = float('-inf');
        for edge in contig.edges:
            if(edge.getCount() > greatest_connect):
                greatest_connect = edge.getCount();
                greatest_edge = edge;

        if(contig.next is not None and greatest_edge is not None and
                   greatest_edge.getEndContig().getName() == contig.next.getName()):
            contig = contig.next;
            continue;


        next_edge = contig.next;

        greatest_edge2 = None;
        greatest_connect2 = float('-inf');

        while next_edge is not None:
            for edge in next_edge.edges:
                if(edge.getCount() > greatest_connect2):
                    greatest_connect2 = edge.getCount();
                    greatest_edge2 = edge;
            next_edge = next_edge.next;

        if(greatest_edge is not None and greatest_edge2 is not None ):
            potential_flips = [greatest_edge2.getStartContig(), greatest_edge.getEndContig()]

        contig = next_edge;

    if len(potential_flips) == 0:
        return potential_flips, False;
    else:
        return potential_flips, True;
