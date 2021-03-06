import glob
import pysam;
import Contig
import Read
from Read import Read
from Edge import Edge
import copy;
import os
import numpy as np;
import gevent;
from pysam import TabixFile
import operator

compliment = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N':'N'};
reverse_dictionary = {'+': '-', '-':'+'}
edge_dict = {};


def get_tbx_array(str):
    case = {'t1': 'tests/t1/', 't2': 'tests/t2/', 't3': 'tests/t3/', 't4': 'tests/t4/', 't5': 'tests/t5/', 't6': 'tests/t6/', 't7': 'tests/t7/', 't8': 'tests/t8/'};
    file_extension = case[str];
    fast = open(file_extension + 'test.fasta');
    test_unit = file_extension + 'test.unit.bed.gz';
    test_unit = pysam.TabixFile(test_unit)


    array = [];
    for file_name in os.listdir(file_extension):
        if file_name.endswith('.gz') and file_name != 'test.unit.bed.gz':
            array.append(file_name);


    tbx_array = [];
    for file_name in array:
        tbx_file = file_extension + file_name;
        tbx = pysam.TabixFile(tbx_file);
        tbx_array.append(tbx);

    return tbx_array;


def getFiles(str):
    case = {'t1': 'tests/t1/', 't2': 'tests/t2/', 't3': 'tests/t3/', 't4': 'tests/t4/', 't5': 'tests/t5/', 't6': 'tests/t6/', 't7': 'tests/t7/', 't8': 'tests/t8/'};

    file_extension = case[str];


    fast = open(file_extension + 'test.fasta');
    test_unit = file_extension + 'test.unit.bed.gz';
    test_unit = pysam.TabixFile(test_unit)
    return fast, test_unit;

def pre_process_fasta(fast):
    sequence = fast.read();
    sequence = sequence.strip(' ');
    for i in range(0, len(sequence)):
        if(sequence[i] == '\n'):
            break;
    ret = sequence[i + 1:].replace('\n','');
    return ret;

def insert_all_reads(head, tbx, fasta, epsilon = 1):
    contig = head;
    total_rows = [];


    for row in tbx.fetch():
            new_row = row.split('\t');

            quality = int(new_row[7]);
            if(epsilon > quality):
                continue;


            first_position = int(new_row[1]);
            second_position = int(new_row[2]);
            third_position = int(new_row[4]);
            fourth_position = int(new_row[5]);
            read_position = new_row[8] + new_row[9];

            if(first_position > third_position):

                new_row[1] = third_position
                new_row[2] = fourth_position
                new_row[4] = first_position
                new_row[5] = second_position

            total_rows.append(new_row);


    for new_row in total_rows:
        first_position = int(new_row[1]);
        second_position = int(new_row[2]);
        third_position = int(new_row[4]);
        fourth_position = int(new_row[5]);

        contig = head;
        first_contig = None;

        while contig is not None:
            if ((contig.getStart() < first_position and contig.getEnd() > second_position) and
                    (contig.getStart() < third_position and contig.getEnd() > fourth_position)):
                    new_read = Read(contig, contig, new_row);
                    contig.add_to_read(new_read);
            elif(contig.getStart() < first_position and contig.getEnd() > second_position):
                first_contig = contig;
            elif(contig.getStart() < third_position and contig.getEnd() > fourth_position and first_contig is not None):
                edge_found = False;

                for edge in first_contig.forward_edges:
                    if edge.getEndContig().getName() == contig.getName():
                        new_read = Read(first_contig, contig, new_row);
                        edge.addRead(new_read);
                        edge_found = True;
                        break;

                if not edge_found:
                    edge = Edge(first_contig, contig);
                    new_read = Read(first_contig, contig, new_row);
                    edge.addRead(new_read);
                    first_contig.add_forward_edge(edge);
                    contig.add_back_edge(edge);

            contig = contig.next;


    calc_orientation(head);
    filter(head);

def fix_orientation(contigs, head_contig, fasta):

    changes = [];

    # contig = head_contig;
    # while contig is not None:
    #     frac = sum_edges(contig);
    #     print(contig.getName(), frac);
    #
    #     for edge in contig.back_edges:
    #         print(edge.getStartContig().getName(), edge.getEndContig().getName(), edge.orientation);
    #
    #     for edge in contig.forward_edges:
    #         print(edge.getStartContig().getName(), edge.getEndContig().getName(), edge.orientation);
    #
    #     contig = contig.next;


    for contig in contigs:


        #Look at the forward edges
        for edge in contig.forward_edges:

            edge_orientation = edge.getOrientation();
            ff = edge_orientation["++"];
            fr = edge_orientation["+-"];

            if(ff > fr):
                endContig = edge.getEndContig();
                frac = sum_edges(endContig);
                if frac < 0:
                    fasta = reverse_compliment(endContig, fasta);
                    changes.append(endContig.getName());
                    position_swap(contigs, head_contig, endContig);

        #Look at the backward edges
        for edge in contig.back_edges:

            edge_orientation = edge.getOrientation();
            rr = edge_orientation["--"];
            fr = edge_orientation["+-"];

            if(rr > fr):
                startContig = edge.getStartContig();
                frac= sum_edges(startContig);

    # contig = head_contig;
    # while contig is not None:
    #     frac = sum_edges(contig);
    #     print(contig.getName(), frac);
    #
    #     for edge in contig.back_edges:
    #         print(edge.getStartContig().getName(), edge.getEndContig().getName(), edge.orientation);
    #
    #     for edge in contig.forward_edges:
    #         print(edge.getStartContig().getName(), edge.getEndContig().getName(), edge.orientation);
    #
    #     contig = contig.next;

    return fasta, changes;


def reformat_head(contig_array, head_contig, changes):

    for change in changes:
        contig = head_contig;

        while(contig is not None):
            if(contig.getName() == change):
                position_swap(contig_array, head_contig, contig);
                break;

            contig = contig.next;

    return head_contig;


def position_swap(contigs_array, head_contig, contig):


    for read in contig.reads:

        read.changeOrientation(read.getOrientation()[1], read.getOrientation()[0]);

        start = contig.getEnd() - read.getStart1();
        start2 = contig.getEnd() - read.getStart2();
        end = contig.getEnd() - read.getEnd1();
        end2 = contig.getEnd() - read.getEnd2();

        read.changeStart1Seq(contig.getStart() + start2);
        read.changeEnd1Seq(contig.getStart() + end2);
        read.changeStart2Seq(contig.getStart() + start)
        read.changeEnd2Seq(contig.getStart() + end);



    for edge in contig.back_edges:
        for read in edge.reads:
            if (read.getOrientation()[1] == "+"):
                    read.changeOrientation(read.getOrientation()[0], "-");
            # else:
                #read.changeOrientation(read.getOrientation()[0], "+");

            start_distance = contig.getEnd() - read.getStart2();
            end_distance = contig.getEnd() - read.getEnd2();
            read.changeStart2Seq(contig.getStart() + start_distance)
            read.changeEnd2Seq(contig.getStart() + end_distance)



    for edge in contig.forward_edges:
        for read in edge.reads:
            if (read.getOrientation()[0] == "-"):
                read.changeOrientation("+", read.getOrientation()[1])
            # else:
            #     read.changeOrientation("-", read.getOrientation()[1]);

            end_distance = contig.getEnd() - read.getStart1();
            start_distance = contig.getEnd() - read.getEnd1();
            read.changeStart1Seq(contig.getStart() + start_distance);
            read.changeEnd1Seq(contig.getStart() + end_distance);

    calc_orientation(head_contig);


def analyze_orientation(tbx):
    wrong_pairs = [];

    for row in tbx.fetch():
        new_row = row.split('\t');
        first_orientation = new_row[8];
        second_orientation = new_row[9];
        if(first_orientation == second_orientation):
            wrong_pairs.append(new_row);

    return wrong_pairs;

def switch_orientation(seq1, seq2):
    reversed_seq1 = seq1[::-1];
    compliment_seq = reversed_seq1;
    for i in range(0, len(reversed_seq1)):
        char = reversed_seq1[i];
        new_char = compliment[char]
        compliment_seq = compliment_seq[:i] + new_char + compliment_seq[i + 1:]

    reversed_seq2 = seq2[::-1];
    compliment_seq2 = reversed_seq2;
    for i in range(0, len(reversed_seq1)):
        char = reversed_seq2[i];
        new_char = compliment[char]
        compliment_seq2 = compliment_seq2[:i] + new_char + compliment_seq2[i + 1:]

    return seq1, compliment_seq2;

def reverse_compliment(contig, fasta):
    start = contig.start;
    end = contig.end;
    reversed_seq = fasta[start:end][::-1];
    compliment_seq = '';

    for i in range(0, len(reversed_seq)):
        char = reversed_seq[i];
        new_char = compliment[char]
        compliment_seq += new_char;

    fasta = fasta[0:start] + compliment_seq + fasta[end:];

    return fasta;


def reset_orientation(head):
    contig = head;
    while contig is not None:
        contig.orientation['++'], contig.orientation['--'], contig.orientation['-+'], contig.orientation['+-'] = 0, 0, 0, 0;
        for edge in contig.forward_edges:
            edge.orientation['++'], edge.orientation['--'], edge.orientation['-+'], edge.orientation[
                '+-'] = 0, 0, 0, 0;
        contig = contig.next;


def calc_orientation(head):
    contig = head;
    reset_orientation(contig);
    lst = np.array([]);

    while contig is not None:

        for read in contig.reads:
            first_position = read.getStart1();
            second_position = read.getEnd1();
            third_position = read.getStart2();
            fourth_position = read.getEnd2();

            if ((contig.getStart() < first_position and contig.getEnd() > second_position) or
                (contig.getStart() < third_position and contig.getEnd() > fourth_position)):
                new_position = read.getOrientation();
                contig.orientation[new_position] += 1;


        for edge in contig.forward_edges:
            for read in edge.reads:
                edge.increment(read.getOrientation());


        contig = contig.next;


def printOrient(head):
    calc_orientation(head);
    contig = head;

    while contig is not None:
        print(contig.getName())
        print(contig.orientation);
        print(' ')
        for edge in contig.edges:
            print(edge.getStartContig().getName(), edge.getEndContig().getName());
            print(edge.getOrientation());
            print(' ');

        contig = contig.next;

def write_fasta(fasta):
    f = open("output.fasta", "w");
    chr = "chr4";
    f.write(">" + chr + "\n");
    for i in range(0, len(fasta), 60):
        if(len(fasta) - i < 60):
            f.write(fasta[i:] + "\n");
        else:
            f.write(fasta[i:i+60] + "\n");


def crossover_detection(contigs):
    potential_flips = [];

    for contig in contigs:
        found = False;

        for edge in contig.edges:
            if(edge.sink.getName() == contig.next.getName()):
                found = True;

        if found:
            continue;

        for edge in contig.edges:
            if(edge.sink.getName() != contig.next.getName() and not found):
                potential_flips.append([contig.next, edge.sink]);

    print(potential_flips)


    actual_flips = [];

    for flip in potential_flips:
        first_contig = flip[0]
        second_contig = flip[1];
        print(first_contig.getName(), second_contig.getName())
        for edge in first_contig.edges:
            if(edge.getEndContig().getName() == second_contig.getName()):
                actual_flips.append([first_contig, second_contig]);


    return actual_flips;


    # for flip in potential_flips:
    #     contig = head;
    #     first_contig = flip[0];
    #     second_contig = flip[1];
    #     print(first_contig.getName());
    #     print(second_contig.getName());
    #     print(' ')
    #
    #     while contig is not None:
    #         for edge in contig.edges:
    #             if (first_contig.getName() == edge.getStartContig().getName()):
    #                 lst = [];
    #                 lst2 = [];
    #                 print(first_contig.getName(), edge.sink.getName())
    #                 print(first_contig.getStart())
    #                 print(first_contig.getEnd())
    #                 print('BELOW IS MEAN')
    #                 print(' ')
    #                 for read in edge.reads:
    #                     lst.append(read.getStart1());
    #                     lst2.append(read.getEnd1());
    #                 print(np.mean(lst))
    #                 print(np.mean(lst2))
    #                 print(' ')
    #         contig = contig.next;

    # #

    # while contig is not None:
    #     for edge in contig.edges:
    #         lst = [];
    #         for read in edge.reads:
    #             lst.append(abs(read.getStart1() - read.getEnd2()));
    #         edge_mean = np.mean(lst)
    #         print(edge.source.getName(), edge.sink.getName(), edge_mean);
    #     contig = contig.next;

    #
    # lst1 = [];
    # lst2 = [];
    #
    # total_lst = [];
    #

    #
    # for flip in potential_flips:
    #     contig1, contig2 = flip[0], flip[1]
    #     main_contig = contig1.prev;
    #
    #     for edge in main_contig.edges:
    #         first_contig = edge.source
    #         second_contig = edge.sink
    #
    #         if(edge.sink.getName() == contig1.getName()):
    #             for read in edge.reads:
    #                 lst1.append(10000 - (abs(first_contig.getEnd() - read.getEnd1())) - (abs(second_contig.getStart() - read.getEnd2())));
    #
    #         if(edge.sink.getName() == contig2.getName()):
    #             for read in edge.reads:
    #                 lst2.append(10000 - (abs(first_contig.getEnd() - read.getEnd1())) - (abs(second_contig.getStart() - read.getEnd2())));
    #

    # contig = head
    # while contig is not None:
    #     for read in contig.reads:
    #         total_lst.append(abs(read.getStart1() - read.getEnd2()));
    #     contig = contig.next;
    #
    # lib_mean = np.mean(total_lst);
    # lib_std = np.std(total_lst);
    #
    #
    # for flip in potential_flips:
    #     contig1, contig2 = flip[0], flip[1]
    #     main_contig = contig1.prev;
    #
    #     for edge in main_contig.edges:
    #         first_contig = edge.source
    #         second_contig = edge.sink
    #
    #         if(edge.sink.getName() == contig1.getName()):
    #             for read in edge.reads:
    #                 lst1.append(lib_mean - (abs(first_contig.getEnd() - read.getStart1())) - (abs(second_contig.getStart() - read.getEnd2())));
    #
    #         if(edge.sink.getName() == contig2.getName()):
    #             for read in edge.reads:
    #                 lst2.append(abs(read.getEnd1() - read.getStart2()));
    #
    #
    #
    #     mean1 = 0;
    #     mean2 = 0;
    #     if(len(lst1) != 0):
    #         mean1 = np.mean(lst1);
    #     if(len(lst2) != 0):
    #         mean2 = np.mean(lst2)
    #
    #     z1 = (mean1 - lib_mean)/(lib_std);
    #     z2 = (mean2 - lib_mean)/((lib_std));
    #     print(z1, z2)

    # return potential_flips;

def switch_contig(contig1, contig2, head, fasta):


    contig = head;
    new_fasta = "";


    if(contig1.prev != None):
        contig1.prev.next = contig2;
    contig2.prev = contig1.prev;

    contig1.next = contig2.next;

    if(contig2.next != None):
        contig2.next.prev = contig1;

    contig1.prev = contig2;
    contig2.next = contig1;


    while(contig is not None):
        new_fasta += fasta[contig.getStart():contig.getEnd()];
        contig = contig.next;


    for read in contig1.reads:
        read.changeStart1Seq(read.getStart1() - (contig1.getEnd() - contig2.getEnd()))
        read.changeEnd1Seq(read.getEnd1() - (contig1.getEnd() - contig2.getEnd()) );
        read.changeStart2Seq(read.getStart2() - (contig1.getEnd() - contig2.getEnd()))
        read.changeEnd2Seq(read.getEnd2() - (contig1.getEnd() - contig2.getEnd()));

    contig = head;
    while contig is not None:
        for edge in contig.edges:
            if(edge.getEndContig().getName() == contig1.getName()):
                for read in edge.reads:
                    read.changeStart2Seq(read.getStart2() + (contig1.getEnd() - contig2.getEnd()));
                    read.changeEnd2Seq(read.getEnd2() + (contig1.getEnd() - contig2.getEnd()));
        contig = contig.next;

    for edge in contig1.edges:
        for read in edge.reads:
            read.changeStart1Seq(read.getStart1() + (contig1.getEnd() - contig2.getEnd()))
            read.changeEnd1Seq(read.getEnd1() + (contig1.getEnd() - contig2.getEnd()) );


    for read in contig2.reads:
        read.changeStart1Seq(read.getStart1() + (contig1.getStart() - contig2.getStart()))
        read.changeEnd1Seq(read.getEnd1() + (contig1.getStart() - contig2.getStart()) );
        read.changeStart2Seq(read.getStart2() + (contig1.getStart() - contig2.getStart()))
        read.changeEnd2Seq(read.getEnd2() + (contig1.getStart() - contig2.getStart()));

    contig = head;
    while contig is not None:
        for edge in contig2.edges:
            if(edge.getEndContig().getName() == contig1.getName()):
                for read in edge.reads:
                    read.changeStart2Seq(read.getStart2() + (contig1.getStart() - contig2.getStart()));
                    read.changeEnd2Seq(read.getEnd2() +(contig1.getStart() - contig2.getStart()));
        contig = contig.next;

    for edge in contig2.edges:
        for read in edge.reads:
            read.changeStart1Seq(read.getStart1() + (contig1.getStart() - contig2.getStart()))
            read.changeEnd1Seq(read.getEnd1() + (contig1.getStart() - contig2.getStart()) );

    distance1 = contig1.getEnd() - contig1.getStart();
    distance2 = contig2.getEnd() - contig2.getStart();
    difference = contig2.getStart() - contig1.getEnd();
    contig2.setStart(contig1.getStart());
    contig2.setEnd(contig1.getStart() + distance2);
    contig1.setStart(contig2.getEnd() + difference);
    contig1.setEnd(contig2.getEnd() + difference + distance1);


    return new_fasta;

def create_spacing(heads, fasta):

    means = [];

    for head in heads:
        lst = [];
        contig = head;
        while contig is not None:
            for read in contig.reads:
                 lst.append(abs(read.getStart1() - read.getEnd2()));
            contig = contig.next;
        uL = np.mean(lst);
        means.append(uL);


    spacings = [];

    while heads.count(None) is not len(heads):
         g_lst = [];
         for i in range(0, len(heads)):
            contig = heads[i];
            for edge in contig.edges:
               if(edge.getEndContig().getName() == contig.next.getName()):
                   first_contig = edge.getStartContig();
                   second_contig = edge.getEndContig();
                   for read in edge.reads:
                       first_dist = abs(read.getStart1() - first_contig.getEnd());
                       second_dist = abs(second_contig.getStart() - read.getEnd2())
                       g = means[i] - first_dist - second_dist;
                       g_lst.append(g);

            heads[i] = heads[i].next;

         if(len(g_lst) != 0):
             spacings.append(max(np.mean(g_lst), 100));
         else:
             spacings.append(100);

    spacings = [int(spacing) for spacing in spacings];
    print(spacings)

    return fasta;


def findDominantOrientation(edge):
    return max(edge.orientation.items(), key=operator.itemgetter(1))[0]

def mode(arr):
    dict = {"++":0 , "+-":0, "--":0, "-+": 0};
    if(len(arr) == 0):
        return None;
    for elem in arr: dict[elem] += 1;
    return max(dict.items(), key=operator.itemgetter(1))[0]

def filter(head_contig, default = 0):
    contig = head_contig;
    while contig is not None:
        contig.forward_edges = [edge for edge in contig.forward_edges if (edge.orientation["++"] + edge.orientation["+-"] + edge.orientation["--"] + edge.orientation["-+"]) > default];
        contig.back_edges = [edge for edge in contig.back_edges if (edge.orientation["++"] + edge.orientation["+-"] + edge.orientation["--"] + edge.orientation["-+"]) > default];



        contig = contig.next;

def sum_edges(contig):
    fr = 0;
    ff = 0;
    rr = 0;

    new_ff = 0;
    new_rr = 0;
    new_fr = 0;

    for edge in contig.back_edges:
        fr += edge.orientation["+-"];
        ff += edge.orientation["++"];

        new_ff += edge.orientation["++"]
        new_rr += edge.orientation["--"]
        new_fr += edge.orientation["+-"]


    for edge in contig.forward_edges:
        fr += edge.orientation["+-"];
        rr += edge.orientation["--"];
        new_ff += edge.orientation["++"]
        new_rr += edge.orientation["--"]
        new_fr += edge.orientation["+-"]


    return fr - ff - rr;


def check_flip(contig, fasta):
    fr1 = 0;
    fr2 = 0;

    for edge in contig.back_edges:
        fr1 += edge.orientation["+-"]

    for edge in contig.forward_edges:
        fr1 += edge.orientation["+-"]

    reverse_compliment(contig, fasta);

    for edge in contig.back_edges:
        fr2 += edge.orientation["+-"]

    for edge in contig.forward_edges:
        fr2 += edge.orientation["+-"]

    reverse_compliment(contig, fasta);

    if (fr2 > fr1):
        return True;
    else:
        return False;


def main():
    string = "t8"
    tbx_array = get_tbx_array(string);
    fast, test_unit = getFiles(string);
    fasta = pre_process_fasta(fast);
    changes = [];
    heads = [];

    for entry in tbx_array:
        fast, test_unit = getFiles(string);
        head_contig = Contig.create_contigs(test_unit);
        contig_array = Contig.create_ordered_contigs(head_contig);
        insert_all_reads(head_contig, entry, fasta);
        reformat_head(contig_array, head_contig, changes);
        fasta, temp = fix_orientation(contig_array, head_contig, fasta);
        changes.extend(temp)
        print(changes)
        # contig_parts = crossover_detection(contig_array)
        #
        # if(len(contig_parts) > 0):
        #
        #     for contig_part in contig_parts:
        #
        #         contig1, contig2 = contig_part[0], contig_part[1];
        #         fasta = switch_contig(contig1, contig2, head_contig, fasta);
        #         changes.append(contig1.getName() + ',' + contig2.getName());
        #         print(contig1.getName() + ',' + contig2.getName());
        #
        #
        # heads.append(copy.copy(head_contig));
        #
        # print(changes)
        break



    fasta = create_spacing(heads, fasta);
    write_fasta(fasta);
    os.system("nucmer --prefix check truth5.fasta output.fasta");
    os.system("delta-filter -1 check.delta >check.1.delta")
    os.system("mummerplot --png check.1.delta")







if __name__ == "__main__": main()





