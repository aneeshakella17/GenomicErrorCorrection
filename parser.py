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

compliment = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N':'N'};
reverse_dictionary = {'+': '-', '-':'+'}

def get_tbx_array(str):
    case = {'t1': 'tests/t1/', 't2': 'tests/t2/', 't3': 'tests/t3/', 't4': 'tests/t4/', 't5': 'tests/t5/', 't6': 'tests/t6/', 't7': 'tests/t7/'};
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
    case = {'t1': 'tests/t1/', 't2': 'tests/t2/', 't3': 'tests/t3/', 't4': 'tests/t4/', 't5': 'tests/t5/', 't6': 'tests/t6/', 't7': 'tests/t7/'};

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
        firstSet = False;
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
                for edge in first_contig.edges:
                    if edge.getEndContig().getName() == contig.getName():
                        new_read = Read(first_contig, contig, new_row);
                        edge.addRead(new_read);
                        edge_found = True;
                        break;

                if not edge_found:
                    edge = Edge(first_contig, contig);
                    new_read = Read(first_contig, contig, new_row);
                    edge.addRead(new_read);
                    first_contig.add_edge(edge);

            contig = contig.next;

    calc_orientation(head);


def fix_orientation(head, fasta):
    contig = head;
    changes = [];

    while contig is not None and contig.next is not None:
        edge_orientation = {"++": 0, "--":0, "+-":0, "-+":0};
        edges = contig.edges;
        for edge in edges:
            edge_orientation = edge.getOrientation();

            rr = edge_orientation["--"];
            ff = edge_orientation["++"];
            fr = edge_orientation["+-"]

            if(rr > fr or ff > fr):


                if(rr > fr):
                    fasta = reverse_compliment(contig, fasta);
                elif(ff > fr):
                    contig = edge.sink;
                    fasta = reverse_compliment(contig, fasta);


                changes.append(contig.getName());


                for read in contig.reads:
                    read.changeOrientation(read.getOrientation()[1], read.getOrientation()[0]);

                    start= contig.getEnd() - read.getStart1();
                    start2 = contig.getEnd() - read.getStart2();
                    end = contig.getEnd() - read.getEnd1();
                    end2 = contig.getEnd() - read.getEnd2();


                    read.changeStart1Seq(contig.getStart() + start2);
                    read.changeEnd1Seq(contig.getStart() + end2);
                    read.changeStart2Seq(contig.getStart() + start)
                    read.changeEnd2Seq(contig.getStart() + end);

                #
                #
                # prev_contig = contig.prev;
                #
                # if(prev_contig is not None):
                #     for edge in prev_contig.edges:
                #         for read in edge.reads:
                #             if(read.getEndContig().getName() == contig.getName()):
                #                 if(read.getOrientation()[1] == "+"):
                #
                #                     read.changeOrientation(read.getOrientation()[0], "-");
                #                     start_distance = contig.getEnd() - read.getStart2();
                #                     end_distance = contig.getEnd() - read.getEnd2();
                #
                #                     read.changeStart2Seq(contig.getStart() + start_distance)
                #                     read.changeEnd2Seq(contig.getStart() + end_distance)
                #
                #
                # #


                for edge in contig.edges:
                    for read in edge.reads:
                        if(read.getOrientation()[0] == "-"):

                            read.changeOrientation("+", read.getOrientation()[1])
                            start_distance = contig.getEnd() - read.getStart1();
                            end_distance = contig.getEnd() - read.getEnd1();

                            read.changeStart1Seq(contig.getStart() + start_distance);
                            read.changeEnd1Seq(contig.getStart() + end_distance);

                contig = edge.source;
                break;

        contig = contig.next;

    return fasta, changes;



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
    reversed_seq = fasta[start:end + 1][::-1];
    compliment_seq = '';
    for i in range(0, len(reversed_seq)):
        char = reversed_seq[i];
        new_char = compliment[char]
        compliment_seq += new_char;

    fasta = fasta[0:start] + compliment_seq + fasta[end + 1:];
    return fasta;


def reset_orientation(head):
    contig = head;
    while contig is not None:
        contig.orientation['++'], contig.orientation['--'], contig.orientation['-+'], contig.orientation['+-'] = 0, 0, 0, 0;
        for edge in contig.edges:
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


        for edge in contig.edges:
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

def write_fasta(fasta, entry):
    f = open("output" + entry + ".fasta", "w");
    chr = "chr4";
    f.write(">" + chr + "\n");
    for i in range(0, len(fasta), 60):
        if(len(fasta) - i < 60):
            f.write(fasta[i:] + "\n");
        else:
            f.write(fasta[i:i+60] + "\n");


def crossover_detection(head):
    contig = head;
    potential_flips = [];
    greatest_edge = None;
    greatest_connect = float('-inf');
    #
    while contig is not None:
        for edge in contig.edges:
                if(edge.sink.getName() == contig.next.getName()):
                    break;
                else:
                    potential_flips.append([contig.next, edge.sink]);

        contig = contig.next;

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

    return potential_flips;


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
        # print(contig.getStart(), contig.getEnd());
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

    tmp = copy.copy(contig2);

    distance1 = contig1.getEnd() - contig1.getStart();
    distance2 = contig2.getEnd() - contig2.getStart();
    difference = contig2.getStart() - contig1.getEnd();
    contig2.setStart(contig1.getStart());
    contig2.setEnd(contig1.getStart() + distance2);
    contig1.setStart(contig2.getEnd() + difference);
    contig1.setEnd(contig2.getEnd() + difference + distance1);


    return new_fasta;


def main():
    string = "t4"
    tbx_array = get_tbx_array(string);
    count = 0;
    before_cross = False;
    changes_array = [];
    fasta_array = [];

    old_contig_parts = [];
    for entry in tbx_array:
        fast, test_unit = getFiles(string);
        new_fast = pre_process_fasta(fast);
        head_contig = Contig.create_contigs(test_unit);
        insert_all_reads(head_contig, entry, new_fast);
        fasta, changes = fix_orientation(head_contig, new_fast);
        contig_parts = crossover_detection(head_contig);
        print(len(contig_parts))
        if(len(contig_parts) > 0 or before_cross):
            if not before_cross or len(contig_parts) > len(old_contig_parts):
                for contig_part in contig_parts:
                    contig1, contig2 = contig_part[0], contig_part[1];
                    fasta = switch_contig(contig1, contig2, head_contig, fasta);
                    changes.append(contig1.getName() + ',' + contig2.getName());
                before_cross = True;
                old_contig_parts = contig_parts;
            else:
                for contig_part in old_contig_parts:
                    contig1, contig2 = contig_part[0], contig_part[1];
                    contig = head_contig;
                    while contig is not None:
                        if(contig1.getName() == contig.getName()):
                            contig1 = contig;
                        elif(contig2.getName() == contig.getName()):
                            contig2 = contig;
                        contig = contig.next;
                    fasta = switch_contig(contig1, contig2, head_contig, fasta);
                    changes.append(contig1.getName() + ',' + contig2.getName());


        write_fasta(fasta, str(count));
        count += 1;
        fasta_array.append(fasta)
        changes_array.append(changes)
        print(changes)

    best_index = changes_array.index(max(changes_array, key = lambda x : len(x)));

    for i in range(0, len(fasta_array)):
        write_fasta(fasta_array[best_index], str(i));




if __name__ == "__main__": main()





