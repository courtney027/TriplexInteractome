#!/usr/bin/python
"""
Main function: Search homopyrimidine:homopurine repeats with mirror symmetry from input fasta file
"""

import argparse
import time
from Bio import SeqIO
from multiprocessing import Pool

my_arp = argparse.ArgumentParser(description="Welcome to use it for searching H-DNA repeats")
group = my_arp.add_mutually_exclusive_group()
group.add_argument("-d", "--DNA", action="store_true")
my_arp.add_argument('SeqFile', help="the input DNA fasta file (Assign -d first!)")
my_arp.add_argument('OutFile', help="the output txt file containing repeat sequence information")
arp = my_arp.parse_args()

my_dict_DNA = SeqIO.index(arp.SeqFile, "fasta")

def mismatch(motif1, motif2):
    mis = sum(1 for c1, c2 in zip(motif1, motif2) if c1 != c2)
    error_rate = mis / len(motif1)
    return error_rate

def reverse(dna):
    rev = ''.join(reversed(dna))
    return rev

def analyze_sequence(sequence, ID):
    results = []
    error_rate = 0 # mistach rate (0~1)
    for seqStart, _ in enumerate(sequence[:-9], start = 1):
        for lenth in range(10, min(101, len(sequence) - seqStart + 1)): # repeats length (10~100)
            seqEnd = seqStart + lenth - 1
            motif = sequence[seqStart - 1:seqEnd]
            AG_count = motif.count('A', 0) + motif.count('G', 0)
            TC_count = motif.count('T', 0) + motif.count('C', 0)
            # define AG strand
            if motif.count('N', 0) == 0 and AG_count >= lenth - lenth * error_rate:
                for spacer, _ in enumerate(range(min(11, len(sequence) - 2 * lenth - seqStart + 2)), start=0):  # spacer length (0~10)
                    tfo_start = seqStart + lenth + spacer
                    tfo_end = seqStart + 2 * lenth + spacer - 1
                    tfo = sequence[tfo_start - 1:tfo_end]
                    triplex = sequence[seqStart - 1:tfo_end]
                    strand = "plus_AG"
                    if tfo.count('N', 0) == 0 and mismatch(reverse(tfo), motif) <= error_rate:
                        tri_type = "mirror"
                        results.append((ID, strand, tri_type, seqStart, tfo_end, lenth, spacer, triplex))
                    spacer += 1
            # define CT strand
            elif motif.count('N', 0) == 0 and TC_count >= lenth - lenth * error_rate:
                for spacer, _ in enumerate(range(min(11, len(sequence) - 2 * lenth - seqStart + 2)), start = 0):
                    tfo_start = seqStart + lenth + spacer
                    tfo_end = seqStart + 2 * lenth + spacer - 1
                    tfo = sequence[tfo_start - 1:tfo_end]
                    triplex = sequence[seqStart - 1:tfo_end]
                    strand = "plus_CT"
                    if tfo.count('N', 0) == 0 and mismatch(reverse(tfo), motif) <= error_rate:
                        tri_type = "mirror"
                        results.append((ID, strand, tri_type, seqStart, tfo_end, lenth, spacer, triplex))
                    spacer += 1
            lenth += 1
        seqStart += 1
    return results

def main():
    start_time = time.time()
    with open(arp.OutFile, 'w') as fw:
        p = Pool(processes = 8)
        if arp.DNA:
            results = [p.apply_async(analyze_sequence, args = (value.seq, ID)) for ID, value in my_dict_DNA.items()]
            print("#ID\tstrand\ttri_type\tseqStart\ttfo_end\tlenth\tspacer\ttriplex", file = fw)
        else:
            print("Please choose -d !")
            exit()
        p.close()
        p.join()
        for result in results:
            for data in result.get():
                print(*data, sep = "\t", file = fw)
    end_time = time.time()
    running_time = end_time - start_time
    # output running time
    print("Time of prediction：{:.5f} second".format(running_time))

if __name__ == '__main__':
    main()
