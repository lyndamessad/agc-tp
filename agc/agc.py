#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""

import argparse
import sys
import os
import gzip
import statistics
from collections import Counter
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw

__author__ = "Your Name"
__copyright__ = "Universite Paris Diderot"
__credits__ = ["Your Name"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Your Name"
__email__ = "your@email.fr"
__status__ = "Developpement"


def isfile(path):
    """Check if path is an existing file.
      :Parameters:
          path: Path to the file
    """
    if not os.path.isfile(path):
        if os.path.isdir(path):
            msg = "{0} is a directory".format(path)
        else:
            msg = "{0} does not exist.".format(path)
        raise argparse.ArgumentTypeError(msg)
    return path


def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True,
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default = 400,
                        help="Minimum sequence length for dereplication")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default = 10,
                        help="Minimum count for dereplication")
    parser.add_argument('-c', '-chunk_size', dest='chunk_size', type=int, default = 100,
                        help="Chunk size for dereplication")
    parser.add_argument('-k', '-kmer_size', dest='kmer_size', type=int, default = 8,
                        help="kmer size for dereplication")
    parser.add_argument('-o', '-output_file', dest='output_file', type=str,
                        default="OTU.fasta", help="Output file")
    return parser.parse_args()


def read_fasta(amplicon_file :str,minseqlen : int):
    """
    This function is used to read the amplican file.
    The file can be zipped or not.
    The first thing to do is to look at the file format, if it's a gzip file
    or normal fasta file.
    From this file, the function select sequences with minimum length >= the
    minseqlen gived as second parametre.

    Parametre:
    ----------
    amplicon_file: gzip fasta file or fasta file// contains sequences
    minseqlen: int// minimum sequence length (by default=400pb)

    Return:
    -------
    sequence yield generator.
    """
    #for gziped files:
    if amplicon_file.endswith(".gz"):
        with gzip.open(amplicon_file, "rb") as file:
            sequence = b""
            for line in file:
                if line.startswith(b">"):
                    if len(sequence) >= minseqlen:
                        yield sequence.decode("ascii")
                    sequence = b""
                else:
                    sequence += line.strip()
        yield sequence.decode("ascii")

    # for no gziped fasta files:
    else:
        with open(amplicon_file, "r") as file:
            sequence = ""
            for line in file:
                if line.startswith(">"):
                    if len(sequence) >= minseqlen:
                        yield sequence
                    sequence = ""
                else:
                    sequence += line.strip()
        yield sequence


def dereplication_fulllength(amplicon_file :str ,minseqlen : int,mincount :int):
    """
    This function is used to .
    First, it uses the read_fasta function to generate sequences generator.
    This generator will be used and modified to select sequences
    with appeared at least mincount time in the previous selected sequences.

    Parametre:
    ----------
    amplicon_file: fasta file// contains sequences
    minseqlen: int// minimum sequence length (by default=400pb)
    mincount: int// minimum recurence of the sequence in the fasta file.

    Return:
    -------
    sequence generator.yield [sequence, count]
    """
    sequences_all = []
    for seq in read_fasta(amplicon_file, minseqlen):
        sequences_all.append(seq)

    for i in Counter(sequences_all).most_common():
        if i[1] > mincount:
            yield i


def get_chunks(sequence, chunk_size:int):
    """
    This function is used to read the sequences file.
    From this file, the function select sequences with minimum length >= the
    minseqlen gived as second parametre.

    Parametre:
    ----------
    sequence: sequence generator.
    chunk_size: int// specify an exact size for chunk
    Return:
    -------
    selected_chunk: list //sous sequence taille I, 4 segment par sequence
    """
    # On verifie que la sequence peut etre découper en 4 segments.
    if len(sequence) < 4*chunk_size:
        raise ValueError

    selected_chunk = []
    for i in range(1,len(sequence)):
        if i*chunk_size<len(sequence):
            selected_chunk.append(sequence[i*chunk_size - chunk_size:i * chunk_size])
    return selected_chunk


def cut_kmer(sequence, kmer_size):
    """
    This function is used to build a kmer generator.
    Using specified kmer size, we cut every sequence to kmers with length equal
    to the kmer size gived as argument to the function.

    Parametre:
    ----------
    sequence: sequence generator.
    kmer_size: int// specify the exact kmer's size.

    Return:
    -------
    kmer generator // yield kmer
    """
    for i in range(len(sequence)-kmer_size+1):
        yield sequence[i:i+kmer_size]


def get_unique_kmer(kmer_dict:dict, sequence, id_seq:int, kmer_size:int):
    """
    This function is used to get kmer dictionnary.
    kmer_dict: {kmer, list(id_seq)}

    Parametre:
    ----------
    kmer_dict: dict // empty dictionnary kmers: {kmer, list(id_seq)}.
    sequence: sequence generator.
    kmer_size: int// longueur de kmer.
    id_seq: int // sequence id
    Return:
    -------
    kmer generator // yield kmer
    """
    # get kmers from the function cut_kmer()
    for kmer in cut_kmer(sequence,kmer_size):
        # verify that the  kmer is already in the dictionary, if not add it.
        if kmer not in kmer_dict:
            kmer_dict[kmer] = [id_seq]

        #verify that each kmer have all id of sequence where it's found.
        elif id_seq not in kmer_dict[kmer]:
            kmer_dict[kmer].append(id_seq)

    return kmer_dict


def search_mates(kmer_dict, sequence, kmer_size):
    """
    takes a  kmerd_dict.
    use Counter from the Collections library and its most_common function
    to identify the 8 most similar sequences to our entered sequence.

    Parametres:
    ----------
    kmer_dict : dictionary {'kmer':[id_seq]}
    sequence :  str//
    kmer_size : int
    """
    return [i[0] for i in Counter([ids for kmer in cut_kmer(sequence, kmer_size)
            if kmer in kmer_dict for ids in kmer_dict[kmer]]).most_common(8)]


def get_identity(alignment_list:list):
    """
    calculates the percentage of identity between the two sequences
    according to the formula:
    id = nb identical nucleotides/length of the alignment.

    Parametres:
    ----------
    alignment_list : list// contains 2 aligned sequences ["seq1", "seq2"]

    Return:
    ------
    id : float.
    """
    nb_nuc_tot = 0
    length_ali = len(alignment_list[0])

    for i in range(0,len(alignment_list)-1):
        for j in range(len(alignment_list[0])):
            if alignment_list[i][j] == alignment_list[i+1][j]:
                nb_nuc_tot =nb_nuc_tot+1

    return (nb_nuc_tot/length_ali)*100


def get_unique(ids):
    """
    function to get unique ids
    """
    return {}.fromkeys(ids).keys()

def common(lst1, lst2):
    """
    common function
    """
    return list(set(lst1) & set(lst2))

def detect_chimera(perc_identity_matrix):
    """
    This function is used to identify chimeric sequences.
    It is based on the calculation of the mean standard
    deviation of the identity percentages.
    If the result is greater than 5 and that at least 2 segments
    of our sequence show a similarity different to one of the two
    parents, we will identify that sequence as chimeric.

    Parametres
    ---------
    perc_identity_matrix: Matrix giving per segment the identity rate
            between the candidate sequence and two parent sequences.

    Return:
    -------
    boolean : - True -> sequence condidate is a chimera.
              - False -> sequence condidate is not a chimera.
    """
    stdevs = []  #liste for standard deviation.
    flag_file = 0
    flag_sim = 0

    for elmt in perc_identity_matrix:
        # standard deviation calcul
        stdevs.append(statistics.stdev([elmt[0], elmt[1]]))
        if flag_file == 0:
            val0 = elmt[0]
            val1 = elmt[1]
            flag_file = 1
        else:
            if flag_sim == 1:
                continue
            if val0 != elmt[0] and val1 != elmt[1]:
                flag_sim = 1
            val0 = elmt[0]
            val1 = elmt[1]
    std_mean = statistics.mean(stdevs)
    if std_mean > 5 and flag_sim == 1:
        return True
    return False


def chimera_removal(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    """
    chimera removal function
    """
    pass


def abundance_greedy_clustering(amplicon_file, minseqlen, mincount, chunk_size, kmer_size):
    """
    fonction abundance greedy clustering
    """
    pass


def write_OTU(OTU_list, output_file):
    """
    This function is used to write the result in output file.

    Parameter:
    ---------
    OTU_list: list of sequence and its occuracy.
    output_file: fasta file.
    Return:
    ------
    file : str // output file in fasta format generated by the functon fill()
    """
    with open(output_file, "w") as file:
        for i in range(0,len(OTU_list)):
            line = ">OTU_"+ str(i+1) + " occurrence:"+ str(OTU_list[i][1])+"\n"
            file.write(line)
            file.write(fill(str(OTU_list[i][0])))
            file.write("\n")

def fill(text, width=80):
    """Split text with a line return to respect fasta format"""
    return os.linesep.join(text[i:i+width] for i in range(0, len(text), width))


#==============================================================
# Main program
#==============================================================
def main():
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    sequences = dereplication_fulllength(args.amplicon_file, args.minseqlen, args.mincount)


if __name__ == '__main__':
    main()
