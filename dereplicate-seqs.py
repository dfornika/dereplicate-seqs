#!/usr/bin/env python

import sys
import argparse
import logging
from itertools import combinations

from Bio import SeqIO
from Bio.SeqIO import FastaIO
from Bio.SeqRecord import SeqRecord
from Bio.Data.CodonTable import TranslationError

def main(args):
    COUNT_SEPARATOR = '_x'
    seqs = {}

    for seq_record in SeqIO.parse(args.input_fasta, "fasta"):
        split_seqid = seq_record.id.split(COUNT_SEPARATOR)
        if len(split_seqid) == 1:
            seqs[split_seqid[0]] = {}
            seqs[split_seqid[0]]['seq'] = seq_record.seq
            count = 1
            seqs[split_seqid[0]]['count'] = count
        elif len(split_seqid) == 2:
            seqs[split_seqid[0]] = {}
            seqs[split_seqid[0]]['seq'] = seq_record.seq
            count = int(split_seqid[1])
            seqs[split_seqid[0]]['count'] = count
        else:
            logging.error("Error parsing: ", seq_record.id)

    # combinations('ABCD', 2) gives:
    # AB AC AD BC BD CD
    # ie. we don't need to compare AB and BA
    for seq1, seq2 in combinations(seqs, 2):
        # Need to skip over sequences that have been removed on previous iterations
        if seq1 not in seqs or seq2 not in seqs:
            continue

        # Translate each pair of seqs
        try:
            seq1_translated = seqs[seq1]['seq'].translate()
        except TranslationError:
            print("Error translating: " + seq1, file=sys.stderr)
        try:
            seq2_translated = seqs[seq2]['seq'].translate()
        except TranslationError:
            print("Error translating: " + seq2, file=sys.stderr)
    
        # Remove seq2 from collection of seqs if it translates to the
        # same amino acid sequence as seq1. Add the counts for seq2 to seq1 before removing.
        if seq1_translated == seq2_translated:
            print(seq1, " translates identicaly to ", seq2, ". Deleting ", seq2)
            seqs[seq1]['count'] += seqs[seq2]['count']
            del seqs[seq2]

    fasta_out = FastaIO.FastaWriter(open('output.fa', 'w'), wrap=None)
    fasta_out.write_file(
        (SeqRecord(seqs[seq]['seq'], id=seq + COUNT_SEPARATOR + str(seqs[seq]['count']), description="") for seq in seqs)
    )

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("input_fasta", help="Input fasta file to dereplicate")
    args = parser.parse_args()
    main(args)
