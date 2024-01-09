import argparse
from Bio import AlignIO, SeqIO, Align
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Align import AlignInfo
import os

def parse_arguments():
    parser = argparse.ArgumentParser(description="Get the pan proteome of multiple protein assemblies in fasta format.")
    parser.add_argument("input_dir", help="Path to the directory containing input protein assembly FASTA files")
    parser.add_argument("output", help="Path to the output pan proteome FASTA file")
    return parser.parse_args()

def read_fasta(file_path):
    records = list(SeqIO.parse(file_path, "fasta"))
    return records

def process_batch(input_files, aligner, all_sequences, identity, coverage):
    batch_sequences = []
    print("Processing batch...")
    for input_file in input_files:
        sequences = read_fasta(input_file)
        batch_sequences += sequences
        all_sequences += sequences
        print(f"Read {len(sequences)} sequences from {input_file}")

    alignments = []
    for seq1 in batch_sequences:
        for seq2 in all_sequences:
            if seq1.id != seq2.id:
                alignment = aligner.align(seq1.seq, seq2.seq)
                if (alignment.score / min(len(seq1), len(seq2))) * 100 >= identity and \
                   (alignment.score / min(len(seq1), len(seq2))) >= coverage:
                    alignments.append(alignment)

    return alignments

def get_pan_proteome(input_dir, output_path, identity=95, coverage=0.9):
    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.match = 1
    aligner.mismatch = 0
    aligner.open_gap_score = -1
    aligner.extend_gap_score = -0.5
    aligner.substitution_matrix = Align.substitution_matrices.load("BLOSUM62")

    all_sequences = []
    all_alignments = []

    for root, dirs, files in os.walk(input_dir):
        for file in files:
            if file.endswith(".fasta"):
                file_path = os.path.join(root, file)
                alignments = process_batch([file_path], aligner, all_sequences, identity, coverage)
                all_alignments.extend(alignments)

    consensus_alignment = AlignInfo.SummaryInfo(Align.MultipleSeqAlignment(all_alignments)).dumb_consensus()
    consensus_sequence = Seq(str(consensus_alignment), consensus_alignment.alphabet)

    pan_proteome_record = SeqRecord(consensus_sequence, id="pan_proteome", description="Pan Proteome")
    SeqIO.write([pan_proteome_record], output_path, "fasta")

    print(f"Pan proteome written to {output_path}")

def main():
    args = parse_arguments()

    print("Starting pan proteome generation...")
    get_pan_proteome(args.input_dir, args.output)

if __name__ == "__main__":
    main()
