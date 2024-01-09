import argparse
from Bio import AlignIO, SeqIO, Align
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Align import AlignInfo

def parse_arguments():
    parser = argparse.ArgumentParser(description="Get the pan proteome of two protein assemblies in fasta format.")
    parser.add_argument("input1", help="Path to the first protein assembly FASTA file")
    parser.add_argument("input2", help="Path to the second protein assembly FASTA file")
    parser.add_argument("output", help="Path to the output pan proteome FASTA file")
    return parser.parse_args()

def read_fasta(file_path):
    records = list(SeqIO.parse(file_path, "fasta"))
    return records

def get_pan_proteome(sequences1, sequences2, output_path, identity=95, coverage=0.9):
    print("Reading sequences from input files...")
    combined_sequences = sequences1 + sequences2

    aligner = Align.PairwiseAligner()
    aligner.mode = 'global'
    aligner.match = 1
    aligner.mismatch = 0
    aligner.open_gap_score = -1
    aligner.extend_gap_score = -0.5
    aligner.substitution_matrix = Align.substitution_matrices.load("BLOSUM62")

    alignments = []
    print("Performing pairwise alignments...")
    for seq1 in combined_sequences:
        for seq2 in combined_sequences:
            if seq1.id != seq2.id:
                alignment = aligner.align(seq1.seq, seq2.seq)
                if (alignment.score / min(len(seq1), len(seq2))) * 100 >= identity and \
                   (alignment.score / min(len(seq1), len(seq2))) >= coverage:
                    alignments.append(alignment)

    consensus_alignment = AlignInfo.SummaryInfo(Align.MultipleSeqAlignment(alignments)).dumb_consensus()
    consensus_sequence = Seq(str(consensus_alignment), consensus_alignment.alphabet)

    pan_proteome_record = SeqRecord(consensus_sequence, id="pan_proteome", description="Pan Proteome")
    SeqIO.write([pan_proteome_record], output_path, "fasta")

    print(f"Pan proteome written to {output_path}")

def main():
    args = parse_arguments()

    print("Starting pan proteome generation...")
    sequences1 = read_fasta(args.input1)
    sequences2 = read_fasta(args.input2)

    get_pan_proteome(sequences1, sequences2, args.output)

if __name__ == "__main__":
    main()
