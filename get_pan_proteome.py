import argparse
from Bio import SeqIO

def get_pan_proteome(file1, file2):
    proteins1 = {}
    proteins2 = {}
    
    for record in SeqIO.parse(file1, "fasta"):
        proteins1[record.id] = record.seq
        
    for record in SeqIO.parse(file2, "fasta"):
        proteins2[record.id] = record.seq
        
    pan_proteins = set(proteins1.keys()) | set(proteins2.keys())
    
    return pan_proteins, proteins1, proteins2

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Get pan proteome of two assemblies')
    parser.add_argument('-f1', '--file1', required=True, help='First assembly fasta file')
    parser.add_argument('-f2', '--file2', required=True, help='Second assembly fasta file')
    parser.add_argument('-o', '--output', required=True, help='Output pan proteome fasta file')
    args = parser.parse_args()
    
    pan_proteins, proteins1, proteins2 = get_pan_proteome(args.file1, args.file2)
    
    with open(args.output, 'w') as f:
        for protein in pan_proteins:
            if protein in proteins1:
                seq = proteins1[protein]
            else:
                seq = proteins2[protein]
                
            f.write('>' + protein + '\n')
            f.write(str(seq) + '\n')