from Bio import pairwise2
from Bio.pairwise2 import format_alignment

def compare_dna(seq1, seq2):
    alignments = pairwise2.align.globalxx(seq1, seq2)
    best_alignment = alignments[0]
    aligned_seq1, aligned_seq2, score, start, end = best_alignment
    
    matches = sum(1 for a, b in zip(aligned_seq1, aligned_seq2) if a == b)
    total_length = max(len(seq1), len(seq2))
    percent_identity = (matches / total_length) * 100
    
    differences = []
    for i, (nuc1, nuc2) in enumerate(zip(aligned_seq1, aligned_seq2), start=1):
        if nuc1 != nuc2:
            if nuc1 == "-":
                diff_type = "Insertion"
            elif nuc2 == "-":
                diff_type = "Deletion"
            else:
                diff_type = "Mismatch"
            differences.append((i, nuc1, nuc2, diff_type))
    
    print("Alignment:")
    print(format_alignment(*best_alignment))
    print(f"Percent Identity: {percent_identity:.2f}%")
    
    print("\nDifferences:")
    print("Position | Seq1 | Seq2 | Type")
    print("--------------------------------")
    for pos, nuc1, nuc2, diff_type in differences:
        print(f"{pos:<8} | {nuc1:<4} | {nuc2:<4} | {diff_type}")
    
    return percent_identity, differences

if __name__ == "__main__":
    seq1 = input("Enter first DNA sequence: ").upper()
    seq2 = input("Enter second DNA sequence: ").upper()
    compare_dna(seq1, seq2)
