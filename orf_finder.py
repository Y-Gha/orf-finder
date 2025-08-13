"""
ORF Finder and Sequence Analyzer

This script allows users to input a DNA or RNA sequence and performs the following tasks:
- Detects whether the input is DNA or RNA
- Transcribes DNA to RNA if necessary
- Finds Open Reading Frames (ORFs) in all 6 reading frames (3 forward, 3 reverse)
- Translates RNA sequences to protein
- Calculates GC content
- Analyzes codon usage frequency
- Optionally saves protein sequences from ORFs to a text file

Intended for use in basic bioinformatics and computational biology analysis.

"""

from collections import Counter

# -----------------------------
# Genetic code for translation
# -----------------------------

GENETIC_CODE = {
    'AUG': 'M', 'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L',
    'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S', 
    'UAU': 'Y', 'UAC': 'Y', 'UAA': '*', 'UAG': '*', 'UGA': '*',
    'UGU': 'C', 'UGC': 'C', 'UGG': 'W',
    'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
    'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 
    'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
    'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

# ------------------------------
# All 64 possible codons for CODON USAGE INITIALIZATION
# ------------------------------
BASES = ['U', 'A', 'C', 'G']
ALL_CODONS = [a + b + c for a in BASES for b in BASES for c in BASES]

# ------------------------------
# Function 1: Identify seq as DNA or RNA
# ------------------------------
def sequence_type (seq):
    seq = seq.upper()
    if 'U' in seq and 'T' not in seq:
        return "RNA"
    elif 'T' in seq and 'U' not in seq:
        return "DNA"
    else:
        return "Unknown or mixed invalid input!"
    
# ------------------------------
# Function 2: Transcription
# ------------------------------
def dna_transcription (dna):
    return dna.upper().replace('T','U')

# ------------------------------
# Function 3: Reverse Complement
# ------------------------------
def reverse_complement (seq):
    seq = seq.upper()
    if 'T' in seq:
        complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
    elif 'U' in seq:
        complement = {'U':'A', 'A':'U', 'C':'G', 'G':'C'}
    else:
        return "Invalid sequence!"
    reserved_seq = seq[::-1]
    return ''.join(complement.get(base, '?') for base in reserved_seq)

# ------------------------------
# Function 4: GC Content
# ------------------------------
def gc_content (seq):
    seq = seq.upper()
    gc = seq.count('G') + seq.count('C')
    return round((gc /len(seq)*100),2) if len(seq) > 0 else 0 

# ------------------------------
# Function 5: Translation
# ------------------------------
def translation (rna):
    rna = rna.upper()
    protein = []
    for i in range(0,len(rna)-2,3):
        codon = rna[i:i+3]
        aminoacid = GENETIC_CODE.get(codon,'?')
        protein.append(aminoacid)
    return ''.join(protein)
    
# ------------------------------
# Function 6: Codon Usage Frequency (all codons included)
# ------------------------------

def codon_usage (rna):
    rna = rna.upper()
    codons = [rna[i:i+3] for i in range(0,len(rna)-2,3)]
    freq = Counter(codons)
    total = sum(freq.values())
    usage = {}
    for codon in ALL_CODONS:
        usage[codon] = round((freq.get(codon,0)/total)*100,2) if total > 0 else 0
    return usage
     
# ------------------------------
# Function 7: ORF Finder in one frame
# ------------------------------
def find_orfs_in_frame(rna,offset = 0, strand = '+',frame = 0, allow_overlap=False, total_length=None):
    rna = rna.upper()
    start_codon = 'AUG'
    stop_codons = {'UAG','UGA','UAA'}
    orfs = []
    i = 0
    while i < len(rna)-2:
        codon = rna[i:i+3]
        if codon == start_codon:
            for j in range(i+3,len(rna)-2,3):
                next_codon = rna [j:j+3]
                if next_codon in stop_codons:
                    orf_seq = rna [i:j+3]
                    protein = translation(orf_seq)
                    start_pos = i + offset if strand == '+' else total_length -(j+3+offset)
                    stop_pos = j + offset + 3 if strand == '+' else total_length - (i+offset)
                    orfs.append({
                        'start': start_pos,
                        'end': stop_pos,
                        'length': j + 3 - i,
                        'protein': protein,
                        'strand': strand,
                        'frame': frame 
                    })
                    i = i + 1 if allow_overlap else j + 3 
                    break
            else:
                i+=3 
        else: 
            i+=3
    return orfs

# ------------------------------
# Function 8: ORF Finder in all 6 frames
# ------------------------------
def find_orfs_in_all_6_reading_frames(rna, allow_overlap = False):
    all_orfs = []
    total_length = len(rna)
    # Forward strand frame (0,1,2)
    for frame in range(3):
        frame_rna = rna[frame:]
        orfs = find_orfs_in_frame(frame_rna, offset = frame, strand ='+', frame = frame,
                                   allow_overlap = allow_overlap, total_length = total_length)
        all_orfs.extend(orfs)
        
    
    # Reverse strande frames (0,1,2)
    rev_rna = reverse_complement(rna)
    for frame in range(3):
        frame_rna = rev_rna[frame:]
        orfs = find_orfs_in_frame(frame_rna, offset = frame, strand = '-', frame = frame,
                                   allow_overlap = allow_overlap, total_length = total_length)
        all_orfs.extend(orfs)
    return all_orfs

# ------------------------------
# Main Driver
# ------------------------------
def main():
    print ("=== ORF Finder & Sequence Analyzer ===\n")
    seq = input("Enter your DNA or RNA sequence: ").strip().upper()

    seq_type = sequence_type(seq)
    if seq_type == "Unknown or mixed invalid input!":
        print ("Invalid sequence. Please input a valid DNA or RNA sequence.")
        return
    print (f"\nSequence type detected: {seq_type}")
    
    rna = dna_transcription(seq) if seq_type == "DNA" else seq

    print(f"GC content: {gc_content(seq)}%")
    allow_overlap = input ("Allow overlapping ORFs? (y/n): ").strip().lower() == 'y'
    orfs = find_orfs_in_all_6_reading_frames(rna, allow_overlap = allow_overlap)
    if not orfs:
        print ("No ORFs found.")
    else:
        print (f"\nFound {len(orfs)} ORFs:\n")
        for i, orf in enumerate(sorted (orfs, key=lambda x: (x['strand'], x['frame'], x['start'])),1):
            print (f"ORF {i}:")
            print (f" Strand : {orf['strand']}")
            print (f" Frame : {orf['frame']}")
            print (f" Start : {orf['start']}")
            print (f" End : {orf['end']}")
            print (f" Length : {orf['length']}nt")
            print (f" Protein : {orf['protein']}\n")
    
    if input("Show codon usage frequncy? (y/n): ").strip().lower() == 'y':
        usage = codon_usage (rna)
        print ("\nCodon usage (%): ")
        for codon in sorted(ALL_CODONS):
            freq = usage[codon]
            if freq > 0:
                print (f"  {codon}:{freq}%")

    if input ("Save ORF protein sequences to file: (y/n) ").strip().lower() == 'y':
        with open("orf_protein.txt", "w") as f:
             for i, orf in enumerate(orfs, 1):
                 f.write(f">ORF_{i} | Strand={orf['strand']} | Frame={orf['frame']} | Start={orf['start']} | End={orf['end']}\n")
                 f.write(orf['protein'] + '\n')
        print ("Saved to orf_protein.txt")

if __name__ == "__main__":
    main()