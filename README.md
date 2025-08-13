# ORF Finder and Sequence Analyzer

This Python script analyzes DNA/RNA sequences. It can:

- Detect sequence type (DNA or RNA)
- Transcribe DNA to RNA
- Find open reading frames (ORFs) in all 6 frames
- Translate RNA to protein
- Show codon usage frequency
- Save ORFs to a file

## How to Run

1. Make sure you have Python installed.  
2. Run the script in your terminal:

```bash
python orf_finder.py
```


## Example
**Input:**
```
Enter your DNA or RNA sequence: ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGA
Allow overlapping ORFs? (y/n): n
Show codon usage frequency? (y/n): n
Save ORF protein sequences to file: (y/n): y

``` 

**Output:**
```
=== ORF Finder & Sequence Analyzer ===

Sequence type detected: DNA
GC content: 56.41%

Found 1 ORFs:

ORF 1:
Strand : +
Frame : 0
Start : 0
End : 36
Length : 36nt
Protein : MAIVMGR*KGAR

Saved to orf_protein.txt
```
