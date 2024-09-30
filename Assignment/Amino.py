# Codon table for mRNA to amino acids
codon_table = {
    'UUU': 'Phe', 'UUC': 'Phe', 'UUA': 'Leu', 'UUG': 'Leu',
    'UCU': 'Ser', 'UCC': 'Ser', 'UCA': 'Ser', 'UCG': 'Ser',
    'UAU': 'Tyr', 'UAC': 'Tyr', 'UGU': 'Cys', 'UGC': 'Cys',
    'UGG': 'Trp', 'CUU': 'Leu', 'CUC': 'Leu', 'CUA': 'Leu',
    'CUG': 'Leu', 'CCU': 'Pro', 'CCC': 'Pro', 'CCA': 'Pro',
    'CCG': 'Pro', 'CAU': 'His', 'CAC': 'His', 'CAA': 'Gln',
    'CAG': 'Gln', 'CGU': 'Arg', 'CGC': 'Arg', 'CGA': 'Arg',
    'CGG': 'Arg', 'AUU': 'Ile', 'AUC': 'Ile', 'AUA': 'Ile',
    'AUG': 'Met', 'ACU': 'Thr', 'ACC': 'Thr', 'ACA': 'Thr',
    'ACG': 'Thr', 'AAU': 'Asn', 'AAC': 'Asn', 'AAA': 'Lys',
    'AAG': 'Lys', 'AGU': 'Ser', 'AGC': 'Ser', 'AGA': 'Arg',
    'AGG': 'Arg', 'GUU': 'Val', 'GUC': 'Val', 'GUA': 'Val',
    'GUG': 'Val', 'GCU': 'Ala', 'GCC': 'Ala', 'GCA': 'Ala',
    'GCG': 'Ala', 'GAU': 'Asp', 'GAC': 'Asp', 'GAA': 'Glu',
    'GAG': 'Glu', 'GGU': 'Gly', 'GGC': 'Gly', 'GGA': 'Gly',
    'GGG': 'Gly', 'UAA': 'Stop', 'UAG': 'Stop', 'UGA': 'Stop'
}

# Function to convert DNA sequence to mRNA and amino acid sequence
def dna_to_protein(dna_sequence):
    dna_sequence = dna_sequence.upper()

    # Step 1: Generate Complementary DNA strand
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    comp_dna_sequence = ''.join([complement[base] for base in dna_sequence])
    
    # Step 2: Convert Complementary DNA to mRNA (T -> U)
    mrna_sequence = comp_dna_sequence.replace('T', 'U')
    
    # Step 3: Translate mRNA to Amino Acids
    protein_sequence = []
    for i in range(0, len(mrna_sequence), 3):
        codon = mrna_sequence[i:i+3]
        if codon in codon_table:
            protein_sequence.append(codon_table[codon])
    
    return comp_dna_sequence, mrna_sequence, protein_sequence

# Example usage
dna_sequence = input("Enter a DNA sequence (multiple of 3): ")

# Ensure that the sequence is valid (multiple of 3 and contains only A, T, C, G)
if len(dna_sequence) % 3 != 0 or not all(base.upper() in 'ATCG' for base in dna_sequence):
    print("Invalid DNA sequence. Please ensure it is a multiple of 3 and contains only A, T, C, and G.")
else:
    comp_dna_sequence, mrna_sequence, protein_sequence = dna_to_protein(dna_sequence)
    
    print(f"Input DNA: {dna_sequence.upper()}")
    print(f"Complementary DNA: {comp_dna_sequence}")
    print(f"mRNA: {mrna_sequence}")
    print(f"Amino Acid Sequence: {' - '.join(protein_sequence)}")
