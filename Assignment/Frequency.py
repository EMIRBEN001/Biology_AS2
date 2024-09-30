# Normal Codon table for mRNA to amino acids
amino_acid_to_codon = {
    'Phe': ['UUU', 'UUC'], 'Leu': ['UUA', 'UUG', 'CUU', 'CUC', 'CUA', 'CUG'],
    'Ile': ['AUU', 'AUC', 'AUA'], 'Met': ['AUG'], 'Val': ['GUU', 'GUC', 'GUA', 'GUG'],
    'Ser': ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'], 'Pro': ['CCU', 'CCC', 'CCA', 'CCG'],
    'Thr': ['ACU', 'ACC', 'ACA', 'ACG'], 'Ala': ['GCU', 'GCC', 'GCA', 'GCG'],
    'Tyr': ['UAU', 'UAC'], 'His': ['CAU', 'CAC'], 'Gln': ['CAA', 'CAG'],
    'Asn': ['AAU', 'AAC'], 'Lys': ['AAA', 'AAG'], 'Asp': ['GAU', 'GAC'],
    'Glu': ['GAA', 'GAG'], 'Cys': ['UGU', 'UGC'], 'Trp': ['UGG'], 'Arg': ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'],
    'Gly': ['GGU', 'GGC', 'GGA', 'GGG'], 'Stop': ['UAA', 'UAG', 'UGA']
}

# Map amino acid single letters to their names the one that looks like a wheel mandala thing
single_letter_to_amino_acid = {
    'F': 'Phe', 'L': 'Leu', 'I': 'Ile', 'M': 'Met', 'V': 'Val', 'S': 'Ser',
    'P': 'Pro', 'T': 'Thr', 'A': 'Ala', 'Y': 'Tyr', 'H': 'His', 'Q': 'Gln',
    'N': 'Asn', 'K': 'Lys', 'D': 'Asp', 'E': 'Glu', 'C': 'Cys', 'W': 'Trp',
    'R': 'Arg', 'G': 'Gly'
}

# Function to get all possible combinations of codons
def generate_mrna_sequences(codon_possibilities):
    if len(codon_possibilities) == 1:
        return codon_possibilities[0]
    else:
        previous_combinations = generate_mrna_sequences(codon_possibilities[:-1])
        new_combinations = []
        for combination in previous_combinations:
            for codon in codon_possibilities[-1]:
                new_combinations.append(combination + codon)
        return new_combinations

# Function to find codon frequency for input amino acid sequence
def amino_acid_to_rna_frequencies(amino_acid_seq):
    # Convert to uppercase to handle both lowercase and uppercase inputs
    amino_acid_seq = amino_acid_seq.upper()
    
    # Step 1: Get all possible codons for the amino acid sequence
    codon_possibilities = []
    for amino_acid in amino_acid_seq:
        if amino_acid in single_letter_to_amino_acid:
            aa_name = single_letter_to_amino_acid[amino_acid]
            codons = amino_acid_to_codon[aa_name]
            codon_possibilities.append(codons)
        else:
            print(f"Invalid amino acid symbol: {amino_acid}")
            return
    
    # Step 2: Generate all possible mRNA sequences
    possible_mrna_sequences = generate_mrna_sequences(codon_possibilities)
    
    # Step 3: Count the codon frequencies for each mRNA sequence and print in a formatted way
    for mrna in possible_mrna_sequences:
        print(f"mRNA = {mrna}")
        codon_frequencies = {}
        for i in range(0, len(mrna), 3):
            codon = mrna[i:i+3]
            if codon in codon_frequencies:
                codon_frequencies[codon] += 1
            else:
                codon_frequencies[codon] = 1
        for codon, frequency in codon_frequencies.items():
            print(f"{codon} = {frequency}")
        print()  # For spacing between sequences

# Example usage
amino_acid_seq = input("Enter an amino acid sequence : ")

amino_acid_to_rna_frequencies(amino_acid_seq)
