import random

CODON_TABLE = {
    "TTT": "Phe", "TTC": "Phe", "TTA": "Leu", "TTG": "Leu",
    "CTT": "Leu", "CTC": "Leu", "CTA": "Leu", "CTG": "Leu",
    "ATT": "Ile", "ATC": "Ile", "ATA": "Ile", "ATG": "Met",
    "GTT": "Val", "GTC": "Val", "GTA": "Val", "GTG": "Val",
    "TCT": "Ser", "TCC": "Ser", "TCA": "Ser", "TCG": "Ser",
    "CCT": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",
    "ACT": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",
    "GCT": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",
    "TAT": "Tyr", "TAC": "Tyr", "TAA": "Stop", "TAG": "Stop",
    "CAT": "His", "CAC": "His", "CAA": "Gln", "CAG": "Gln",
    "AAT": "Asn", "AAC": "Asn", "AAA": "Lys", "AAG": "Lys",
    "GAT": "Asp", "GAC": "Asp", "GAA": "Glu", "GAG": "Glu",
    "TGT": "Cys", "TGC": "Cys", "TGA": "Stop", "TGG": "Trp",
    "CGT": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg",
    "AGT": "Ser", "AGC": "Ser", "AGA": "Arg", "AGG": "Arg",
    "GGT": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly",
}

BASES = "ACGT"

TEST_SEQUENCE = (
    "ATGGCCTTTAAAGGCTGGCTCGAGAAACCC"
    "GGGTACCATCGTAGCGATCAAATTGTGTAG"
)

def translate(seq):
    """Translate DNA into a list of amino acids, stopping at the first stop codon.
    Leftover bases at the end (fewer than 3) are ignored."""
    protein = []
    for i in range(0, len(seq) - 2, 3):
        amino_acid = CODON_TABLE[seq[i:i + 3]]
        if amino_acid == "Stop":
            break
        protein.append(amino_acid)
    return protein

##mutator
def substitute(seq):
    pos = random.randrange(len(seq))
    new_base = random.choice([b for b in BASES if b != seq[pos]])
    return seq[:pos] + new_base + seq[pos + 1:], pos


def insert(seq):
    pos = random.randrange(len(seq) + 1)
    new_base = random.choice(BASES)
    return seq[:pos] + new_base + seq[pos:], pos


def delete(seq):
    pos = random.randrange(len(seq))
    return seq[:pos] + seq[pos + 1:], pos

mutations = ["substitution", "insertion", "deletion"]

def random_mutation(seq):
    mut = random.choice(mutations)
    if mut == "substitution":
        new_seq, pos = substitute(seq)
    elif mut == "insertion":
        new_seq, pos = insert(seq)
    elif mut == "deletion":
        new_seq, pos = delete(seq)

    return (new_seq, mut, pos)
    pass


##classifier
def classify(original_protein, mutated_protein, kind = ["substitution", "insertion", "deletion"]):
    if mutated_protein == original_protein:
        effect = "Mutation type: silent"
    elif kind in ("insertion", "deletion"):
        effect = "Mutation type: frameshift"
    elif len(mutated_protein) < len(original_protein):
        effect = "Mutation type: nonsense"
    else:
        effect = "Mutation type: missense"
    return(original_protein, mutated_protein, kind, effect)
    pass

##simulator
def run_simulation(seq, trials, mode="any"):
    counts = {
        "silent": 0,
        "missense": 0,
        "nonsense": 0,
        "frameshift": 0
    }
    original_protein = translate(seq)
    for i in range(trials):
        if mode == "substitution":
            new_seq, position = substitute(seq)
            kind = "substitution"
        else:
            new_seq, kind, position = random_mutation(seq)
        mutated_protein = translate(new_seq)
        result = classify(original_protein, mutated_protein, kind)
        counts[result] += 1
    return counts
    pass

##reporter
def print_report(counts):
    """Print each category with its count and percentage of all trials."""
    total = sum(counts.values())

    print(f"Total mutations: {total}")
    print(f"{'Effect':<15} {'Count':>8} {'Percentage':>12}")

    for effect, count in counts.items():
        percentage = count / total * 100 if total > 0 else 0
        print(f"{effect:<15} {count:>8} {percentage:>11.2f}%")


def plot_results(counts, title="Effects of random mutations"):
    import matplotlib.pyplot as plt

    plt.figure(figsize=(8, 5))
    plt.bar(list(counts.keys()), list(counts.values()))

    plt.title(title)
    plt.xlabel("Effect on protein")
    plt.ylabel("Number of mutations")

    plt.tight_layout()
    plt.show()


def demo_one_mutation(seq):
    original_protein = translate(seq)

    new_seq, kind, position = random_mutation(seq)
    mutated_protein = translate(new_seq)

    effect = classify(original_protein, mutated_protein, kind)

    print("=== Single mutation demo ===")
    print("Original DNA:", seq)
    print("Original protein:", original_protein)
    print("Mutation type:", kind)
    print("Position (0-based index):", position)

    if kind == "substitution":
        print(f"Change: {seq[position]} -> {new_seq[position]}")
    elif kind == "insertion":
        print(f"Change: inserted {new_seq[position]}")
    elif kind == "deletion":
        print(f"Change: deleted {seq[position]}")

    print("Mutated DNA:", new_seq)
    print("Mutated protein:", mutated_protein)
    print("Classification:", effect)
    print()


##demo
if __name__ == "__main__":
    demo_one_mutation(TEST_SEQUENCE)

    print("=== Substitutions only ===")
    sub_counts = run_simulation(TEST_SEQUENCE, 1000, mode="substitution")
    print_report(sub_counts)
    print()

    print("=== All mutation types ===")
    all_counts = run_simulation(TEST_SEQUENCE, 1000, mode="any")
    print_report(all_counts)

    plot_results(sub_counts, "Effects of 1,000 random substitutions")
