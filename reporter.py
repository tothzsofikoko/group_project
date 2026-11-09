# ---------------------------------------------------------------------------
# Role 4: Reporter
# ---------------------------------------------------------------------------

def print_report(counts):
    """Print each category with its count and percentage of all trials."""
    total = sum(counts.values())

    print(f"Total mutations: {total}")
    print(f"{'Effect':<15} {'Count':>8} {'Percentage':>12}")

    for effect, count in counts.items():
        percentage = count / total * 100 if total > 0 else 0
        print(f"{effect:<15} {count:>8} {percentage:>11.2f}%")


def plot_results(counts, title="Effects of random mutations"):
    """Draw a bar chart of the counts."""
    import matplotlib.pyplot as plt

    plt.figure(figsize=(8, 5))
    plt.bar(list(counts.keys()), list(counts.values()))

    plt.title(title)
    plt.xlabel("Effect on protein")
    plt.ylabel("Number of mutations")

    plt.tight_layout()
    plt.show()


def demo_one_mutation(seq):
    """Show one mutation step by step."""
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