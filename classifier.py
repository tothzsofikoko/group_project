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