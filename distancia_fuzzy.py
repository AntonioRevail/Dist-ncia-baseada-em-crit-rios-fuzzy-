from Bio import SeqIO
import math
import csv
from scipy.optimize import fsolve
from datasketch import MinHash

# ---------------------------------------------
# Funções de pré-processamento
# ---------------------------------------------

def filtrar_bases_validas(seq):
    return ''.join(base for base in seq.upper() if base in 'ACGT')

# ---------------------------------------------
# Critérios fuzzy
# ---------------------------------------------

def gc_content(seq):
    seq = filtrar_bases_validas(seq)
    if not seq:
        return 0
    gc_count = seq.count('G') + seq.count('C')
    return gc_count / len(seq)

def shannon_entropy(seq):
    seq = filtrar_bases_validas(seq)
    total = len(seq)
    if total == 0:
        return 0
    freqs = {base: seq.count(base)/total for base in 'ACGT'}
    entropy = -sum(p * math.log2(p) for p in freqs.values() if p > 0)
    return entropy / 2  # normalizado

def cg_pattern(seq):
    seq = seq.upper()
    digrams = [seq[i:i+2] for i in range(len(seq) - 1)]
    valid_digrams = [d for d in digrams if set(d) <= set('ACGT')]
    if not valid_digrams:
        return 0
    return sum(1 for d in valid_digrams if d == 'CG') / len(valid_digrams)

def kmers(seq, k=2):
    seq = seq.upper()
    return set(seq[i:i+k] for i in range(len(seq) - k + 1) if set(seq[i:i+k]) <= set('ACGT'))

# ---------------------------------------------
# MinHash para estimar similaridade de Jaccard
# ---------------------------------------------

def create_minhash(kmer_set, num_perm=128):
    mh = MinHash(num_perm=num_perm)
    for kmer in kmer_set:
        mh.update(kmer.encode('utf8'))
    return mh

def minhash_similarity(mh1, mh2):
    return mh1.jaccard(mh2)

# ---------------------------------------------
# Integral de Sugeno
# ---------------------------------------------

def solve_lambda(measures):
    def equation(lmb):
        prod = 1
        for mu in measures:
            prod *= (1 + lmb * mu)
        return prod - (1 + lmb)
    sol = fsolve(equation, x0=0.1)
    return sol[0]

def sugeno_measure(subset, lambd, mu_values):
    if not subset:
        return 0
    if len(subset) == 1:
        return mu_values[subset[0]]
    else:
        A = subset[:-1]
        B = subset[-1:]
        muA = sugeno_measure(A, lambd, mu_values)
        muB = mu_values[B[0]]
        return muA + muB + lambd * muA * muB

# ---------------------------------------------
# Função principal
# ---------------------------------------------

def fuzzy_distance_fasta(input_fasta, output_dist_csv, output_feats_csv=None):
    records = list(SeqIO.parse(input_fasta, "fasta"))
    seqs = {rec.id: str(rec.seq) for rec in records}

    kmer_sets = {name: kmers(seq) for name, seq in seqs.items()}
    minhashes = {name: create_minhash(kmer_sets[name]) for name in seqs}

    features = {}
    for name, seq in seqs.items():
        features[name] = {
            'GC': gc_content(seq),
            'Entropy': shannon_entropy(seq),
            'CG': cg_pattern(seq),
        }

    for name in seqs:
        others = [n for n in seqs if n != name]
        sims = [minhash_similarity(minhashes[name], minhashes[o]) for o in others]
        features[name]['Kmer'] = sum(sims) / len(sims) if sims else 0

    # Exporta CSV com os critérios fuzzy
    if output_feats_csv:
        with open(output_feats_csv, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['ID', 'GC_Content', 'Shannon_Entropy', 'CG_Pattern', 'Kmer_Similarity'])
            for name in features:
                row = [name] + [features[name][key] for key in ['GC', 'Entropy', 'CG', 'Kmer']]
                writer.writerow(row)

    # Matriz de distância fuzzy com Sugeno
    seq_ids = list(seqs.keys())
    distance_matrix = []

    for i in seq_ids:
        row = []
        for j in seq_ids:
            if i == j:
                row.append(0.0)
                continue

            mu1 = features[i]
            mu2 = features[j]
            h = {k: 1 - abs(mu1[k] - mu2[k]) for k in ['GC', 'Entropy', 'CG', 'Kmer']}
            h_sorted = sorted(h.items(), key=lambda x: x[1])

            mu_C = {'GC': 0.2, 'Entropy': 0.3, 'CG': 0.2, 'Kmer': 0.3}
            lambd = solve_lambda(list(mu_C.values()))

            sugeno_vals = []
            crits = list(h_sorted)
            for idx in range(len(crits)):
                subset_keys = [k for k, _ in crits[idx:]]
                hi = crits[idx][1]
                mu_sub = sugeno_measure(subset_keys, lambd, mu_C)
                sugeno_vals.append(min(hi, mu_sub))

            sugeno_integral = max(sugeno_vals)
            distance = 1 - sugeno_integral
            row.append(round(distance, 5))
        distance_matrix.append(row)

    # Exporta matriz de distância
    with open(output_dist_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow([''] + seq_ids)
        for name, row in zip(seq_ids, distance_matrix):
            writer.writerow([name] + row)

# ---------------------------------------------
# Exemplo de uso
# ---------------------------------------------

fuzzy_distance_fasta("input.fasta", "output.csv", "caracteristicas_fuzzy.csv")
