"""
To get how many Pattern in DNA sequence
(string Text, string Pattern)
"""
import random


def pattern_count(text, pattern):
    # pattern is the most frequent k-mer
    p_count = 0
    n = len(text)
    k = len(pattern)
    # go through all Text to get how many Pattern one by one in text
    for i in range(n - k + 1):
        if text[i:i + k] == pattern:
            p_count += 1
    return p_count


"""
The Frequent Words Problem
get all k-mer in text, and its count appear
"""


def frequency_map(text, k):
    freq = {}
    n = len(text)
    # go through each k-mers existed
    for i in range(n - k + 1):
        pattern = text[i:i + k]
        # if there is existed pattern, add by 1
        if pattern in freq:
            freq[pattern] += 1
        # if not, add it to freq as 1
        else:
            freq[pattern] = 1
    return freq


"""
find the most frequent k-mers in text
"""


def most_frequent(text, k):
    words = []
    freq = frequency_map(text, k)
    # find max value in freq
    m = max(freq.values())
    # find word same count as max value
    for word in freq:
        if freq[word] == m:
            words.append(word)
    return words


"""
each DNA strand is read in the 5' --> 3'
Reverse Complement of a DNA string Pattern, For example, "AGTCGCATAGT" -> "TCAGCGTATCA" -> "ACTATGCGACT"
"""


def reverse_complement(pattern):
    # reverse from end to begin of pattern
    pattern = pattern[::-1]
    pattern = complement(pattern)
    return pattern


def complement(pattern):
    # make complement_string and match them one by one char
    complement_string = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    string = ''.join([complement_string[char] for char in pattern])
    pattern = string
    return pattern


"""
find position of pattern find in DNA strand
"""


def pattern_matching(pattern, genome):
    positions = []
    k = len(pattern)
    n = len(genome)
    # find the all positions
    for i in range(n - k + 1):
        if genome[i:i + k] == pattern:
            positions.append(i)
    return positions


"""
The DNA replication is semiconservative hypothesis, two pairs of strand, each pair 
have one new and one old 

DNA replication process: 1. To separate the two strands. (helicase enzyme)
                         2. Primase enzyme make small piece of RNA, called a primer
                            make starting point for the new strand of DNA.
                         3. DNA polymerase make new strand of DNA from 5'end to 3'end.
                         4. opposite strand, the lagging strand cannot made in continuous way,
                            so DNA polymerase can only make this strand in a series of small
                            chunks, called Okazaki fragments. Each fragment is started with 
                            an RNA primer, and then DNA polymerase add chunks in 5'ene to 3'end 
                         5. DNA ligase seals up the fragments of DNA.                          
"""

"""
tracking of the total number of occurrences of symbol in half of the sequences 
"""


def symbol_array(genome, symbol):
    array = {}
    n = len(genome)
    # extend genome with more half of genome
    extend_genome = genome + genome[0:n // 2]
    # count symbol will appear in rest of extend_genome
    for i in range(n):
        array[i] = pattern_count(extend_genome[i:i + (n // 2)], symbol)
    return array


"""
faster way to count total number of occurrences of symbol in half of the sequences
determine where is ori, maximum (forward half-strand), and the minimum (reverse half-strand)
"""


def faster_symbol_array(genome, symbol):
    array = {}
    n = len(genome)
    # extend genome with more half of genome
    extend_genome = genome + genome[0:n // 2]
    # make first array count first
    array[0] = pattern_count(extend_genome[0:(n // 2)], symbol)
    for i in range(1, n):
        array[i] = array[i - 1]
        # if previous or the char after half of genome sequence is that symbol, at most 1
        if extend_genome[i - 1] == symbol:
            array[i] -= 1
        if extend_genome[i + (n // 2) - 1] == symbol:
            array[i] += 1

    return array


"""
searching for ori by G-C > 0 (forward half-strand) G-C < 0 (reverse half-strand)
tracking difference between C and G 
"""


def skew_array(genome):
    skew_list = [0]
    for char in genome:
        if char == 'C':
            skew_list.append(skew_list[-1] - 1)
        elif char == 'G':
            skew_list.append(skew_list[-1] + 1)
        else:
            skew_list.append(skew_list[-1])
    return skew_list


"""
find the minimum value position 
"""


def minimum_skew(genome):
    skew_list = skew_array(genome)
    # find all min in skew list, so that we can find ori
    min_value = min(skew_list)
    positions = []
    for i in range(len(skew_list)):
        if skew_list[i] == min_value:
            positions.append(i)
    return positions


"""
Find difference between two k-mers
"""


def hamming_distance(strand1, strand2):
    d_count = 0
    for i in range(len(strand1)):
        if strand1[i] != strand2[i]:
            d_count += 1
    return d_count


"""
Find difference of d by compare with pattern to genome
"""


def approximate_pattern_matching(text, pattern, d):
    positions = []
    # compared one by one to whole text, only remind at most d position
    for i in range(len(text) - len(pattern) + 1):
        d_count = hamming_distance(text[i:i + len(pattern)], pattern)
        if d_count <= d:
            positions.append(i)
    return positions


"""
Gene generation: DNA transcribed into a strand of RNA composed of four ribonucleotides:
                 A,T,C,U.
                 Then, RNA transcript is translated into the amino acid sequence of protein

plant cell: LHY, CCA1, TOC1, the clock's master. light activates LHY, CCA1, and triggering 
            TOC1. with the light diminishes, so does the production of LHY, CCA1, and they 
            do not repress TOC1, the the peaks of TOC1 at night, and starts promoting LHY, CCA1
                  
"""

"""
Scoring Motifs: count how many char have in each column 
"""


def count(motifs):
    motif_count = {}
    # initially the motif_count with same row length
    for char in "ACGT":
        motif_count[char] = [0] * len(motifs[0])
    col = len(motifs)

    # count them with go through each column
    for i in range(len(motifs[0])):
        for j in range(col):
            motif_count[motifs[j][i]][i] += 1
    return motif_count


"""
count the percent value after we count scoring motifs
"""


def profile(motifs):
    motif_count = count(motifs)
    profile_array = {}
    t = len(motifs)
    # count how many percent of each nucleotide in each column they were in
    for nucleotide in motif_count.keys():
        profile_array[nucleotide] = [base_count / t for base_count in motif_count[nucleotide]]
    return motif_count


"""
find the highest appearance of nucleotide in each column
"""


def consensus(motifs):
    motif_count = count(motifs)
    motif_consensus = ""
    # find the highest nucleotide each column
    for i in range(len(motifs[0])):
        max_value = 0
        symbol = ""
        for nucleotide in motif_count.keys():
            if motif_count[nucleotide][i] > max_value:
                max_value = motif_count[nucleotide][i]
                symbol = nucleotide
        motif_consensus += symbol
    return motif_consensus


"""
sum the time exception of highest nucleotides in each column
for example "AAACC" = 2 because C is differ about highest "A"
"""


def score(motifs):
    motif_count = count(motifs)
    motif_consensus = consensus(motifs)
    total = 0
    for i in range(len(motif_consensus)):
        total += motif_count[motif_consensus[i]][i]
    dismatch = len(motifs) * len(motifs[0]) - total

    return dismatch


"""
Greedy Motif Search: "most attractive", often fast heuristics that trade accuracy
for speed in order to find an approximate solution
"""

"""
calculate NF-kB, time all percentage depend on profile that existed 
"""


def pr(text, profile):
    # Cromwell's rule
    p = 1
    for i in range(len(text)):
        p *= profile[text[i]][i]
    return p


"""
find the k-mers that have highest Pr (NF-kB) by looking through the whole text 
"""


def profile_most_probable(text, k, profile):
    most_p = -1.0
    most_pattern = ""
    # go through all k-mers that in text, find highest pr, and its pattern
    for i in range(len(text) - k + 1):
        if most_p < pr(text[i: i + k], profile):
            most_p = pr(text[i: i + k], profile)
            most_pattern = text[i: i + k]
    return most_pattern


"""
Find the most possible k-mers in motifs by greedy search
    calculate profile of dna[0] by k-mers, and then find k-mers in rest of row in DNA
    Find the smallest score of motifs format
"""


def greedy_motif_search(dna, k, t):
    # initialize the best motifs as first k-mer each row
    best_motifs = [dna[i][0:k] for i in range(t)]
    for i in range(len(dna[0]) - k + 1):
        motifs = [dna[0][i:i + k]]
        for j in range(1, t):
            # dna[0] is starting positions in the first sequence, so we use dna[0]
            p = profile(motifs[0:j])
            # base on profile above to find most probable motifs each row
            motifs.append(profile_most_probable(dna[j], k, p))

        if score(motifs) < score(best_motifs):
            best_motifs = motifs
    return best_motifs


"""
pseudocounts:
often add 1 to each element of count(motifs), so we did not have 0 probability
"""


def count_with_pseudocounts(motifs):
    t = len(motifs)
    m_count = {}
    # different with count, initial all with 1 instead of 0
    for symbol in "ACGT":
        m_count[symbol] = [1] * len(motifs[0])

    for i in range(len(motifs[0])):
        for j in range(t):
            m_count[motifs[j][i]][i] += 1

    return m_count


"""
profile with pseudocounts, divided (t + 4)
"""


def profile_pseudocounts(motifs):
    t = len(motifs)
    m_count = count_with_pseudocounts(motifs)
    profile_array = {}

    # differ with profile, need to divide by (t + 4)
    for symbol in m_count.keys():
        profile_array[symbol] = [char / (t + 4) for char in m_count[symbol]]

    return profile_array


def greedy_with_pseudocounts(dna, k, t):
    best_motifs = [dna[i][0:k] for i in range(t)]
    n = len(dna[0])

    for i in range(n - k + 1):
        motifs = [dna[0][i:i + k]]
        for j in range(1, t):
            # make profile for each updated motifs after it been appended new row
            p = profile_pseudocounts(motifs[0:j])
            motifs.append(profile_most_probable(dna[j], k, p))

        if score(motifs) < score(best_motifs):
            best_motifs = motifs
    return best_motifs


"""
Monte Carlo algorithm: not return exact solutions, but it can quickly find 
approximate solutions. (speed higher)
"""


def motifs(profile, dna):
    k = len(profile['A'])
    best_motifs = []
    # find each row best k-mers output based on the profile given
    for i in range(len(dna)):
        best_motifs.append(profile_most_probable(dna[i], k, profile))

    return best_motifs


"""
randomly pick up any k-mers each row
"""


def random_motifs(dna, k, t):
    r_motifs = []
    for i in range(t):
        start = random.randint(0, len(dna[0]) - k)
        r_motifs.append(dna[i][start:start + k])

    return r_motifs


def randomized_motif_search(dna, k, t):
    r_motifs = random_motifs(dna, k, t)
    best_motifs = r_motifs

    # random pick up a motifs, and then do the loop, until best motifs has the smallest score
    while True:
        profile = profile_pseudocounts(r_motifs)
        r_motifs = motifs(profile, dna)
        if score(r_motifs) < score(best_motifs):
            best_motifs = r_motifs
        else:
            return best_motifs


"""
GibbsSampler: each time, single k-mer decides to either keep or replace, move 
with more caution
"""

"""
make sum of probabilities as 1
"""


def normalize(probabilities):
    total = sum(probabilities.values())

    for char in probabilities.keys():
        probabilities[char] /= total

    return probabilities


def weighted_die(probabilities):
    p = random.uniform(0, 1)
    # get floating point number to get k_mer in range
    for k_mer, prob in probabilities.items():
        p -= prob
        if p <= 0:
            return k_mer
    return None


def profile_generated_string(text, profile, k):
    n = len(text)
    probabilities = {}

    # get all possible k-mers probabilities in the text (only one row)
    for i in range(n - k + 1):
        probabilities[text[i:i + k]] = pr(text[i:i + k], profile)
    # then, normalize this p, and get most possible k-mer randomly
    probabilities = normalize(probabilities)
    return weighted_die(probabilities)


def gibbs_sampler(dna, k, t, N):
    r_motifs = random_motifs(dna, k, t)
    best_motifs = r_motifs
    # each time, we pick up one row randomly, and remove it from motifs
    for i in range(1, N):
        r_row = random.randint(1, t)
        except_row = r_motifs.pop(r_row-1)
        # removed motifs as profile, and get new k-mer in that removed row based on this profile
        r_profile = profile_pseudocounts(r_motifs)
        except_row = profile_generated_string(except_row, r_profile, k)
        # add new row to the motifs again, and compare with best_motifs
        r_motifs.insert(r_row-1, except_row)

        if score(r_motifs) < score(best_motifs):
            best_motifs = r_motifs
    return best_motifs


dna = [
    "CGCCCCTCTCGGGGGTGTTCAGTAAACGGCCA",
    "GGGCGAGGTATGTGTAAGTGCCAAGGTGCCAG",
    "TAGTACCGAGACCGAAAGAAGTATACAGGCGT",
    "TAGATCAAGTTTCAGGTGCACGTCGGTGAACC",
    "AATCCACCAGCTCCACGTGCAATGTTGGCCTA"
]

print(gibbs_sampler(dna, 8, 5, 100))
