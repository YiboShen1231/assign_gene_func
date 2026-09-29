def global_alignment(seq1, seq2, scoring_function):
    """Global sequence alignment using the Needleman–Wunsch algorithm.

    Indels should be denoted with the "-" character.

    Parameters
    ----------
    seq1: str
        First sequence to be aligned.
    seq2: str
        Second sequence to be aligned.
    scoring_function: Callable

    Returns
    -------
    str
        First aligned sequence.
    str
        Second aligned sequence.
    float
        Final score of the alignment.

    Examples
    --------
    >>> global_alignment("abracadabra", "dabarakadara", lambda x, y: [-1, 1][x == y])
    ('-ab-racadabra', 'dabarakada-ra', 5.0)

    Other alignments are not possible.

    """

    n = len(seq1)
    m = len(seq2)

    score = [[0.0] * (m + 1) for _ in range(n + 1)]
    trace = [[None] * (m + 1) for _ in range(n + 1)]

    for i in range(1, n + 1):
        score[i][0] = (
            score[i - 1][0]
            + scoring_function(seq1[i - 1], "-")
        )
        trace[i][0] = "U"

    for j in range(1, m + 1):
        score[0][j] = (
            score[0][j - 1]
            + scoring_function("-", seq2[j - 1])
        )
        trace[0][j] = "L"

    for i in range(1, n + 1):
        for j in range(1, m + 1):

            diagonal = (
                score[i - 1][j - 1]
                + scoring_function(seq1[i - 1], seq2[j - 1])
            )

            up = (
                score[i - 1][j]
                + scoring_function(seq1[i - 1], "-")
            )

            left = (
                score[i][j - 1]
                + scoring_function("-", seq2[j - 1])
            )

            best_score = max(diagonal, up, left)
            score[i][j] = best_score

            if best_score == diagonal:
                trace[i][j] = "D"
            elif best_score == up:
                trace[i][j] = "U"
            else:
                trace[i][j] = "L"

    final_score = float(score[n][m])

    aligned_seq1 = []
    aligned_seq2 = []

    i = n
    j = m

    while i > 0 or j > 0:

        direction = trace[i][j]

        if direction == "D":
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append(seq2[j - 1])
            i -= 1
            j -= 1

        elif direction == "U":
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append("-")
            i -= 1

        else:
            aligned_seq1.append("-")
            aligned_seq2.append(seq2[j - 1])
            j -= 1

    aligned_seq1 = "".join(reversed(aligned_seq1))
    aligned_seq2 = "".join(reversed(aligned_seq2))

    return aligned_seq1, aligned_seq2, final_score


def local_alignment(seq1, seq2, scoring_function):
    """Local sequence alignment using the Smith-Waterman algorithm.

    Indels should be denoted with the "-" character.

    Parameters
    ----------
    seq1: str
        First sequence to be aligned.
    seq2: str
        Second sequence to be aligned.
    scoring_function: Callable

    Returns
    -------
    str
        First aligned sequence.
    str
        Second aligned sequence.
    float
        Final score of the alignment.

    Examples
    --------
    >>> local_alignment("pending itch", "unending glitch", lambda x, y: [-1, 1][x == y])
    ('ending --itch', 'ending glitch', 9.0)

    Other alignments are not possible.

    """
    raise NotImplementedError()


## This is an example scoring function, you should implement a version which uses a scoring matrix 
def scoring_function_simple(aa_i,aa_j):
    score = [-1, 1][aa_i == aa_j]
    return (score)

from Bio.Align import substitution_matrices

blosum62 = substitution_matrices.load("BLOSUM62")

def scoring_function_blosum62(aa_i, aa_j):
    if aa_i == "-" or aa_j == "-":
        return -4
    return blosum62[aa_i, aa_j]
