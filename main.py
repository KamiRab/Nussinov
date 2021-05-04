import numpy as np


def score(valueCG, valueAU, base1, base2):
    if (base1 == 'A' and base2 == 'U') or (base1 == 'U' and base2 == 'A'):
        return valueAU
    if (base1 == 'C' and base2 == 'G') or (base1 == 'G' and base2 == 'C'):
        return valueCG
    else:
        return 0


def nussinov(seq, valueCG, valueAU):
    N = len(seq)
    scoreMatrix = np.zeros((N, N), dtype=int)
    for k in range(N):
        for i, j in zip(range(0, N - k), range(k+1, N)):
            max_k = (
                lambda i, j: max(
                    scoreMatrix[i, m] + scoreMatrix[m + 1, j] for m in range(i + 1, j)
                )
                if j-i > 1
                else 0
            )
            scoreMatrix[i][j] = max(
                scoreMatrix[i + 1, j],
                scoreMatrix[i, j - 1],
                scoreMatrix[i + 1, j - 1] + score(valueCG, valueAU, seq[i], seq[j]),
                max_k(i, j)
            )

    return scoreMatrix
def traceback(i, j, structure, DP, sequence, valueCG, valueAU):
    #in this case we've gone through the whole sequence. Nothing to do.
    if j <= i:
        return
    #if j is unpaired, there will be no change in score when we take it out, so we just recurse to the next index
    elif DP[i][j] == DP[i][j-1]:
        traceback(i, j-1, structure, DP, sequence, valueCG, valueAU)
    elif DP[i][j] == DP[i+1][j]:
        traceback(i+1, j, structure, DP, sequence, valueCG, valueAU)
    elif DP[i][j] == DP[i+1][j-1] + score(valueCG, valueAU, seq[i], seq[j]):
        structure.append((i, j))
        traceback(i + 1, j - 1, structure, DP, sequence, valueCG, valueAU)
    else:
        val, k = max([(DP[i][k] + DP[k + 1][j],k) for k in range(i + 1, j)])
        if DP[i, j] == val:
            traceback(i, k, structure, DP, sequence, valueCG, valueAU)
            traceback(k+1, j, structure, DP, sequence, valueCG, valueAU)

def write_structure(sequence, structure):
    dot_bracket = ["." for _ in range(len(sequence))]
    for s in structure:
        dot_bracket[min(s)] = "("
        dot_bracket[max(s)] = ")"
    return "".join(dot_bracket)
def parenthesingNussinov(seq, scoreMatrix, valueCG, valueAU,):
    struct = []
    traceback(0, len(seq)-1, struct, scoreMatrix, seq, valueCG, valueAU)
    return (write_structure(seq, struct), struct)
seq = "AGUCUGA"
print(nussinov(seq, 2, 1))
print(parenthesingNussinov(seq,nussinov(seq, 2, 1),2,1))
