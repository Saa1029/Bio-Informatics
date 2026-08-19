#!/usr/bin/env python3

def smith_waterman(seq1, seq2, match_score=2, mismatch_penalty=-1, gap_penalty=-2):
    """
    Smith-Waterman local alignment algorithm.

    Parameters:
        seq1 (str): First sequence
        seq2 (str): Second sequence
        match_score (int): Score for a match
        mismatch_penalty (int): Penalty for a mismatch
        gap_penalty (int): Penalty for a gap

    Returns:
        tuple: best local alignment of seq1 and seq2, alignment score
    """

    rows = len(seq1) + 1
    cols = len(seq2) + 1

    # Score matrix
    score_matrix = [[0 for _ in range(cols)] for _ in range(rows)]

    # Traceback matrix
    # Values: "D" = diagonal, "U" = up, "L" = left, "0" = stop
    traceback_matrix = [["0" for _ in range(cols)] for _ in range(rows)]

    max_score = 0
    max_position = (0, 0)

    # Fill matrices
    for i in range(1, rows):
        for j in range(1, cols):

            if seq1[i - 1] == seq2[j - 1]:
                diagonal_score = score_matrix[i - 1][j - 1] + match_score
            else:
                diagonal_score = score_matrix[i - 1][j - 1] + mismatch_penalty

            up_score = score_matrix[i - 1][j] + gap_penalty
            left_score = score_matrix[i][j - 1] + gap_penalty

            best_score = max(0, diagonal_score, up_score, left_score)

            score_matrix[i][j] = best_score

            if best_score == 0:
                traceback_matrix[i][j] = "0"
            elif best_score == diagonal_score:
                traceback_matrix[i][j] = "D"
            elif best_score == up_score:
                traceback_matrix[i][j] = "U"
            else:
                traceback_matrix[i][j] = "L"

            if best_score > max_score:
                max_score = best_score
                max_position = (i, j)

    # Traceback from highest-scoring cell
    aligned_seq1 = []
    aligned_seq2 = []

    i, j = max_position

    while traceback_matrix[i][j] != "0":
        direction = traceback_matrix[i][j]

        if direction == "D":
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append(seq2[j - 1])
            i -= 1
            j -= 1

        elif direction == "U":
            aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append("-")
            i -= 1

        elif direction == "L":
            aligned_seq1.append("-")
            aligned_seq2.append(seq2[j - 1])
            j -= 1

    aligned_seq1 = "".join(reversed(aligned_seq1))
    aligned_seq2 = "".join(reversed(aligned_seq2))

    return aligned_seq1, aligned_seq2, max_score, score_matrix


def print_alignment(aligned_seq1, aligned_seq2):
    """
    Print alignment with match indicators.
    """

    match_line = []

    for a, b in zip(aligned_seq1, aligned_seq2):
        if a == b:
            match_line.append("|")
        elif a == "-" or b == "-":
            match_line.append(" ")
        else:
            match_line.append(".")

    print(aligned_seq1)
    print("".join(match_line))
    print(aligned_seq2)


if __name__ == "__main__":

    seq1 = "ACACACTGA"
    seq2 = "AGCACACA"

    seq1 = "GGATCAGTACGTTACCGGAT"
    seq2 = "TTACGATACCGTTAAGG"

    aligned_seq1, aligned_seq2, score, score_matrix = smith_waterman(
        seq1,
        seq2,
        match_score=2,
        mismatch_penalty=-1,
        gap_penalty=-2
    )

    print("Sequence 1:", seq1)
    print("Sequence 2:", seq2)
    print("\nBest Local Alignment:")
    print_alignment(aligned_seq1, aligned_seq2)

    print("\nAlignment Score:", score)
