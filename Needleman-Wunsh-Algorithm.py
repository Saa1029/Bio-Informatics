#!/usr/bin/env python
# coding: utf-8

# In[1]:


def needleman_wunsch(seq1, seq2, match_score=1, mismatch_penalty=-1, gap_penalty=-2):
    m = len(seq1)
    n = len(seq2)

    score_matrix = [[0 for _ in range(n + 1)] for _ in range(m + 1)]

    for i in range(1, m + 1):
        score_matrix[i][0] = score_matrix[i - 1][0] + gap_penalty

    for j in range(1, n + 1):
        score_matrix[0][j] = score_matrix[0][j - 1] + gap_penalty

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            if seq1[i - 1] == seq2[j - 1]:
                diagonal_score = score_matrix[i - 1][j - 1] + match_score
            else:
                diagonal_score = score_matrix[i - 1][j - 1] + mismatch_penalty

            up_score = score_matrix[i - 1][j] + gap_penalty
            left_score = score_matrix[i][j - 1] + gap_penalty

            score_matrix[i][j] = max(diagonal_score, up_score, left_score)

    aligned_seq1 = ""
    aligned_seq2 = ""

    i = m
    j = n

    while i > 0 and j > 0:
        current_score = score_matrix[i][j]

        if seq1[i - 1] == seq2[j - 1]:
            score = match_score
        else:
            score = mismatch_penalty

        if current_score == score_matrix[i - 1][j - 1] + score:
            aligned_seq1 = seq1[i - 1] + aligned_seq1
            aligned_seq2 = seq2[j - 1] + aligned_seq2
            i -= 1
            j -= 1

        elif current_score == score_matrix[i - 1][j] + gap_penalty:
            aligned_seq1 = seq1[i - 1] + aligned_seq1
            aligned_seq2 = "-" + aligned_seq2
            i -= 1

        else:
            aligned_seq1 = "-" + aligned_seq1
            aligned_seq2 = seq2[j - 1] + aligned_seq2
            j -= 1

    while i > 0:
        aligned_seq1 = seq1[i - 1] + aligned_seq1
        aligned_seq2 = "-" + aligned_seq2
        i -= 1

    while j > 0:
        aligned_seq1 = "-" + aligned_seq1
        aligned_seq2 = seq2[j - 1] + aligned_seq2
        j -= 1

    return aligned_seq1, aligned_seq2, score_matrix[m][n], score_matrix


def print_matrix(matrix, seq1, seq2):
    print("      ", end="")
    for char in seq2:
        print(f"{char:4}", end="")
    print()

    for i in range(len(matrix)):
        if i == 0:
            print(" ", end=" ")
        else:
            print(seq1[i - 1], end=" ")

        for j in range(len(matrix[0])):
            print(f"{matrix[i][j]:4}", end="")
        print()


seq1 = "GATTACA"
seq2 = "GCATGCU"

aligned_seq1, aligned_seq2, score, matrix = needleman_wunsch(seq1, seq2)

print("Original sequence 1:", seq1)
print("Original sequence 2:", seq2)
print()
print("Aligned sequence 1:", aligned_seq1)
print("Aligned sequence 2:", aligned_seq2)
print("Alignment score:", score)
print()
print("Scoring matrix:")
print_matrix(matrix, seq1, seq2)


# In[ ]:




