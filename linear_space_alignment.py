#!/usr/bin/env python3

import sys

def LinearSpaceAlignment(match, mismatch, indel, v, w):
    def nw_score(s1, s2):
        """Compute the Needleman-Wunsch score matrix using linear space.
        Returns the last row of the score matrix."""
        previous = [j * (-indel) for j in range(len(s2) + 1)]
        for i, c1 in enumerate(s1, 1):
            current = [i * (-indel)]
            for j, c2 in enumerate(s2, 1):
                if c1 == c2:
                    score = previous[j-1] + match
                else:
                    score = previous[j-1] - mismatch
                delete = previous[j] - indel
                insert = current[j-1] - indel
                current.append(max(score, delete, insert))
            previous = current
        return previous

    def hirschberg(s1, s2):
        """Recursively compute the optimal alignment using Hirschberg's algorithm."""
        if len(s1) == 0:
            return '-' * len(s2), s2
        elif len(s2) == 0:
            return s1, '-' * len(s1)
        elif len(s1) == 1 or len(s2) == 1:
            return needleman_wunsch(s1, s2)
        else:
            i = len(s1) // 2
            score_l = nw_score(s1[:i], s2)
            score_r = nw_score(s1[i:][::-1], s2[::-1])
            max_j = 0
            max_score = float('-inf')
            for j in range(len(s2)+1):
                total = score_l[j] + score_r[len(s2)-j]
                # Change the condition to '>=', preferring the rightmost j in case of ties
                if total >= max_score:
                    max_score = total
                    max_j = j
            s1_left, s2_left = hirschberg(s1[:i], s2[:max_j])
            s1_right, s2_right = hirschberg(s1[i:], s2[max_j:])
            return s1_left + s1_right, s2_left + s2_right

    def needleman_wunsch(s1, s2):
        """Perform Needleman-Wunsch alignment for small strings."""
        m = len(s1)
        n = len(s2)
        score = [[0] * (n +1) for _ in range(m +1)]
        for i in range(m +1):
            score[i][0] = -indel * i
        for j in range(n +1):
            score[0][j] = -indel * j
        for i in range(1, m +1):
            for j in range(1, n +1):
                if s1[i-1] == s2[j-1]:
                    diag = score[i-1][j-1] + match
                else:
                    diag = score[i-1][j-1] - mismatch
                delete = score[i-1][j] - indel
                insert = score[i][j-1] - indel
                score[i][j] = max(diag, delete, insert)
        # Traceback
        align1 = ''
        align2 = ''
        i = m
        j = n
        while i >0 and j >0:
            current = score[i][j]
            if s1[i-1] == s2[j-1]:
                score_diag = score[i-1][j-1] + match
            else:
                score_diag = score[i-1][j-1] - mismatch
            if current == score_diag:
                align1 = s1[i-1] + align1
                align2 = s2[j-1] + align2
                i -=1
                j -=1
            elif current == score[i-1][j] - indel:
                align1 = s1[i-1] + align1
                align2 = '-' + align2
                i -=1
            else:
                align1 = '-' + align1
                align2 = s2[j-1] + align2
                j -=1
        while i >0:
            align1 = s1[i-1] + align1
            align2 = '-' + align2
            i -=1
        while j >0:
            align1 = '-' + align1
            align2 = s2[j-1] + align2
            j -=1
        return align1, align2

    aligned_v, aligned_w = hirschberg(v, w)
    # Calculate score
    score = 0
    for a, b in zip(aligned_v, aligned_w):
        if a == b and a != '-':
            score += match
        elif a == '-' or b == '-':
            score -= indel
        else:
            score -= mismatch
    return f"{score}\n{aligned_v}\n{aligned_w}"

if __name__ == "__main__":
    # Read input from stdin
    input_lines = sys.stdin.read().strip().split('\n')
    
    if len(input_lines) < 3:
        print("Invalid input format. Expected 3 lines.")
        sys.exit(1)
    
    # Parse the first line for match, mismatch, and indel penalties
    penalties = input_lines[0].split()
    if len(penalties) != 3:
        print("Invalid penalties line. Expected 3 integers.")
        sys.exit(1)
    try:
        match_reward = int(penalties[0])
        mismatch_penalty = int(penalties[1])
        indel_penalty = int(penalties[2])
    except ValueError:
        print("Penalties must be integers.")
        sys.exit(1)
    
    # Read the two sequences
    v = input_lines[1].strip()
    w = input_lines[2].strip()
    
    # Compute the alignment
    result = LinearSpaceAlignment(match_reward, mismatch_penalty, indel_penalty, v, w)
    
    # Print the result
    print(result)
