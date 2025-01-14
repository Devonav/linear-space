README: Linear Space Alignment Tool
Overview
This Python program implements the Linear Space Alignment algorithm for solving the Global Alignment Problem. It uses Hirschberg's algorithm for memory-efficient sequence alignment and integrates Needleman-Wunsch for small subsequences. The program is designed to compute the optimal alignment and score between two input sequences, commonly used in bioinformatics for DNA, RNA, or protein sequence analysis.

Features
Linear Space Alignment (Hirschberg's Algorithm): Efficiently aligns long sequences using divide-and-conquer with reduced memory usage.
Needleman-Wunsch Alignment: Used for small subproblems, providing a detailed traceback of the alignment.
Customizable Penalties: Allows users to specify match reward, mismatch penalty, and indel (gap) penalty.
Dynamic Scoring: Computes the optimal alignment score based on the provided sequences and penalties.
Input Format
The program reads input from standard input (stdin) with the following format:

A single line containing three integers: match reward, mismatch penalty, and indel penalty.
A line containing the first sequence (v).
A line containing the second sequence (w).
Example:

Copy code
2 1 2
ACGT
ACCT
Output Format
The program outputs:

The alignment score.
The aligned version of the first sequence (v).
The aligned version of the second sequence (w).
Example Output:

r
Copy code
4
ACGT
AC-T
Usage
Running the Program
Save the script as linear_space_alignment.py.
Run the program using the command line:
bash
Copy code
./linear_space_alignment.py < input.txt
Replace input.txt with a file containing the formatted input as described above.
