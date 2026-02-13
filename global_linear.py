from Bio import SeqIO
import numpy as np

def read_fasta(filename: str, num_sequences: int = 2):
    with open(filename) as handle:
        sequences = []
        for (i, record) in enumerate(SeqIO.parse(handle, "fasta")):
            if i < num_sequences:
                sequences.append(record.seq)
    return sequences

def read_cost_matrix_and_gap_cost(cost_filename: str):
    """
    Function to read a file with the cost matrix and the gapcost, 
    The input file must be a .txt file formated as follow:
    1 line:
    n g
    where n is the size of the alphabet and g is the linear gapcost

    next n lines:
    letter cost1 cost2 ... costn
    
    Where the letter is the i'th letter of the alphabet and each row defines the cost of a substitution
    of alphabet[i] to the rest of the alphabet
    Notice that there is just one space between values of the same line
    """
    with open(cost_filename) as f:
        lines = f.readlines()
        n, gapcost = map(int,lines[0].split(' '))
        alphabet = []
        cost_matrix = []
        for i in range(1, n+1):
            line = lines[i].split(' ')
            alphabet.append(line[0])
            cost_matrix.append(list(map(int, line[1:])))
        #print(n,gapcost)
        #print(alphabet)
        #print(cost_matrix)
    return gapcost, alphabet, cost_matrix
     
def check_sequences_validity(sequences, alphabet):
    """
    Function to check the validity of all the sequences
    Verify that they don't have characters outside of our alphabet

    (This is the lazy implementation, could be more efficient if we have the time to do so)
    """
    for sequence in sequences:
        for i, site in enumerate(sequence):
            if site not in alphabet:
                raise Exception(f"'{site}' of sequence {i} is not part of the defined alphabet {alphabet}")

def backtrack_from_score_matrix(i: int, j: int, sequences: list, score_matrix: list, cost_matrix: list, gapcost: int, alignment1: str = '', alignment2: str = ''):
    #print(i,j)
    if i > 0 and j > 0 and score_matrix[i][j] == score_matrix[i-1][j-1] + cost_matrix[alphabet.index(sequences[0][j-1])][alphabet.index(sequences[1][i-1])]:
        alignment1 += sequences[0][j-1]
        alignment2 += sequences[1][i-1]
        return backtrack_from_score_matrix(i-1, j-1, sequences, score_matrix, cost_matrix, gapcost, alignment1, alignment2)
    elif i > 0 and j >= 0 and score_matrix[i][j] == score_matrix[i-1][j] + gapcost:
        alignment1 += '-'
        alignment2 += sequences[1][i-1]
        return backtrack_from_score_matrix(i-1, j, sequences, score_matrix, cost_matrix, gapcost, alignment1, alignment2)
    elif i >= 0 and j > 0 and score_matrix[i][j] == score_matrix[i][j-1] + gapcost:
        alignment1 += sequences[0][j-1]
        alignment2 += '-'
        return backtrack_from_score_matrix(i, j-1, sequences, score_matrix, cost_matrix, gapcost, alignment1, alignment2)
    return alignment1[::-1], alignment2[::-1] #The alignments were filled from back to front, so we must invert the string of each


def global_linear(sequences: list, alphabet: list, cost_matrix: list, gap_cost: int, backtracking: bool = False):
    """
    Function to get the global pairwise alignment using linear gap cost
    It recieves:
        - sequences: a list of two sequences that we want to align
        - alphabet: a list of the current alphabet, aligned with the index of the cost matrix
        - cost_matrix: a matrix with the scores of each substitution
        - gap_cost: an integer representing the gap cost on the alignment
        - backtracking: a boolean variable to return an optimal alignment using backtracking, if desired 
    """
    #print(sequences)
    n = len(sequences[0])
    m = len(sequences[1])
    score_matrix = [list(range(n+1))] + [[(i+1)* gap_cost] + [0] * n for i in range(m)]
    score_matrix[0] = [i*gap_cost for i in score_matrix[0]] 
    #print(n,m)
    #Columns [1:n] seq 0 - j
    #Rows [1:m] seq 1 - i
    for i in range(1, m + 1):
        for j in range(1, n + 1):
                score_matrix[i][j] = min(score_matrix[i-1][j-1] + cost_matrix[alphabet.index(sequences[0][j-1])][alphabet.index(sequences[1][i-1])],
                                         score_matrix[i-1][j] + gap_cost,
                                         score_matrix[i][j-1] + gap_cost)
    #print(sequences[0])
    #print(sequences[1])
    #print(np.array(score_matrix)) #Convert to numpy just to view it pretty in terminal
    if backtracking:
        alignment = backtrack_from_score_matrix(m,n,sequences, score_matrix, cost_matrix, gap_cost)
        return score_matrix, alignment
    else:
        return score_matrix, ["", ""]

def write_alignment_in_fasta(alignment: list, filename: str):
    """
    Function to write the output alignment in a fasta file. Receives:
        - alignment: A list with the two sequences already aligned
        - filename: the name of the file where we will write our output
    (I'm sure that there is a more elegant way of doing this, but I don't know if I want to look at it)
    """
    print(filename)
    with open(filename, "w") as f:
        f.writelines([">seq1\n", alignment[0]+'\n', "\n", ">seq2\n", alignment[1]+'\n'])


tests_directory = 'tests/'
filename = "test4.fasta"
cost_matrix_filename = 'cost_matrix.txt'
gapcost, alphabet, cost_matrix = read_cost_matrix_and_gap_cost(tests_directory + cost_matrix_filename)
sequences = read_fasta(tests_directory + filename)
backtracking = True
check_sequences_validity(sequences, alphabet)
score_matrix, alignment = global_linear(sequences, alphabet, cost_matrix, gapcost, backtracking)
#print(np.array(score_matrix))
print(f"The cost of a global alignment is {score_matrix[-1][-1]}")

print(">seq1")
print(alignment[0])

print()
print(">seq2")
print(alignment[1])
write_alignment_in_fasta(alignment, tests_directory + "output_" + filename)
