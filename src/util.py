import numpy as np

def fastaToDict(fasta_file):
    """
    Reads a fasta file and returns a dictionary {sequence: reference}.
    
    Args:
        fasta_file (str): The path to the input fasta file.
    Returns:
        refSeq (dict): A dictionary with reference as key and sequence as value.
        
    """

    refSeq = {}

    with open(fasta_file, 'r') as handle:
        ref = 'placeholder'
        while True: # while not EOF
            ref, seq = handle.readline(), handle.readline()
            if ref == '' or seq == '':
                break
            assert ref[0] == '>', 'The reference should start with ">"'
            
            refSeq[ref[1:].strip()] = {'sequence': seq.strip().upper()}
    
    return refSeq


def addBinomialNoise(signal, n, p):
    """Add binomial noise to a signal.

    Args:
        signal (list): The signal to add noise to.
        n (int): The number of trials.
        p (float): The probability of success.        
    """
    
    return list(np.array(signal) + np.random.binomial(n, p, len(signal)) / n)


def get_pair_from_ct_line(line):
    while '  ' in line:
        line = line.replace('  ',' ')
    line = [l for l in line.split(' ') if l != '']
    return int(line[0]), int(line[4])

def ct2list(ctfile):
    pairs = []
    f = open(ctfile, 'r')
    L = int([l for l in f.readline().split(' ') if l != ''][0])
        
    for _ in range(L):
        l = f.readline()
        b1, b2 = get_pair_from_ct_line(l)
        if b2: 
            pairs.append((b1,b2))

    return pairs