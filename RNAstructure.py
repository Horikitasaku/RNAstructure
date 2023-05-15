import src.rnastructure as rnastructure
from tqdm import tqdm
import src.util as util
import os
import numpy as np

def predictFromFasta(fasta_file, rnastructure_path='', temp_dir = 'temp', predict_structure = True, predict_pairing_probability = True, sequencer_noise=0):
    """
    Reads a fasta file and outputs a dict with reference, sequence and base-pairing prediction.

    Args:
        fasta_file (str): The path to the input fasta file.
        rnastructure_path (str): The path to the RNAstructure executable. Default is '' (i.e. RNAstructure is in the PATH).
        predict_structure (bool): Add structure prediction to the output.
        predict_pairing_probability (bool): Add pairing prediction to the output.
        sequencer_noise (float): The amount of sequencer noise to add to the base-pairing prediction. Default is 0. The noise follows a binomial distribution B(n=3000, p=sequencer_noise).
    
    Returns:
        data (dict): A dictionary with reference as key and sequence / main structure / pairing probability as value.
    
    Example:
        >>> predictFasta('testData/refs.fasta')
        Predicting RNA structures: 100%|███████████████████████████████████| 3/3 [00:01<00:00,  1.92seq/s]
        {'3042-O-flank_1=hp1-DB': 
            {
            'sequence':  'TTAAACCGGCCAACATACCGCATATGAGGATCACCCATATGCTCAAGATATTCGAAAGAATATCTTTCCACAGTCGAAAGACTGTGTCTCTCTCTTCCTTTTTCTCTTCCTCTTTCTCTTTCTCTTTCTCTTCTCTTCTGTATTACGAGTTCGCTACTCGTTCCTTTCGA', 
            'structure': '..............((((.(((((((.((....)))))))))..(((((((((....)))))))))..(((((((....))))))).....................................................)))).((((((.....)))))).........',
            'pairing': [0.0008763192815798796, 0.003282160575810705, 0.019867824327687585, 0.02829582753478791, 0.031266714341378524, 0.010288700618344545, 0.22755049833180016, 0.5741991237121613, 0.5018175337184233, 0.004588974386334777, 0.0027875607249180545, 0.010719083819803172, 0.0311800339547883, 0.026734132953685468, 0.5171140931132598, 0.526320420155328, 0.5302629606137425, 0.5973646221995781, 0.06787639095534494, 0.996380706717648, 0.9988325172245213, 0.9986708144405785, 0.9984979885954683, 0.9984931638770012, 0.9984703346038281, 0.9986003515173016, 0.0005489978550856028,
            },
        '3043-CC-flank_1=hp1-DB_2':
            ...
    """
    
    # make temp folder
    os.makedirs(temp_dir, exist_ok=True)
    
    output = util.fastaToDict(fasta_file)
    rna = rnastructure.RNAstructure(rnastructure_path, temp_dir)
    for ref in tqdm(output, desc='Predicting RNA structures'):
        
        # predict structure and pairing probability
        sequence = output[ref]['sequence']
        output[ref]['structure'] = rna.predictStructure(sequence)
        output[ref]['pairing'] = rna.predictPairingProbability(sequence)
        
        # add sequencer noise
        if sequencer_noise > 0:
            output[ref]['pairing'] = util.addBinomialNoise(signal=output[ref]['pairing'], n=3000, p=sequencer_noise)

    return output


def predictFromSequence(sequence, rnastructure_path='',  temp_dir = 'temp', predict_structure = True, predict_pairing_probability = True, constraints = [], dms = None, sequencer_noise=0, matrix = False):
    """
    Reads a RNA sequence and a list of constraints, and outputs a dict with reference, sequence and base-pairing prediction.

    Args:
        sequence (str): The RNA sequence.
        rnastructure_path (str): The path to the RNAstructure executable. Default is '' (i.e. RNAstructure is in the PATH).
        predict_structure (bool): Add structure prediction to the output.
        predict_pairing_probability (bool): Add pairing prediction to the output.
        constraints (list): A list of 1-based indexes of nucleotides that are constrained to be unpaired. Default is [].
        dms (list): A list of DMS reactivities, with the same length as the sequence. Used as an input to RNAstructure. Default is None.
        sequencer_noise (float): The amount of sequencer noise to add to the base-pairing prediction. Default is 0. The noise follows a binomial distribution B(n=3000, p=sequencer_noise).
        matrix (bool): make predict_pairing_probability a 2D matrix
    
    Returns:
        data (dict): A dictionary with sequence / main structure / pairing probability.
    
    Example:
        >>> predictSequence('TTAAACCGGCCAACATACCGCATATGAGGATCACCCATATGCTCAAGATATTCGAAAGAATATCTTTCCACAGTCGAAAGACTGTGTCTCTCTCTTCCTTTTTCTCTTCCTCTTTCTCTTTCTCTTTCTCTTCTCTTCTGTATTACGAGTTCGCTACTCGTTCCTTTCGA')
       {
        'sequence':  'TTAAACCGGCCAACATACCGCATATGAGGATCACCCATATGCTCAAGATATTCGAAAGAATATCTTTCCACAGTCGAAAGACTGTGTCTCTCTCTTCCTTTTTCTCTTCCTCTTTCTCTTTCTCTTTCTCTTCTCTTCTGTATTACGAGTTCGCTACTCGTTCCTTTCGA', 
        'structure': '..............((((.(((((((.((....)))))))))..(((((((((....)))))))))..(((((((....))))))).....................................................)))).((((((.....)))))).........',
        'pairing': [0.0008763192815798796, 0.003282160575810705, 0.019867824327687585, 0.02829582753478791, 0.031266714341378524, 0.010288700618344545, 0.22755049833180016, 0.5741991237121613, 0.5018175337184233, 0.004588974386334777, 0.0027875607249180545, 0.010719083819803172, 0.0311800339547883, 0.026734132953685468, 0.5171140931132598, 0.526320420155328, 0.5302629606137425, 0.5973646221995781, 0.06787639095534494, 0.996380706717648, 0.9988325172245213, 0.9986708144405785, 0.9984979885954683, 0.9984931638770012, 0.9984703346038281, 0.9986003515173016, 0.0005489978550856028,
        }

            ...
    """
    # make temp folder
    os.makedirs(temp_dir, exist_ok=True)

    # predict structure and pairing probability
    rna = rnastructure.RNAstructure(rnastructure_path, temp_dir)
    output = {'sequence': sequence}

    # make sequence lower when there are constraints => RNAstructure will use the constraints
    for c in constraints:
        sequence = sequence[:c-1] + sequence[c-1].lower() + sequence[c:]
        
    if predict_structure:
        output['structure'] = rna.predictStructure(sequence, dms=dms)
        
    if predict_pairing_probability:
        output['pairing'] = rna.predictPairingProbability(sequence, dms=dms, matrix=matrix)
        # add sequencer noise
        if sequencer_noise > 0:
            output['pairing'] = util.addBinomialNoise(signal=output['pairing'], n=3000, p=sequencer_noise)
        
    return output


if __name__=='__main__':
    print(predictFromSequence(
        'TTAAACCGGCCAACATACCGCATATGAGGATCACCCATATGCTCAAGATATTCGAAAGAATATCTTTCCACAGTCGAAAGACTGTGTCTCTCTCTTCCTTTTTCTCTTCCTCTTTCTCTTTCTCTTTCTCTTCTCTTCTGTATTACGAGTTCGCTACTCGTTCCTTTCGA',
        rnastructure_path='/Users/ymdt/src/RNAstructure/exe'))
    print(predictFromSequence(
        'TTAAACCGGCCAACATACCGCATATGAGGATCACCCATATGCTCAAGATATTCGAAAGAATATCTTTCCACAGTCGAAAGACTGTGTCTCTCTCTTCCTTTTTCTCTTCCTCTTTCTCTTTCTCTTTCTCTTCTCTTCTGTATTACGAGTTCGCTACTCGTTCCTTTCGA',
        dms = np.random.random(170),
        rnastructure_path='/Users/ymdt/src/RNAstructure/exe',
        matrix = True))
    print(predictFromFasta(
        rnastructure_path='/Users/ymdt/src/RNAstructure/exe',
        fasta_file='testData/refs.fasta')
    )