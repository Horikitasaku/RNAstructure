# RNAstructure

A python wrapper for [RNAstructure (Matthews lab @Rochester)](https://rna.urmc.rochester.edu/RNAstructure.html).

The function `predictFasta` reads a fasta file and outputs a dict with the reference, the sequence, the structure and a base-pairing (bp) probability prediction array.

The function `predictSequence` reads a fasta file and a list of constraints, and outputs a dict with the reference, the sequence, the structure and a base-pairing (bp) probability prediction array.

A base-pairing (bp) probability prediction array that can be used as a synthetic normalized DMS signal.


## Key functions

```Python
def predictFromFasta(fasta_file, rnastructure_path='', predict_structure = True, predict_pairing_probability = True, sequencer_noise=0.001):
    """
    Reads a fasta file and outputs a dict with reference, sequence and base-pairing prediction.

    Args:
        fasta_file (str): The path to the input fasta file.
        rnastructure_path (str): The path to the RNAstructure executable. Default is '' (i.e. RNAstructure is in the PATH).
        predict_structure (bool): Add structure prediction to the output.
        predict_pairing_probability (bool): Add pairing prediction to the output.
        sequencer_noise (float): The amount of sequencer noise to add to the base-pairing prediction. Default is 0.001. The noise follows a binomial distribution B(n=3000, p=sequencer_noise).
    
    Returns:
        data (dict): A dictionary with reference as key and sequence / main structure / pairing probability as value.
    
    Example:
        >>> predictFasta('testData/refs.fasta')
        Predicting RNA structures: 100%|███████████████████████████████████| 3/3 [00:01<00:00,  1.92seq/s]
        {'3042-O-flank_1=hp1-DB': 
            {
            'sequence': 'TTAAACCGGCCAACATACCGCATATGAGGATCACCCATATGCTCAAGATATTCGAAAGAATATCTTTCCACAGTCGAAAGACTGTGTCTCTCTCTTCCTTTTTCTCTTCCTCTTTCTCTTTCTCTTTCTCTTCTCTTCTGTATTACGAGTTCGCTACTCGTTCCTTTCGA',        
            'structure': '..............((((.(((((((.((....)))))))))..(((((((((....)))))))))..(((((((....))))))).....................................................)))).((((((.....)))))).........',
            'pairing': [0.0008763192815798796, 0.003282160575810705, 0.019867824327687585, 0.02829582753478791, 0.031266714341378524, 0.010288700618344545, 0.22755049833180016, 0.5741991237121613, 0.5018175337184233, 0.004588974386334777, 0.0027875607249180545, 0.010719083819803172, 0.0311800339547883, 0.026734132953685468, 0.5171140931132598, 0.526320420155328, 0.5302629606137425, 0.5973646221995781, 0.06787639095534494, 0.996380706717648, 0.9988325172245213, 0.9986708144405785, 0.9984979885954683, 0.9984931638770012, 0.9984703346038281, 0.9986003515173016, 0.0005489978550856028,
            },
        '3043-CC-flank_1=hp1-DB_2':
            ...
    """

def predictFromSequence(sequence, rnastructure_path='', predict_structure = True, predict_pairing_probability = True, constraints = [], sequencer_noise=0.001):
    """
    Reads a RNA sequence and a list of constraints, and outputs a dict with reference, sequence and base-pairing prediction.

    Args:
        sequence (str): The RNA sequence.
        rnastructure_path (str): The path to the RNAstructure executable. Default is '' (i.e. RNAstructure is in the PATH).
        predict_structure (bool): Add structure prediction to the output.
        predict_pairing_probability (bool): Add pairing prediction to the output.
        constraints (list): A list of 1-based indexes of nucleotides that are constrained to be unpaired. Default is [].
        sequencer_noise (float): The amount of sequencer noise to add to the base-pairing prediction. Default is 0.001. The noise follows a binomial distribution B(n=3000, p=sequencer_noise).
    
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
```


## Installation

### Install RNAstructure

If RNAstructure is already installed, skip this step.

**Linux and MacOS:**

```
git clone https://github.com/rouskinlab/RNAstructure
cd RNAstructure
unzip RNAstructure.zip
cd RNAstructure
make all
export PATH=$PATH:$(pwd)/exe
```

If you want to make RNAstructure available in all your terminal sessions, add the line `export PATH=$PATH:$(pwd)/exe` to your `.bashrc` file.

```
echo "export PATH=$PATH:$(pwd)/exe" >> ~/.bashrc
```

**Windows:**

```
git clone https://github.com/rouskinlab/RNAstructure
cd RNAstructure
Expand-Archive RNAstructure.zip
cd RNAstructure
nmake /f Makefile all
setx PATH "%PATH%;%cd%\exe"
```

If you want to make RNAstructure available in all your terminal sessions, you can add the following line to your System Environment Variables:

1. Press the Windows key + R to open the Run box.

2. Type sysdm.cpl and press Enter to open the System Properties window.

3. Click on the "Advanced" tab and click on the "Environment Variables" button.

4. Under "System Variables", click on the "New" button and enter the following:

5. Variable name: RNASTRUCTURE_PATH

6. Variable value: <path to RNAstructure>/exe

Note: Replace <path to RNAstructure> with the actual path to the RNAstructure folder.

7. Click on "OK" to close all windows.
