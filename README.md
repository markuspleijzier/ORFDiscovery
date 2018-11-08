# ORFDiscovery

ORFDiscovery is a miniprogram for searching and analysing Open Reading Frames (ORFs) within long non-coding RNAs (although it will work for any gene sequence (FASTA)). lncRNAs are stretches of RNA that, by definition, do not code for protein. However, this definition comes from genome sequencing studies and is slightly abitrary. lncRNAs are defined as small, mRNA-like stretches of RNA that are > 200 nucleutides long. This numerical length threshold exists because many sequences less than 200 nucleuotides do not code for protein, but rather have transcriptional/regulatory functions. However, this definition is based on older protein detection techniques like Western Blotting, where it is difficult to detect proteins < 100 amino acids long (300 nucleotides = 100 amino acids). _In silico_ analysis offers one avenue to test whether these lncRNAs are have coding potential, by examining the length of the polypeptide chains produced by the ORFs. ORFs are defined as sections of the sequences of interest where a start codon (AUG) is separated from the stop codon (UAA, UAG, UGA) by intermediary codons. Examining the length of the polypeptide is equivalent to measuring the length of the codon sequence (as 3 RNA nucleotides code for 1 amino acid).

The purpose of this code is not to show that lncRNAs all have coding potential - that requires wet-lab testing. The purpose is to serve as a platform for examining whether gene sequences of interest, described as lncRNAs, could potentially code for protein. Examining the distributions of ORF length within the sequence offers an indication of whether protein detection should be performed.

## Installation
### Prerequisites

ORFDiscovery is dependent upon the following packages:
    1. Pandas
    2. Numpy
    3. Matplotlib
    4. Collections
### Install

To install ORFDiscovery, simply run the following code in your terminal: 

```pip install git+git://github.com/markuspleijzier/ORFDiscovery@master```
    
## Tutorial

Now in a Python environment/notebook, we should import ORFDiscovery and our dependencies:

```python
import ORFDiscovery
import numpy as np
import pandas as ps
import matplotlib.pyplot as plt
import collections
```

Included in the install is an example gene sequence, called `example.txt`.

### 1. Extracting the Gene Sequence

Included in ORFDiscovery is a function, called `open_fasta`, which allows one to extract the genome sequence information of our example gene.


```python
gene = ORFDiscovery.open_fasta('example.txt')
```

The object `gene` has two attributes, which we can see if we run gene on its own: `info` and `data`.
While the `info` about the gene (its chromosomal location, which animal its from) is important, we are only interested in the actual gene sequence itself. To isolate this: 

```python
sequence = gene['data']
```

### 2. Getting Sense and Antisense Info

ORFDiscovery has two separate functions for getting the sense RNA and the antisense RNA from a gene. These are `get_sense_strand` and `get_antisense_strand`.

For this tutorial, let's focus on the antisense strand.

```python
antisense = ORFDiscovery.get_antisense_strand(sequence)
```

### 3. Codons

The function `get_codon` allows us to separate three letter codons based upon an index position. 

For example, if we want to get the first codon of the antisense strand, then we would run:

```python
ORFDiscovery.get_codon(antisense, position = 0)
```

* This function uses python based indexing: 0 is the first letter of the antisense sequence

Having just one codon isn't that interesting. However, the following two functions expand upon this `get_codon` function so that all of the codon sequences of the antisense strand are produced
_______________

The next two important functions are:
```
    1. rna_to_codon_sFrame()
    2. rna_to_codon_aFrame() 
```

__The `s` and `a` stand for Single and All, respectively__

The reason we have __single__ or __all__ is to increase the usability of ORFDiscovery. The frame position where gene translation begins will have very different consequences on the codon sequence and therefore amino acid sequence (polypeptide chain or primary structure) produced. For this reason we have `rna_to_codon_sFrame()` which takes an RNA sequence (sense or antisense) and finds the codon sequence produced based on a defined frame number.

* Here we do not have python based positional indexing, so the first frame is frame_position = 1

The `rna_to_codon_aFrame()` does not need a frame_position argument, as it will find all of the codon sequences on each of the three frame numbers (there are three frame numbers because the DNA code is a triplet code, where three DNA bases code for an amino acid). We can still access the individual codon sequences in `rna_to_codon_aFrame()` by subscripting. For example,

```
rna_to_codon_aFrame(antisense)[0]
```

is the codon sequence if we start from the first frame position. __Note that subscripting is python based indexing__


### 4. Amino Acid sequences (polypeptide chains)

To find the amino acid (AA) sequence of the codon in question, we use `codon_sequence_to_polypeptide()`. This function has one special argument, `frames` which will be __False__ if we used the single function (`rna_to_codon_sFrame`) but __True__ if we used the all function (`rna_to_codon_aFrame`).

However, there is one function that allows one to go straight from the gene sequence to polypeptides, which is built upon all the functions mentioned so far. This function is called
```
DNA_sequence_to_polypeptide()
```
and will take the original gene sequence and find all of the amino acid chains (all frame positions) for both sense and antisense strands, which can be accessed via subscripting. The `PRINT` argument for this function will just print out the sequences and indicate which frame number - this is useful for checking, but is only really beneficial for very small sequences.
```python
DNA_sequence_to_polypeptide(gene)[0]
```
is the amino acid sequences of all three frame positions for the sense strand. The first/ second/ third can be accessed via
```python
DAN_sequence_to_polypeptide(gene)[0][0]
DAN_sequence_to_polypeptide(gene)[0][1]
DAN_sequence_to_polypeptide(gene)[0][2]
```
whereas the antisense strand, for all three frame positions is:
```python
DAN_sequence_to_polypeptide(gene)[0][0]
DAN_sequence_to_polypeptide(gene)[1][1]
DAN_sequence_to_polypeptide(gene)[2][2]
```

### 5. Open Reading Frame Presence

Now that we have our amino acid sequences, we can check whether there are any ORFs. We use the `ORF_presence` function to return a list of Open Reading Frames that are present within that polypeptide sequence. As a reminder, an ORF is a sequence of amino acids that could potentially code for protein - all that is required of an ORF is that it has a start codon and a stop codon. As we have already converted the codon sequences into potential amino acid sequences, all we have to do is to look for sequences containing Methionine (AUG: __M__) and a stop codon (UAA, UAG, UGA: __\*__ ).

The `ORF_presence()` function will return a list of all of the possible ORFs. 

### 6. Open Reading Frame Analysis

The `ORF_lengths()` function looks at the ORFs produced by `ORF_presence` returns a `list` of the lengths, or it can do the counting already. This first is useful for creating histograms, but the second is useful for bar plots. For now, lets just create a list of the ORFs contained in the amino acid sequence from the second frame of the anitsense strand of the gene. 
```python
frame_2_antisense = ORFDiscovery.rna_to_codon_sFrame(antisense, frame_number = 2)
polypep = ORFDiscovery.codon_sequence_to_polypeptide(frame_2_antisense)
orf_presence = ORFDiscovery.ORF_presence(polypep)
LENGTHS= ORFDiscovery.lengths(orf_presence, counter = False)
```
__HOWEVER__

Before creating any charts, we can do our own analyses where the length of the ORF is set by us. To do this there is a `class` object called `ORFDiscovery.thresholds()` which allows one to define the thresholds. By manipulating classes (see [this link](https://www.learnpython.org/en/Classes_and_Objects) if unsure) we can set custom thresholds. But the default thresholds are 10, 25, 50, 100, 200, 300 AAs. So, to set up the threshold object, all one has to code is: 

```python
mythresh = ORFDiscovery.thresholds() 
mythresh.threshold_print()
```
The second line of code will tell you what the functions are. It is important to use the same object name for mythresh. 

_____________________

Now we can sort out the ORFs by their length, using the `measuringORFs()`, and get a dataframe of how many are contained within our thresholds

```python
Measuring_DF = ORFDiscovery.measuringORFs(LENGTHS, dataframe = True)
```

#### Plotting

Once we have the dataframe, now stored in `Measuring_DF` we can easily plot a bar chart using

```python
Measuring_DF.plot.bar()
```

![first-bar](https://github.com/markuspleijzier/ORFDiscovery/blob/master/images/Bar%20chart%20from%20measuringORFs.png)


If we want to look at the histogram of the length distribution then simply code

```python
fig, ax = plt.subplots()
n, bins, patches = plt.hist(LENGTHS, bins = 'auto')
plt.show()
```
![histogram](https://github.com/markuspleijzier/ORFDiscovery/blob/master/images/Histogram%20of%20ORF%20amino%20acid%20polypeptides.png)


Finally, if we want to see a bar chart where the x axis is continuous, showing the individual peaks for ORF length frequencies then we have to rerun our `ORF_lengths()` function but this time set the `counter` argument equal to True. 

```python
LENGTHS= ORFDiscovery.lengths(orf_presence, counter = True)
```
This returns a tuple of all of the ORF lengths found in the polypeptide sequence and their corresponding frequencies. (For more information, google collections.Counter python). We can then plot the bar chart with individual x-value peaks shown: 

```python
fig, ax = plt.subplots()
plt.bar(LENGTHS[0], LENGTHS[1])
plt.show()
```

![2bar](https://github.com/markuspleijzier/ORFDiscovery/blob/master/images/Bar%20chart%20of%20ORF%20length%20frequency%20sizes.png)

___________________________________


# Acknowledgements
This code was developed for Markus William Pleijzier's final year MSci project in the Neural Development, Plasticity & Repair Wing, Wolfson Institute of Biomedical Research, University College London, under Dr Huiliang Li. The title of the thesis which used this code was '_Expression of a novel lncRNA in Mouse Spinal Cord Oligodendrocytes_'.

# Verioning

* 01/11/2018 = Version 1.0

# License

GNU General Public License v3.0

