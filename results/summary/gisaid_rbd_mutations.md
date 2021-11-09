# Align spike sequences in GISAID and count RBD mutations
This Python Jupyter notebook reads in a file of all spike sequences from GISAID, parses for "high-quality" sequences, builds a RBD alignment, and then makes a file that gives the counts of each mutation at each site.

## Set up analysis
Import Python modules:


```python
import io
import lzma
import os
import re
import subprocess

from Bio.Data.IUPACData import protein_letters
import Bio.SeqIO

from IPython.display import display, HTML

import matplotlib.pyplot as plt

import pandas as pd

from plotnine import *

import yaml
```

Read the configuration file:


```python
with open('config.yaml') as f:
    config = yaml.safe_load(f)
```

Create output directory:


```python
os.makedirs(config['gisaid_mutations_dir'], exist_ok=True)
```

## Parse full-length human human spike sequences

Read the spikes from the file downloaded from GISAID:


```python
print(f"Reading GISAID spikes in {config['gisaid_spikes']}")
# file is `xz` compressed
with lzma.open(config['gisaid_spikes'], 'rt') as f:
    spikes = list(Bio.SeqIO.parse(f, 'fasta'))   
print(f"Read {len(spikes)} spike sequences.")
```

    Reading GISAID spikes in data/spikeprot1030.tar.tar.xz
    Read 4607422 spike sequences.


Make a data frame that has the BioPython SeqRecord, length, host, and geographic location (country) for each spike.
Also determine whether sequences have ambiguous amino acids or are all valid amino acids:


```python
spikes_df = (
    pd.DataFrame({'seqrecord': spikes})
    .assign(description=lambda x: x['seqrecord'].map(lambda rec: rec.description),
            country=lambda x: x['description'].str.split('|').str[-1],
            host=lambda x: x['description'].str.split('|').str[6].str.strip(),
            length=lambda x: x['seqrecord'].map(len),
            n_ambiguous=lambda x: x['seqrecord'].map(lambda rec: rec.seq.count('X') + rec.seq.count('x')),
            )
    )
```

Show number of sequences from different hosts, then keep only human ones:


```python
display(HTML(
    spikes_df
    .groupby('host')
    .aggregate(n_sequences=pd.NamedAgg('seqrecord', 'count'))
    .sort_values('n_sequences', ascending=False)
    .to_html()
    ))

spikes_df = spikes_df.query('host == "Human"')
```


<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>n_sequences</th>
    </tr>
    <tr>
      <th>host</th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>Human</th>
      <td>4604437</td>
    </tr>
    <tr>
      <th>Environment</th>
      <td>1696</td>
    </tr>
    <tr>
      <th>Neovison vison</th>
      <td>961</td>
    </tr>
    <tr>
      <th>Felis catus</th>
      <td>83</td>
    </tr>
    <tr>
      <th>Canis lupus familiaris</th>
      <td>44</td>
    </tr>
    <tr>
      <th>Panthera leo</th>
      <td>37</td>
    </tr>
    <tr>
      <th>Mustela lutreola</th>
      <td>23</td>
    </tr>
    <tr>
      <th>Manis javanica</th>
      <td>19</td>
    </tr>
    <tr>
      <th>Odocoileus virginianus</th>
      <td>14</td>
    </tr>
    <tr>
      <th>Neovision vision</th>
      <td>14</td>
    </tr>
    <tr>
      <th>Panthera tigris jacksoni</th>
      <td>13</td>
    </tr>
    <tr>
      <th>Environmental</th>
      <td>10</td>
    </tr>
    <tr>
      <th>Rhinolophus malayanus</th>
      <td>7</td>
    </tr>
    <tr>
      <th>unknown</th>
      <td>6</td>
    </tr>
    <tr>
      <th>Aonyx cinereus</th>
      <td>5</td>
    </tr>
    <tr>
      <th>P1 culture</th>
      <td>5</td>
    </tr>
    <tr>
      <th>Snow Leopard</th>
      <td>4</td>
    </tr>
    <tr>
      <th>Gorilla gorilla gorilla</th>
      <td>3</td>
    </tr>
    <tr>
      <th>Panthera uncia</th>
      <td>3</td>
    </tr>
    <tr>
      <th>Panthera tigris</th>
      <td>3</td>
    </tr>
    <tr>
      <th>P2 culture</th>
      <td>3</td>
    </tr>
    <tr>
      <th>Mus musculus</th>
      <td>3</td>
    </tr>
    <tr>
      <th>Tiger</th>
      <td>2</td>
    </tr>
    <tr>
      <th>Dog</th>
      <td>2</td>
    </tr>
    <tr>
      <th>Mesocricetus auratus</th>
      <td>2</td>
    </tr>
    <tr>
      <th>Control</th>
      <td>2</td>
    </tr>
    <tr>
      <th>Rhinolophus stheno</th>
      <td>2</td>
    </tr>
    <tr>
      <th>Mustela putorius furo</th>
      <td>2</td>
    </tr>
    <tr>
      <th>Rhinolophus shameli</th>
      <td>2</td>
    </tr>
    <tr>
      <th>Rhinolophus pusillus</th>
      <td>2</td>
    </tr>
    <tr>
      <th>North America / USA / Louisiana / New Orleands</th>
      <td>1</td>
    </tr>
    <tr>
      <th>Chlorocebus sabaeus</th>
      <td>1</td>
    </tr>
    <tr>
      <th>Rhinolophus sinicus</th>
      <td>1</td>
    </tr>
    <tr>
      <th>Rhinolophus marshalli</th>
      <td>1</td>
    </tr>
    <tr>
      <th>Rhinolophus affinis</th>
      <td>1</td>
    </tr>
    <tr>
      <th>Rhinolophus bat</th>
      <td>1</td>
    </tr>
    <tr>
      <th>Host</th>
      <td>1</td>
    </tr>
    <tr>
      <th>Puma concolor</th>
      <td>1</td>
    </tr>
    <tr>
      <th>Prionailurus bengalensis euptilurus</th>
      <td>1</td>
    </tr>
    <tr>
      <th>Mus musculus (BALB/c mice)</th>
      <td>1</td>
    </tr>
    <tr>
      <th>Manis pentadactyla</th>
      <td>1</td>
    </tr>
    <tr>
      <th>P3 culture</th>
      <td>1</td>
    </tr>
    <tr>
      <th>Panthera tigris sondaica</th>
      <td>1</td>
    </tr>
  </tbody>
</table>


Plot distribution of lengths and only keep sequences that are full-length (or near full-length) spikes:


```python
print('Distribution of length for all spikes:')
p = (ggplot(spikes_df) +
     aes('length') +
     geom_bar() +
     ylab('number of sequences') +
     theme(figure_size=(10, 2))
     )
fig = p.draw()
display(fig)
plt.close(fig)

min_length, max_length = 1260, 1276
print(f"\nOnly keeping spikes with lengths between {min_length} and {max_length}")
spikes_df = (
    spikes_df
    .assign(valid_length=lambda x: (min_length <= x['length']) & (x['length'] <= max_length))
    )

print('Here are number of sequences with valid and invalid lengths:')
display(HTML(spikes_df
             .groupby('valid_length')
             .aggregate(n_sequences=pd.NamedAgg('seqrecord', 'count'))
             .to_html()
             ))

print('\nDistribution of lengths for sequences with valid and invalid lengths; '
      'dotted red lines delimit valid lengths:')
p = (ggplot(spikes_df
            .assign(valid_length=lambda x: x['valid_length'].map({True: 'valid length',
                                                                  False: 'invalid length'}))
            ) +
     aes('length') +
     geom_bar() +
     ylab('number of sequences') +
     theme(figure_size=(10, 2), subplots_adjust={'wspace': 0.2}) +
     facet_wrap('~ valid_length', scales='free') +
     geom_vline(xintercept=min_length - 0.5, color='red', linetype='dotted') +
     geom_vline(xintercept=max_length + 0.5, color='red', linetype='dotted')
     )
fig = p.draw()
display(fig)
plt.close(fig)

spikes_df = spikes_df.query('valid_length')
```

    Distribution of length for all spikes:



    
![png](gisaid_rbd_mutations_files/gisaid_rbd_mutations_14_1.png)
    


    
    Only keeping spikes with lengths between 1260 and 1276
    Here are number of sequences with valid and invalid lengths:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>n_sequences</th>
    </tr>
    <tr>
      <th>valid_length</th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>False</th>
      <td>83275</td>
    </tr>
    <tr>
      <th>True</th>
      <td>4521162</td>
    </tr>
  </tbody>
</table>


    
    Distribution of lengths for sequences with valid and invalid lengths; dotted red lines delimit valid lengths:



    
![png](gisaid_rbd_mutations_files/gisaid_rbd_mutations_14_5.png)
    


Finally, we get rid of spikes with **lots** of ambiguous residues as they may hinder the alignment below.
We will then do more detailed filtering for ambiguous residues just in the RBD region after alignment:


```python
max_ambiguous = 100
print(f"Filtering sequences with > {max_ambiguous} ambiguous residues")
spikes_df = (
    spikes_df
    .assign(excess_ambiguous=lambda x: x['n_ambiguous'] > max_ambiguous)
    )
display(HTML(
    spikes_df
    .groupby('excess_ambiguous')
    .aggregate(n_sequences=pd.NamedAgg('seqrecord', 'count'))
    .to_html()
    ))
```

    Filtering sequences with > 100 ambiguous residues



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>n_sequences</th>
    </tr>
    <tr>
      <th>excess_ambiguous</th>
      <th></th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>False</th>
      <td>4173868</td>
    </tr>
    <tr>
      <th>True</th>
      <td>347294</td>
    </tr>
  </tbody>
</table>


## Align the RBD region of the spikes
We now align the RBD regions of the spikes.
We do this **before** we filter sequences with too many ambiguous residues so that we can do that filtering just on the RBD region.

We align with `mafft` using the `--addfragments` and `--keeplength` options (see [here](https://mafft.cbrc.jp/alignment/software/closelyrelatedviralgenomes.html) and [here](https://mafft.cbrc.jp/alignment/software/addsequences.html)) to align relative to a reference that is just the RBD; these options also clip the sequences relative to the reference.
Note that these options make sense if the following conditions are met:
  1. Sequences are all very similar.
  2. We are not worried about insertions.
For now, both of these appear to be true, but this choice should be kept in mind if there is a lot more divergence or insertions.

We align relative to the reference that is the wildtype sequence used for the experiments:


```python
print(f"Reading reference nucleotide sequence in {config['wildtype_sequence']}")
refseq = Bio.SeqIO.read(config['wildtype_sequence'], 'fasta')

refprotfile = os.path.join(config['gisaid_mutations_dir'], 'reference_RBD.fasta')
print(f"Writing protein translation of reference sequence to {refprotfile}")
refseq.seq = refseq.seq.translate()
_ = Bio.SeqIO.write(refseq, refprotfile, 'fasta')
```

    Reading reference nucleotide sequence in data/wildtype_sequence.fasta
    Writing protein translation of reference sequence to results/GISAID_mutations/reference_RBD.fasta


Write all the other spikes to a file:


```python
spikes_file = os.path.join(config['gisaid_mutations_dir'],
                           'human_full-length_spikes.fasta')
print(f"Writing the spikes to {spikes_file}")
_ = Bio.SeqIO.write(spikes_df['seqrecord'].tolist(), spikes_file, 'fasta')
```

    Writing the spikes to results/GISAID_mutations/human_full-length_spikes.fasta


Now make the alignment.
Note that we use multiple threads to speed things up, and also align the spikes in chunks.
The reason that we have to the chunkwise alignment is that some unclear `mafft` error was arising if we tried to align them all at once:


```python
chunksize = 50000

aligned_rbds = []

for i in range(0, len(spikes_df), chunksize):
    spikes_file = os.path.join(config['gisaid_mutations_dir'],
                               f"human_full-length_spikes_{i + 1}-to-{i + chunksize}.fasta")
    print(f"Writing spikes {i + 1} to {i + chunksize} to {spikes_file}")
    _ = Bio.SeqIO.write(spikes_df['seqrecord'].tolist()[i: i + chunksize], spikes_file, 'fasta')
    print('Now aligning these sequences...')
    cmds = ['mafft', '--auto', '--thread', str(config['max_cpus']),
            '--keeplength', '--addfragments', spikes_file, refprotfile]
    res = subprocess.run(cmds, capture_output=True)
    if res.returncode:
        raise RuntimeError(f"Error in alignment:\n{res.stderr}")
    else:
        print('Alignment complete.\n')
        with io.StringIO(res.stdout.decode('utf-8')) as f:
            iseqs = list(Bio.SeqIO.parse(f, 'fasta'))
            # remove reference sequence, which should be first in file
            assert iseqs[0].seq == refseq.seq and iseqs[0].description == refseq.description
            iseqs = iseqs[1:]
            assert len(iseqs) == min(chunksize, len(spikes_df) - i)
            aligned_rbds += iseqs
            
assert len(aligned_rbds) == len(spikes_df)
```

    Writing spikes 1 to 50000 to results/GISAID_mutations/human_full-length_spikes_1-to-50000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 50001 to 100000 to results/GISAID_mutations/human_full-length_spikes_50001-to-100000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 100001 to 150000 to results/GISAID_mutations/human_full-length_spikes_100001-to-150000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 150001 to 200000 to results/GISAID_mutations/human_full-length_spikes_150001-to-200000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 200001 to 250000 to results/GISAID_mutations/human_full-length_spikes_200001-to-250000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 250001 to 300000 to results/GISAID_mutations/human_full-length_spikes_250001-to-300000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 300001 to 350000 to results/GISAID_mutations/human_full-length_spikes_300001-to-350000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 350001 to 400000 to results/GISAID_mutations/human_full-length_spikes_350001-to-400000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 400001 to 450000 to results/GISAID_mutations/human_full-length_spikes_400001-to-450000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 450001 to 500000 to results/GISAID_mutations/human_full-length_spikes_450001-to-500000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 500001 to 550000 to results/GISAID_mutations/human_full-length_spikes_500001-to-550000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 550001 to 600000 to results/GISAID_mutations/human_full-length_spikes_550001-to-600000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 600001 to 650000 to results/GISAID_mutations/human_full-length_spikes_600001-to-650000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 650001 to 700000 to results/GISAID_mutations/human_full-length_spikes_650001-to-700000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 700001 to 750000 to results/GISAID_mutations/human_full-length_spikes_700001-to-750000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 750001 to 800000 to results/GISAID_mutations/human_full-length_spikes_750001-to-800000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 800001 to 850000 to results/GISAID_mutations/human_full-length_spikes_800001-to-850000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 850001 to 900000 to results/GISAID_mutations/human_full-length_spikes_850001-to-900000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 900001 to 950000 to results/GISAID_mutations/human_full-length_spikes_900001-to-950000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 950001 to 1000000 to results/GISAID_mutations/human_full-length_spikes_950001-to-1000000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 1000001 to 1050000 to results/GISAID_mutations/human_full-length_spikes_1000001-to-1050000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 1050001 to 1100000 to results/GISAID_mutations/human_full-length_spikes_1050001-to-1100000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 1100001 to 1150000 to results/GISAID_mutations/human_full-length_spikes_1100001-to-1150000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 1150001 to 1200000 to results/GISAID_mutations/human_full-length_spikes_1150001-to-1200000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 1200001 to 1250000 to results/GISAID_mutations/human_full-length_spikes_1200001-to-1250000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 1250001 to 1300000 to results/GISAID_mutations/human_full-length_spikes_1250001-to-1300000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 1300001 to 1350000 to results/GISAID_mutations/human_full-length_spikes_1300001-to-1350000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 1350001 to 1400000 to results/GISAID_mutations/human_full-length_spikes_1350001-to-1400000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 1400001 to 1450000 to results/GISAID_mutations/human_full-length_spikes_1400001-to-1450000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 1450001 to 1500000 to results/GISAID_mutations/human_full-length_spikes_1450001-to-1500000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 1500001 to 1550000 to results/GISAID_mutations/human_full-length_spikes_1500001-to-1550000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 1550001 to 1600000 to results/GISAID_mutations/human_full-length_spikes_1550001-to-1600000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 1600001 to 1650000 to results/GISAID_mutations/human_full-length_spikes_1600001-to-1650000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 1650001 to 1700000 to results/GISAID_mutations/human_full-length_spikes_1650001-to-1700000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 1700001 to 1750000 to results/GISAID_mutations/human_full-length_spikes_1700001-to-1750000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 1750001 to 1800000 to results/GISAID_mutations/human_full-length_spikes_1750001-to-1800000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 1800001 to 1850000 to results/GISAID_mutations/human_full-length_spikes_1800001-to-1850000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 1850001 to 1900000 to results/GISAID_mutations/human_full-length_spikes_1850001-to-1900000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 1900001 to 1950000 to results/GISAID_mutations/human_full-length_spikes_1900001-to-1950000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 1950001 to 2000000 to results/GISAID_mutations/human_full-length_spikes_1950001-to-2000000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 2000001 to 2050000 to results/GISAID_mutations/human_full-length_spikes_2000001-to-2050000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 2050001 to 2100000 to results/GISAID_mutations/human_full-length_spikes_2050001-to-2100000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 2100001 to 2150000 to results/GISAID_mutations/human_full-length_spikes_2100001-to-2150000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 2150001 to 2200000 to results/GISAID_mutations/human_full-length_spikes_2150001-to-2200000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 2200001 to 2250000 to results/GISAID_mutations/human_full-length_spikes_2200001-to-2250000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 2250001 to 2300000 to results/GISAID_mutations/human_full-length_spikes_2250001-to-2300000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 2300001 to 2350000 to results/GISAID_mutations/human_full-length_spikes_2300001-to-2350000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 2350001 to 2400000 to results/GISAID_mutations/human_full-length_spikes_2350001-to-2400000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 2400001 to 2450000 to results/GISAID_mutations/human_full-length_spikes_2400001-to-2450000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 2450001 to 2500000 to results/GISAID_mutations/human_full-length_spikes_2450001-to-2500000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 2500001 to 2550000 to results/GISAID_mutations/human_full-length_spikes_2500001-to-2550000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 2550001 to 2600000 to results/GISAID_mutations/human_full-length_spikes_2550001-to-2600000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 2600001 to 2650000 to results/GISAID_mutations/human_full-length_spikes_2600001-to-2650000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 2650001 to 2700000 to results/GISAID_mutations/human_full-length_spikes_2650001-to-2700000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 2700001 to 2750000 to results/GISAID_mutations/human_full-length_spikes_2700001-to-2750000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 2750001 to 2800000 to results/GISAID_mutations/human_full-length_spikes_2750001-to-2800000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 2800001 to 2850000 to results/GISAID_mutations/human_full-length_spikes_2800001-to-2850000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 2850001 to 2900000 to results/GISAID_mutations/human_full-length_spikes_2850001-to-2900000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 2900001 to 2950000 to results/GISAID_mutations/human_full-length_spikes_2900001-to-2950000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 2950001 to 3000000 to results/GISAID_mutations/human_full-length_spikes_2950001-to-3000000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 3000001 to 3050000 to results/GISAID_mutations/human_full-length_spikes_3000001-to-3050000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 3050001 to 3100000 to results/GISAID_mutations/human_full-length_spikes_3050001-to-3100000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 3100001 to 3150000 to results/GISAID_mutations/human_full-length_spikes_3100001-to-3150000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 3150001 to 3200000 to results/GISAID_mutations/human_full-length_spikes_3150001-to-3200000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 3200001 to 3250000 to results/GISAID_mutations/human_full-length_spikes_3200001-to-3250000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 3250001 to 3300000 to results/GISAID_mutations/human_full-length_spikes_3250001-to-3300000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 3300001 to 3350000 to results/GISAID_mutations/human_full-length_spikes_3300001-to-3350000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 3350001 to 3400000 to results/GISAID_mutations/human_full-length_spikes_3350001-to-3400000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 3400001 to 3450000 to results/GISAID_mutations/human_full-length_spikes_3400001-to-3450000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 3450001 to 3500000 to results/GISAID_mutations/human_full-length_spikes_3450001-to-3500000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 3500001 to 3550000 to results/GISAID_mutations/human_full-length_spikes_3500001-to-3550000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 3550001 to 3600000 to results/GISAID_mutations/human_full-length_spikes_3550001-to-3600000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 3600001 to 3650000 to results/GISAID_mutations/human_full-length_spikes_3600001-to-3650000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 3650001 to 3700000 to results/GISAID_mutations/human_full-length_spikes_3650001-to-3700000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 3700001 to 3750000 to results/GISAID_mutations/human_full-length_spikes_3700001-to-3750000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 3750001 to 3800000 to results/GISAID_mutations/human_full-length_spikes_3750001-to-3800000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 3800001 to 3850000 to results/GISAID_mutations/human_full-length_spikes_3800001-to-3850000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 3850001 to 3900000 to results/GISAID_mutations/human_full-length_spikes_3850001-to-3900000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 3900001 to 3950000 to results/GISAID_mutations/human_full-length_spikes_3900001-to-3950000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 3950001 to 4000000 to results/GISAID_mutations/human_full-length_spikes_3950001-to-4000000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 4000001 to 4050000 to results/GISAID_mutations/human_full-length_spikes_4000001-to-4050000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 4050001 to 4100000 to results/GISAID_mutations/human_full-length_spikes_4050001-to-4100000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 4100001 to 4150000 to results/GISAID_mutations/human_full-length_spikes_4100001-to-4150000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 4150001 to 4200000 to results/GISAID_mutations/human_full-length_spikes_4150001-to-4200000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 4200001 to 4250000 to results/GISAID_mutations/human_full-length_spikes_4200001-to-4250000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 4250001 to 4300000 to results/GISAID_mutations/human_full-length_spikes_4250001-to-4300000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 4300001 to 4350000 to results/GISAID_mutations/human_full-length_spikes_4300001-to-4350000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 4350001 to 4400000 to results/GISAID_mutations/human_full-length_spikes_4350001-to-4400000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 4400001 to 4450000 to results/GISAID_mutations/human_full-length_spikes_4400001-to-4450000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 4450001 to 4500000 to results/GISAID_mutations/human_full-length_spikes_4450001-to-4500000.fasta
    Now aligning these sequences...
    Alignment complete.
    
    Writing spikes 4500001 to 4550000 to results/GISAID_mutations/human_full-length_spikes_4500001-to-4550000.fasta
    Now aligning these sequences...
    Alignment complete.
    


## Parse / filter aligned RBDs

Now put all of the aligned RBDs into a data frame to filter and parse:


```python
rbd_df = (
    pd.DataFrame({'seqrecord': aligned_rbds})
    .assign(description=lambda x: x['seqrecord'].map(lambda rec: rec.description),
            country=lambda x: x['description'].str.split('|').str[-1],
            host=lambda x: x['description'].str.split('|').str[6].str.strip(),
            length=lambda x: x['seqrecord'].map(len),
            n_ambiguous=lambda x: x['seqrecord'].map(lambda rec: rec.seq.count('X') + rec.seq.count('x')),
            n_gaps=lambda x: x['seqrecord'].map(lambda rec: rec.seq.count('-')),
            all_valid_aas=lambda x: x['seqrecord'].map(lambda rec: re.fullmatch(f"[{protein_letters}]+",
                                                                                str(rec.seq)) is not None),
            )
    )

assert all(rbd_df['length'] == len(refseq))
```

Plot number of gaps and ambiguous nucleotides among sequences:


```python
for prop in ['n_ambiguous', 'n_gaps']:
    p = (ggplot(rbd_df) +
         aes(prop) +
         ylab('number of sequences') +
         theme(figure_size=(10, 2.5)) +
         geom_bar()
         )
    _ = p.draw()
```


    
![png](gisaid_rbd_mutations_files/gisaid_rbd_mutations_26_0.png)
    



    
![png](gisaid_rbd_mutations_files/gisaid_rbd_mutations_26_1.png)
    


Based on above plots, we will retain just RBDs with no ambiguous amino acids and no gaps:


```python
rbd_df = rbd_df.query('n_ambiguous == 0').query('n_gaps == 0')
assert rbd_df['all_valid_aas'].all()
print(f"Retained {len(rbd_df)} RBDs.")
```

    Retained 4187147 RBDs.


Now get and plot the number of amino-acid mutations per RBD relative to the reference sequence, plotting on both a linear and log scale.
We then filter all RBDs that have more than some maximum number of mutations, based on the idea that ones that are extremely highly mutated probably are erroneous.
**Note that this maximum number of mutations will change over time, so should be re-assessed periodically by looking at below plots.**


```python
max_muts = 8

refseq_str = str(refseq.seq)
rbd_df = (
    rbd_df
    .assign(seq=lambda x: x['seqrecord'].map(lambda rec: str(rec.seq)),
            n_mutations=lambda x: x['seq'].map(lambda s: sum(x != y for x, y in zip(s, refseq_str))))
    )

p = (ggplot(rbd_df) +
     aes('n_mutations') +
     geom_bar() +
     theme(figure_size=(10, 2.5)) +
     geom_vline(xintercept=max_muts + 0.5, color='red', linetype='dotted')
     )
_ = p.draw()
_ = (p + scale_y_log10()).draw()

rbd_df = rbd_df.query('n_mutations <= @max_muts')
```


    
![png](gisaid_rbd_mutations_files/gisaid_rbd_mutations_30_0.png)
    



    
![png](gisaid_rbd_mutations_files/gisaid_rbd_mutations_30_1.png)
    


Write RBD sequences that pass filtering to a file:


```python
print(f"Overall, there are {len(rbd_df)} aligned RBDs that passed filters.")

rbd_alignment_file = os.path.join(config['gisaid_mutations_dir'], 'RBD_alignment.fasta')
print(f"Writing alignment to {rbd_alignment_file}")
_ = Bio.SeqIO.write(rbd_df['seqrecord'].tolist(), rbd_alignment_file, 'fasta')
```

    Overall, there are 4186099 aligned RBDs that passed filters.
    Writing alignment to results/GISAID_mutations/RBD_alignment.fasta


## Get counts of each mutation
Now we get a data frame that gives the count of each mutation at each site:


```python
records = []
for tup in rbd_df[['seq', 'country']].itertuples():
    for isite, (mut, wt) in enumerate(zip(tup.seq, refseq_str), start=1):
        if mut != wt:
            records.append((isite, isite + config['site_number_offset'], wt, mut, tup.country))
            
muts_df = (pd.DataFrame.from_records(records,
                                     columns=['isite', 'site', 'wildtype', 'mutant', 'country'])
           .groupby(['isite', 'site', 'wildtype', 'mutant'])
           .aggregate(count=pd.NamedAgg('country', 'count'),
                      n_countries=pd.NamedAgg('country', 'nunique'))
           .reset_index()
           .sort_values('count', ascending=False)
           .assign(frequency=lambda x: x['count'] / len(rbd_df))
           )

print('Here are first few lines of mutation counts data frame:')
display(HTML(muts_df.head(n=15).to_html(index=False)))
```

    Here are first few lines of mutation counts data frame:



<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th>isite</th>
      <th>site</th>
      <th>wildtype</th>
      <th>mutant</th>
      <th>count</th>
      <th>n_countries</th>
      <th>frequency</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td>122</td>
      <td>452</td>
      <td>L</td>
      <td>R</td>
      <td>2026360</td>
      <td>627</td>
      <td>0.484069</td>
    </tr>
    <tr>
      <td>148</td>
      <td>478</td>
      <td>T</td>
      <td>K</td>
      <td>1959540</td>
      <td>565</td>
      <td>0.468106</td>
    </tr>
    <tr>
      <td>171</td>
      <td>501</td>
      <td>N</td>
      <td>Y</td>
      <td>1212453</td>
      <td>854</td>
      <td>0.289638</td>
    </tr>
    <tr>
      <td>154</td>
      <td>484</td>
      <td>E</td>
      <td>K</td>
      <td>197486</td>
      <td>554</td>
      <td>0.047177</td>
    </tr>
    <tr>
      <td>87</td>
      <td>417</td>
      <td>K</td>
      <td>T</td>
      <td>87611</td>
      <td>269</td>
      <td>0.020929</td>
    </tr>
    <tr>
      <td>147</td>
      <td>477</td>
      <td>S</td>
      <td>N</td>
      <td>65953</td>
      <td>256</td>
      <td>0.015755</td>
    </tr>
    <tr>
      <td>87</td>
      <td>417</td>
      <td>K</td>
      <td>N</td>
      <td>33120</td>
      <td>241</td>
      <td>0.007912</td>
    </tr>
    <tr>
      <td>109</td>
      <td>439</td>
      <td>N</td>
      <td>K</td>
      <td>27995</td>
      <td>208</td>
      <td>0.006688</td>
    </tr>
    <tr>
      <td>160</td>
      <td>490</td>
      <td>F</td>
      <td>S</td>
      <td>12694</td>
      <td>220</td>
      <td>0.003032</td>
    </tr>
    <tr>
      <td>164</td>
      <td>494</td>
      <td>S</td>
      <td>P</td>
      <td>12147</td>
      <td>143</td>
      <td>0.002902</td>
    </tr>
    <tr>
      <td>16</td>
      <td>346</td>
      <td>R</td>
      <td>K</td>
      <td>11947</td>
      <td>165</td>
      <td>0.002854</td>
    </tr>
    <tr>
      <td>154</td>
      <td>484</td>
      <td>E</td>
      <td>Q</td>
      <td>10317</td>
      <td>132</td>
      <td>0.002465</td>
    </tr>
    <tr>
      <td>122</td>
      <td>452</td>
      <td>L</td>
      <td>Q</td>
      <td>8647</td>
      <td>136</td>
      <td>0.002066</td>
    </tr>
    <tr>
      <td>110</td>
      <td>440</td>
      <td>N</td>
      <td>K</td>
      <td>8299</td>
      <td>105</td>
      <td>0.001983</td>
    </tr>
    <tr>
      <td>190</td>
      <td>520</td>
      <td>A</td>
      <td>S</td>
      <td>6441</td>
      <td>113</td>
      <td>0.001539</td>
    </tr>
  </tbody>
</table>


Plot how many mutations are observed how many times:


```python
p = (ggplot(muts_df) +
     aes('count') +
     geom_histogram(bins=20) +
     scale_x_log10() +
     ylab('number of sequences') +
     xlab('times mutation observed')
     )

_ = p.draw()
```


    
![png](gisaid_rbd_mutations_files/gisaid_rbd_mutations_36_0.png)
    


Write the mutation counts to a file:


```python
print(f"Writing mutation counts to {config['gisaid_mutations_dir']}")
muts_df.to_csv(config['gisaid_mutation_counts'], index=False)
```

    Writing mutation counts to results/GISAID_mutations

