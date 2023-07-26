# Lateral inhibition in relaxation oscillators provides a basis for computation

<!---
[![Paper DOI : https://doi.org/10.1371/journal.pcbi.1006977](https://badgen.net/badge/PLOS%20Comput%20Bio%20DOI/10.1371%2Fjournal.pcbi.1006977)](https://doi.org/10.1371/journal.pcbi.1006977)
--->

This repository contains data for the manuscript:

> Shamim AP, Menon SN and Sinha S (2022) Lateral inhibition in relaxation oscillators provides a basis for computation.
_arXiv_ :2208.09187
<!---
> https://doi.org/10....
--->

The folders contain program files used to generate the data displayed in the figures of the above manuscript, as well as associated ```.dat``` files.

## Contents of folder **time-series_array**

These codes simulate the dynamics of an array (ring/chain) of N FHN oscillators, using random initial conditions. The code will output the complete time series without removing any transience.

## Contents of folder **steady-state_ensemble**

These codes run an ensemble of simulations on an array (ring/chain) of N FHN oscillators, for different random initial conditions and saves the final u and v values in each case. The data is also binarized using the given thresholds in the paper.

## Contents of folder **unique_binary**

These codes identify the unique binary strings in an array (ring/chain) of N=10 FHN oscillators, starting from the full set of data obtained using folder ```steady_state_ensemble'''.

## Contents of folder **check_stability**

These codes checks the stability of each of the 1024 potential SPOD states in an array of N=10 oscillators. The results obtained here validate the findings in the folder ```unique_binary'''.

## Contents of folder **perturb_SPOD**

These codes simulate the effect of pertubations of specified amplitude and duration to an initial specified SPOD state.

## Contents of folder **transition_matrix**

These codes create a matrix that details the transitions between SPOD states that arise upon perturbations, in a ring/chain of oscillators.

## Contents of folder **compute_stability**

These codes compute the probabilities of obtaining a given SPOD state, starting from random initial conditions (local stability measure), and that obtained upon perturbation (global stability measure). The data used to obtain the local and global stabilities are generated in folders ```unique_binary''' and ```perturb_SPOD'''', respectively.

## Contents of folder **data_files**

1. ```unique_binary_chain_10.dat```: List of unique binary strings for a chain of length N=10.

2. ```unique_binary_ring_perm_10.dat```: List of unique binary strings for a ring of length N=10, under circular permutations.

3. ```stimulation_sites_10.dat```: List of 1024 stimulation protocols used for perturbing the SPOD states in an array of length N=10.

4. ```StrxStr_chain_0.0008_50.dat```: Matrix of size 68x68 containing the frequency of transitions between the SPOD states in a chain.

5. ```StrxStr_ring_0.0008_50.dat```: Matrix of size 14x14 containing the frequency of transitions between the SPOD states in a ring.

6. ```StrxStr_transition_chain_0.0008_50.dat```: Matrix of size 68x68 containing the transition probabilities between the SPOD states in a chain.

7. ```StrxStr_transition_ring_0.0008_50.dat```: Matrix of size 68x68 containing the transition probabilities between the SPOD states in a ring.

