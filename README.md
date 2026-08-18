# On the Geometry of Sixers on the Fermat Cubic Surface – computational files

This repository contains supplementary **Python** and **Singular** code for the article:
**Giuseppe Favacchio, Grzegorz Malara**  
*On the Geometry of Sixers on the Fermat Cubic Surface*

The scripts provide the computational verification of the results concerning sixers on the Fermat cubic surface, their associated groups, and the corresponding plane blow-up models.
A mathematical description of the main computational procedures is given in the appendices of the paper. The source files in this repository contain short comments intended to clarify the main computational steps.

## Contents

### Python scripts

- [`fermat_sixers.py`](fermat_sixers.py)  
  Constructs the 27 lines on the Fermat cubic, determines their skewness graph, enumerates the 72 sixers, computes the action of the automorphism group of the Fermat cubic, and determines the two orbits of sizes 54 and 18.

- [`fermat_normalization.py`](fermat_normalization.py)  
  Normalizes representatives of the two sixer orbits and computes the corresponding matrices over $\mathbb{Q}(\omega)$.

- [`fermat_determinants.py`](fermat_determinants.py)  
  Computes the determinants of the pairwise differences of the normalized matrices and provides the data used to determine the images of the determinant square-class character.

- [`fermat_finite_reductions.py`](fermat_finite_reductions.py)  
  Computes finite reductions of the associated projective matrix groups modulo selected primes and determines the orders of the resulting subgroups of $PGL_2(F_p)$.

- [`fermat_mod13_normalizations.py`](fermat_mod13_normalizations.py)  
  Checks the determinant square-class behavior modulo 13 for all 120 ordered normalization triples of each representative sixer.

### Singular scripts

The Singular computations are divided into three scripts:

- [`fermat_blowup.txt`](fermat_blowup.txt)
  
  Constructs the anticanonical blow-up model of the Fermat cubic, identifies the 27 lines in the standard notation
  $E_i$, $L_{ij}$, $Q_i$, generates the 72 sixers, constructs the second plane blow-down model, and verifies the projective equivalence between the two anticanonical models.

- [`fermat_aut_Z18.txt`](fermat_aut_Z18.txt) 
  Tests all 720 permutations of the six-point configuration $Z_{18}$ and computes its projective automorphisms.

- [`fermat_aut_Z54.txt`](fermat_aut_Z54.txt) 
  Tests all 720 permutations of the six-point configuration $Z_{54}$ and computes its projective automorphisms.

## Requirements

### Python

The Python scripts use only standard Python libraries unless stated otherwise in the individual source files.

### Singular

The Singular computations were carried out and tested in **Singular 4.4.1**.
Where necessary, the scripts contain separate commands intended to accommodate the output format used by earlier versions of Singular.


## Reference

If you use the code contained in this repository, please cite the following paper:
Giuseppe Favacchio, Grzegorz Malara  
*On the Geometry of Sixers on the Fermat Cubic Surface*


## Funding

The work of Giuseppe Favacchio was supported by the funding **PREMIO_SINGOLI_RIC_[2025]** from the Department of Engineering, University of Palermo.
