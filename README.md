# Mock Alexander Polynomial Computations

This repository contains Python code used to check the permanent computations appearing in the examples of the paper

**Connected Sums of Starred Knotoids in S^2 and the Mock Alexander Polynomial**

by Neslihan Gügümcü and Selçuk İlbeyli.

## Description

The mock Alexander polynomial of a starred knotoid can be computed as the permanent of a matrix associated with the diagram. The permanent is similar to the determinant, but without the signs of permutations:

per(A) = sum over all permutations sigma of the product of a_{i,sigma(i)}.

The script in this repository computes the permanents of the matrices used in the examples of the paper and verifies the multiplicativity identities for the connected sums.

In particular, it checks that

nabla(K_1* # K_2*; W) = nabla(K_1*; W) nabla(K_2*; W)

and

nabla(K_3* # K_2*; W) = nabla(K_3*; W) nabla(K_2*; W).

## Files

* `mock_alexander_permanent.py`
  Python script for computing permanents and checking the example computations.

## Requirements

The code requires:

* Python 3
* SymPy

Install SymPy with:

`pip install sympy`

## How to run

Run the script with:

`python mock_alexander_permanent.py`

The script prints the mock Alexander polynomials for the example starred knotoids and verifies the connected-sum multiplicativity identities.

## Output

The script computes:

* nabla(K_1*; W)
* nabla(K_2*; W)
* nabla(K_3*; W)
* nabla(K_1* # K_2*; W)
* nabla(K_3* # K_2*; W)

It also checks whether the permanents of the connected-sum matrices agree with the products of the permanents of the summand matrices.

## Citation

If you use this code, please cite the paper:

N. Gügümcü and S. İlbeyli,
*Connected Sums of Starred Knotoids in S^2 and the Mock Alexander Polynomial*.

## Repository

This repository is intended to make the computations in the examples transparent and reproducible.
