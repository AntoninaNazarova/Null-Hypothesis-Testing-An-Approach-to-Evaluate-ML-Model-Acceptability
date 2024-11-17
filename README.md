ESOL Dataset & Sphere-Exclusion and Kennard-Stone-Based Dataset Partitioning for Training, Test, and Evaluation Sets


This repository contains datasets and source codes submitted as part of a JCIM submission. The provided files are associated with the dataset partitioning approaches and the feedforward neural network design used in the study.

Contents
1. Sphere-Exclusion-Based Partitioning
SE_based_ESOL_partitioning_training_test_evaluation.cpp:
Source code for sphere-exclusion-based dataset partitioning into training, test, and evaluation sets. This code is a fragment of the designed feedforward neural network.
2. Kennard-Stone-Based Partitioning
KS_based_ESOL_partitioning_1.cpp and KS_based_ESOL_partitioning_2.cpp:
Source codes for Kennard-Stone-based dataset partitioning into training, test, and evaluation sets. These codes are fragments of the designed feedforward neural network.

Kennard-Stone_based_ESOL_ranking.dat:
The Kennard-Stone ranked ESOL dataset used for partitioning.

3. ESOL Dataset and Supporting Files
ESOL_smiles_eps.txt:
The ESOL dataset obtained from the MoleculeNet portal, used in this study.

s_eps.dat:
A dataset containing solubility values of ESOL compounds.

s_smiles.dat:
A dataset of VLA-SMILES descriptors (with d=4) for ESOL compounds. The VLA-SMILES descriptors were generated using the methodology described in the following publication:
Mach. Learn. Knowl. Extr. 2022, 4(3), 715-737.
https://doi.org/10.3390/make4030034
