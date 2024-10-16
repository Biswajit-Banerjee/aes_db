# Guide: Building Covariance and Hidden Markov Models

This guide will walk you through the steps to build Covariance Models (CMs) and Hidden Markov Models (HMMs) using the tools `cmbuild` and `hmmbuild` respectively.

## Prerequisites
Make sure you have the necessary tools installed:
- **cmbuild**: For building Covariance Models. Two versions of Infernal is being used, infernal-1.0.2 (for cmcompare) and infernal-1.1.5 (for cm_search and cm_scan).
- **hmmbuild**: For building Hidden Markov Models. The [hmmer-3.4](http://eddylab.org/software/hmmer/) is being used. 

Both of these tools typically come with the `Infernal` and `HMMER` software packages, respectively. 

## Input
- **Multiple Sequence Alignment (MSA)**
- **Species name of Anchor Structure**
- **PDB id of anchor structure**
- **Chain Number**
- **Ancestral Expansion Segments(AES) definition on the Anchor Structure**

## Workflow 
- Extract The gapped anchor sequence from MSA
- Extract BasePairs using FR3D API call. [FR3D for 3Q1Q chain B](https://rnacentral.org/api/internal/proxy?url=http://rna.bgsu.edu/rna3dhub/rest/getChainSequenceBasePairs?pdb_id=3q1q&chain=B&only_nested=True)
- Map the Base pairs and AES definitions from Anchor Index to MSA index
- Break The MSA with respect to AES definition
- Convert the AES level MSA to Stockholm with Secondary Structure derived from FR3D

## Building a Covariance Model

To build a Covariance Model (CM) from input alignment, use the following command:

```bash
cmbuild -F <output\cm_file path> <input\stockholm_file path>
```

The cmbuild binary can be found inside src/ of the program base directory.
We only use the -F flag.

## Building a Hidden Markov Model

To build a Hidden Markov Model (HMM) from input alignment, use the following command:

```bash
hmmbuild <output\hmm_file path> <input\stockholm_file path>
```

The hmmbuild binary can be found inside src/ of the program base directory 
