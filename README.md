# 🧬 Guide: Ribosome Building Blocks analysis

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Infernal](https://img.shields.io/badge/Infernal-1.1.5-blue)](http://eddylab.org/infernal/)
[![HMMER](https://img.shields.io/badge/HMMER-3.4-green)](http://eddylab.org/software/hmmer/)

This comprehensive guide walks you through the process of building Covariance Models (CMs) and Hidden Markov Models (HMMs) using the powerful tools `cmbuild` and `hmmbuild`. These models are crucial for various bioinformatics applications, particularly in RNA sequence analysis.

## 📋 Table of Contents

- [Prerequisites](#prerequisites)
- [Input Requirements](#input-requirements)
- [Workflow Overview](#workflow-overview)
- [Building Models](#building-models)
  - [Covariance Models (CMs)](#covariance-models-cms)
  - [Hidden Markov Models (HMMs)](#hidden-markov-models-hmms)
- [Additional Tools](#additional-tools)
- [Contributing](#contributing)
- [License](#license)

## 🛠 Prerequisites

Ensure you have the following tools installed:

- **cmbuild**: For building Covariance Models
  - Infernal-1.0.2 (for `cmcompare`)
  - Infernal-1.1.5 (for `cm_search` and `cm_scan`)
- **hmmbuild**: For building Hidden Markov Models
  - [HMMER-3.4](http://eddylab.org/software/hmmer/)

These tools are typically included in the `Infernal` and `HMMER` software packages, respectively.

## 📥 Input Requirements

To build your models, you'll need the following inputs:

1. **Multiple Sequence Alignment (MSA)**
2. **Species name of Anchor Structure**
3. **PDB ID of anchor structure**
4. **Chain Number**
5. **Ancestral Expansion Segments (AES) definition on the Anchor Structure**

## 🔄 Workflow Overview

1. Extract the gapped anchor sequence from MSA
2. Extract BasePairs using FR3D API call
   - Example: [FR3D for 3Q1Q chain B](https://rnacentral.org/api/internal/proxy?url=http://rna.bgsu.edu/rna3dhub/rest/getChainSequenceBasePairs?pdb_id=3q1q&chain=B&only_nested=True)
3. Map the Base pairs and AES definitions from Anchor Index to MSA index
4. Break the MSA with respect to AES definition
5. Convert the AES level MSA to Stockholm format with Secondary Structure derived from FR3D

## 🏗 Building Models

### Covariance Models (CMs)

To build a Covariance Model from your input alignment:

```bash
cmbuild -F <output/cm_file_path> <input/stockholm_file_path>
```

> **Note**: The `cmbuild` binary can be found in the `src/` directory of the program base directory. We only use the `-F` flag.

### Hidden Markov Models (HMMs)

To build a Hidden Markov Model from your input alignment:

```bash
hmmbuild <output/hmm_file_path> <input/stockholm_file_path>
```

> **Note**: The `hmmbuild` binary can be found in the `src/` directory of the program base directory.

## 🔬 Additional Tools

### Running CM Compare

To compare two Covariance Models:

```bash
cmcompare <input/cm1_file_path> <input/cm2_file_path>
```

> **Note**: Comparing cm1 to cm2 and cm2 to cm1 are different. They are not associative (A . B != B . A).

## 🤝 Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

For more information on Infernal and HMMER, visit their official websites:
- [Infernal](http://eddylab.org/infernal/)
- [HMMER](http://hmmer.org/)
- [CMcompare](https://github.com/choener/CMCompare)