# atomID

[![Continuous Integration](https://github.com/Materials-Data-Science-and-Informatics/atomID/actions/workflows/ci.yml/badge.svg)](https://github.com/Materials-Data-Science-and-Informatics/atomID/actions/workflows/ci.yml)
[![Docs](https://github.com/Materials-Data-Science-and-Informatics/atomID/actions/workflows/docs.yml/badge.svg)](https://github.com/Materials-Data-Science-and-Informatics/atomID/actions/workflows/docs.yml)
[![PyPI - Version](https://img.shields.io/pypi/v/atomid)](https://pypi.org/project/atomID/)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/atomid)](https://img.shields.io/pypi/v/atomid)

Welcome to the `atomID` package! This README will guide you through the initial steps required to start using the `AnnotateCrystal` class for annotating crystal structures and defects. Follow the steps outlined below to get started.

## Overview

"AtomID" is an open-source Python package designed to annotate crystal structure data, enhancing data standardisation and interoperability in the field of computational materials science.


## Key Features

- **Multi-format Data Input**: Supports data import from multiple commonly used formats, such as CIF, POSCAR, and LAMMPS, ensuring compatibility across various modelling tools and workflows.
- **Crystal Structure Identification**: Utilises Common Neighbour Analysis to accurately identify crystal structures within different lattices.
- **Defect Detection and Characterisation**: Applies Wigner-Seitz analysis to detect and categorise various types of lattice defects.
- **Defect Concentration Calculation**: Computes the concentration of different point defects.
- **Knowledge Graph Creation**: Generates a knowledge graph for crystal structures using the Computational Material Sample Ontology (CMSO), which can be exported and stored in Turtle (TTL) format, enabling sharing and complex querying on the data.

## Installation

To begin, you need to install the `atomID` package. This can be done using pip:

```bash
pip install atomid
```

## Usage

This section provides a step-by-step guide on how to utilise the `AnnotateCrystal` class within the `atomID` package.

### 1. Import the Required Class

Start by importing the `AnnotateCrystal` class from the `atomid.annotate` module:

```python
from atomid.annotate import AnnotateCrystal
```

### 2. Create an Instance of `AnnotateCrystal`

Next, create an instance of the `AnnotateCrystal` class:

```python
crystal = AnnotateCrystal()
```

### 3. Read the Crystal Structure File

Read the crystal structure file by using the `read_crystal_structure_file` method. Make sure to replace `crystal_data_file_path` with the actual path to your crystal structure file:

```python
crystal.read_crystal_structure_file(crystal_data_file_path, "vasp")
```

### 4. Annotate the Crystal Structure

You can now annotate the crystal structure with the `annotate_crystal_structure` method:

```python
crystal.annotate_crystal_structure()
```

### 5. Annotate Defects

To annotate defects, provide a reference file path. Replace `ref_file_path` with the actual path to your reference file:

```python
crystal.annotate_defects(ref_file_path, "vasp")
```

### 6. Write to File

Finally, write the annotated data to a file using the `write_to_file` method. Specify the output file name and format:

```python
crystal.write_to_file(output_file_path, "ttl")
```

## Example

Here is a complete example that combines all the steps:

```python
from atomid.annotate import AnnotateCrystal

# Create an instance of AnnotateCrystal
crystal = AnnotateCrystal()

# Read the crystal structure file
crystal.read_crystal_structure_file("path/to/your/interstitial_file.poscar", "vasp")

# Annotate the crystal structure
crystal.annotate_crystal_structure()

# Annotate defects using a reference file
crystal.annotate_defects("path/to/your/reference_file.poscar", "vasp")

# Write the annotated data to a file
crystal.write_to_file("Al_inter.ttl", "ttl")
```

## Documentation

Detailed documentation for the `atomID` package can be found in the [docs](https://github.com/Materials-Data-Science-and-Informatics/atomID/tree/main/docs) folder.

## Contributing

Please refer to the [CONTRIBUTING GUIDE](https://github.com/Materials-Data-Science-and-Informatics/atomID/blob/main/docs/contributing.md) on contributing to the `atomID` package.

## Contact

For any queries or feedback, kindly create an issue on the [GitHub repository](https://github.com/Materials-Data-Science-and-Informatics/atomID).
