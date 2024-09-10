# atomID Quickstart Guide

Welcome to the atomID package quickstart guide! This guide will help you get started with using the `AnnotateCrystal` class for annotating crystal structures and defects. Follow the steps below to begin.

## Key Features

- **Multi-format Data Input**: Supports data import from multiple commonly used formats, such as CIF, POSCAR, and LAMMPS, ensuring compatibility across various modelling tools and workflows.
- **Crystal Structure Identification**: Utilises Common Neighbour Analysis to accurately identify crystal structures within different lattices.
- **Defect Detection and Characterisation**: Applies Wigner-Seitz analysis to detect and categorise various types of lattice defects.
- **Defect Concentration Calculation**: Computes the concentration of different point defects.
- **Knowledge Graph Creation**: Generates a knowledge graph for crystal structures using the Computational Material Sample Ontology (CMSO), which can be exported and stored in Turtle (TTL) format, enabling sharing and complex querying on the data.

## Installation

First, ensure you have the `atomID` package installed. You can install it via pip:

```bash
pip install atomid
```

## Usage

Here's a step-by-step guide on how to use the `AnnotateCrystal` class in the atomID package.

### 1. Import the Required Class

Begin by importing the `AnnotateCrystal` class from the `atomid.annotate` module.

```python
from atomid.annotate import AnnotateCrystal
```

### 2. Create an Instance of `AnnotateCrystal`

Create an instance of the `AnnotateCrystal` class.

```python
crystal = AnnotateCrystal()
```

### 3. Read the Crystal Structure File

Use the `read_crystal_structure_file` method to read the crystal structure file. Replace `crystal_data_file_path` with the path to your crystal structure file.

```python
crystal.read_crystal_structure_file(crystal_data_file_path, "vasp")
```

### 4. Annotate the Crystal Structure

Annotate the crystal structure using the `annotate_crystal_structure` method.

```python
crystal.annotate_crystal_structure()
```

### 5. Annotate Defects

Annotate defects by providing a reference file path. Replace `ref_file_path` with the path to your reference file.

```python
crystal.annotate_defects(ref_file_path, "vasp")
```

### 6. Write to File

Finally, write the annotated data to a file using the `write_to_file` method. Specify the output file name and format.

```python
crystal.write_to_file(output_file_path, "ttl")
```

## Example

Below is a complete example putting all the steps together:

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

With these steps, you should be able to start annotating your crystal structures and defects using the atomID package. For more detailed information and advanced usage, please refer to the official documentation.
