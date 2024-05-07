"""Main module for the atomid package."""

import logging
from typing import Optional, Tuple

import ase
from ase.io import read as ase_read
from atomrdf import KnowledgeGraph, System
from rdflib import Namespace

from atomid.crystal.structure_identification import (
    find_lattice_parameter,
    get_crystal_structure_using_cna,
)
from atomid.point_defect_analysis.wigner_seitz_method import analyze_defects

CMSO = Namespace("http://purls.helmholtz-metadaten.de/cmso/")
PODO = Namespace("http://purls.helmholtz-metadaten.de/podo/")


def read_crystal_structure_file(
    filename: str, format: Optional[str] = None
) -> Tuple[ase.Atoms, System, KnowledgeGraph]:
    """
    Read a crystal file and return the pyscal atoms object.

    Parameters
    ----------
    filename : str
        The name of the file to read
    format : str
        The format of the file. If None, the format is guessed from the file extension

    Returns
    -------
    crystal_structure : ase.Atoms
        The crystal structure
    system : pyscal.System
        The pyscal system object
    kg : pyscal.KnowledgeGraph
        The pyscal knowledge graph object

    """
    crystal_structure = ase_read(filename, format=format)
    kg = KnowledgeGraph()
    # Convert to pyscal atoms object
    system = System.read.file(filename=crystal_structure, format="ase", graph=kg)

    return crystal_structure, system, kg


def identify_defects_in_crystal_structure(
    pyscal_system: System, reference_data_file: str, ref_format: Optional[str] = None
) -> dict[str, dict[str, float]]:
    """Annotates defects in the crystal structure using the reference data file."""
    actual_positions = pyscal_system.atoms.positions
    ref_ase = ase_read(reference_data_file, format=ref_format)
    ref_positions = ref_ase.positions

    defects: dict[str, dict[str, float]] = analyze_defects(
        reference_positions=ref_positions, actual_positions=actual_positions
    )

    return defects


def annotate_crystal_structure(
    data_file: str,
    format: str,
    output_file: str,
    reference_data_file: Optional[str] = None,
    ref_format: Optional[str] = None,
) -> None:
    """Annotate the crystal structure using pyscal and save the results to a knowledge graph.

    Parameters
    ----------
    data_file : str
        The path to the crystal structure file
    format : str
        The format of the crystal structure file
    output_file : str
        The path to save the knowledge graph

    Returns
    -------
    None
    """
    crystal_structure, system, kg = read_crystal_structure_file(
        data_file, format=format
    )
    system.to_graph()

    crystal_type = get_crystal_structure_using_cna(system)

    if crystal_type != "others":
        kg = KnowledgeGraph()
        lattice_constants = find_lattice_parameter(system, crystal_type)
        System.read.file(
            crystal_structure,
            format="ase",
            graph=kg,
            lattice=crystal_type,
            lattice_constant=lattice_constants,
        )

    if reference_data_file:
        logging.info("Annotating defects in the crystal structure")
        defects = identify_defects_in_crystal_structure(
            system, reference_data_file, ref_format
        )
        vacancies = defects.get("Vacancies", {"count": 0, "fraction": 0})

        if vacancies["count"] > 0:
            system.add_vacancy(
                concentration=vacancies["fraction"], number=vacancies["count"]
            )
    kg.write(output_file, format="ttl")
