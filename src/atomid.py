"""Main module for the atomid package."""

import logging
from typing import Optional, Tuple

import ase
from ase.io import read as ase_read
from atomrdf import KnowledgeGraph, System
from crystal_analysis.structure_identification import (
    find_lattice_parameter,
    get_crystal_structure_using_cna,
)
from general_utils import convert_plural_to_singular
from point_defect_analysis.wigner_seitz_method import analyze_defects
from rdflib import RDF, XSD, Literal, Namespace, URIRef

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


def add_defects_to_graph(kg: KnowledgeGraph, system_name: str, defects: dict) -> None:
    """
    Add defect data to the knowledge graph.

    Parameters
    kg (KnowledgeGraph): The knowledge graph instance.
    system_name (str): The name of the system being modelled.
    defects (dict): A dictionary with defect types as keys and dictionaries with 'count' and 'fraction' as values.
    """
    for defect_type, defect_info in defects.items():
        if defect_info["count"] > 0:
            add_defect_relations(kg, system_name, defect_type, defect_info)


def add_defect_relations(
    kg: KnowledgeGraph, system_name: str, defect_type: str, defect_info: dict
) -> None:
    """Add defects to the knowledge graph."""
    singular_defect_type: str = convert_plural_to_singular(defect_type)
    kg.graph.add(
        (
            URIRef(f"{system_name}_Material"),
            CMSO["hasDefect"],
            URIRef(f"{system_name}_{defect_type}"),
        )
    )
    kg.graph.add(
        (
            URIRef(f"{system_name}_SimulationCell"),
            PODO[f"hasNumberOf{defect_type}"],
            Literal(defect_info["count"], datatype=XSD.integer),
        )
    )
    kg.graph.add(
        (
            URIRef(f"{system_name}_SimulationCell"),
            PODO[f"has{singular_defect_type}Concentration"],
            Literal(defect_info["fraction"], datatype=XSD.float),
        )
    )
    kg.graph.add(
        (
            URIRef(f"{system_name}_{singular_defect_type}"),
            RDF.type,
            PODO[singular_defect_type],
        )
    )


def annotate_defects(
    input_turtle_file: str, reference_data_file: str, output_file: str
) -> None:
    """Annotates defects in the material sample described in the input turtle file using the reference data file."""
    kg_input = KnowledgeGraph(graph_file=input_turtle_file)
    sample = kg_input.samples[0]
    input_system = kg_input.get_system_from_sample(sample)

    actual_positions = input_system.atoms.positions
    _, ref_system, _ = read_crystal_structure_file(reference_data_file, format="cif")
    ref_positions = ref_system.atoms.positions

    defects = analyze_defects(
        reference_positions=ref_positions, actual_positions=actual_positions
    )
    vacancies = defects.get("Vacancies", {"count": 0, "fraction": 0})

    input_system._name = str(sample)
    print(input_system.simulation_cell)
    if vacancies["count"] > 0:
        input_system.add_vacancy(
            concentration=vacancies["fraction"], number=vacancies["count"]
        )
    print(f"Vacancies: {vacancies}")
    print(input_system.graph)
    # print(kg_input.graph.serialize(format="turtle"))
    kg_input.write(output_file, format="ttl")


def annotate_defects_in_crystal_structure(
    pyscal_system: System, reference_data_file: str, ref_format: str
) -> KnowledgeGraph:
    """Annotates defects in the crystal structure using the reference data file."""
    actual_positions = pyscal_system.atoms.positions
    _, ref_system, _ = read_crystal_structure_file(
        reference_data_file, format=ref_format
    )
    ref_positions = ref_system.atoms.positions

    defects = analyze_defects(
        reference_positions=ref_positions, actual_positions=actual_positions
    )
    vacancies = defects.get("Vacancies", {"count": 0, "fraction": 0})

    if vacancies["count"] > 0:
        print("graph", pyscal_system.graph)
        pyscal_system.add_vacancy(
            concentration=vacancies["fraction"], number=vacancies["count"]
        )
        pyscal_system.graph.write("def_out.ttl", format="ttl")

    return pyscal_system.graph
    print(f"Vacancies: {vacancies}")


def annotate_crystal_structure(
    data_file: str,
    format: str,
    output_file: str,
    annotate_defects: Optional[bool] = None,
    reference_data_file: Optional[str] = None,
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
    crystal_structure, system, kg = read_crystal_structure_file(data_file, format="cif")
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

    if annotate_defects:
        logging.info("Annotating defects in the crystal structure")
        kg = annotate_defects_in_crystal_structure(system, reference_data_file, format)

    kg.write(output_file, format="ttl")
