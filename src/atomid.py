from ase.io import read as ase_read
from atomrdf import System, KnowledgeGraph
from typing import Tuple, Optional
import ase
from scipy.signal import find_peaks
from math import sqrt
from rdflib import Literal, XSD, URIRef, Namespace, RDF
from point_defect_analysis.wigner_seitz_method import analyze_defects
import logging

CMSO = Namespace("http://purls.helmholtz-metadaten.de/cmso/")
PODO = Namespace("http://purls.helmholtz-metadaten.de/podo/")


def read_crystal_structure_file(
    filename: str, format: Optional[str] = None
) -> Tuple[System, KnowledgeGraph]:
    """
    Read a crystal file and return the pyscal atoms object

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
    system = System(filename=crystal_structure, format="ase", graph=kg)

    return system, kg


def get_crystal_structure_using_cna(pyscal_system: System) -> str:
    """
    Get the crystal structure using adaptive common neighbour analysis.

    Parameters
    ----------
    pyscal_system : pyscal.System
        The pyscal system object.

    Returns
    -------
    str
        The identified crystal structure type.
    """
    cna_results = pyscal_system.analyze.common_neighbor_analysis()
    logging.info("Adaptive common neighbour analysis results: %s", cna_results)

    # Find the most frequent crystal structure from CNA results.
    crystal_type = max(cna_results, key=cna_results.get)

    if crystal_type == "others":
        # Further analyse if the crystal type could be diamond related.
        crystal_type = analyse_diamond_structures(pyscal_system)

    logging.info("Selected crystal structure type: %s", crystal_type)
    return str(crystal_type)


def analyse_diamond_structures(pyscal_system: System) -> str:
    """
    Analyse diamond structures and identify the dominant type if any.

    Parameters
    ----------
    pyscal_system : pyscal.System
        The pyscal system object.

    Returns
    -------
    str
        The dominant diamond structure type or 'others' if none found.
    """
    diamond_results = pyscal_system.analyze.diamond_structure()
    logging.info("Initial diamond structure analysis results: %s", diamond_results)

    # Aggregate the diamond structure counts.
    diamond_results["cubic diamond"] = sum(
        diamond_results.get(key, 0)
        for key in ["cubic diamond", "cubic diamond 1NN", "cubic diamond 2NN"]
    )
    diamond_results["hex diamond"] = sum(
        diamond_results.get(key, 0)
        for key in ["hex diamond", "hex diamond 1NN", "hex diamond 2NN"]
    )

    # Clean up the dictionary.
    for key in [
        "cubic diamond 1NN",
        "cubic diamond 2NN",
        "hex diamond 1NN",
        "hex diamond 2NN",
    ]:
        diamond_results.pop(key, None)

    logging.info("Consolidated diamond structure analysis results: %s", diamond_results)

    # Choose the most frequent diamond structure.
    diamond_type = max(diamond_results, key=diamond_results.get, default="others")
    return str(diamond_type)


def calculate_volume(crystal_structure: ase.Atoms) -> float:
    """
    Calculate the volume of the crystal structure

    Parameters
    ----------
    crystal_structure : ase.Atoms
        The crystal structure

    Returns
    -------
    volume : float
        The volume of the crystal structure
    """

    volume = crystal_structure.get_volume()

    return float(volume)


def find_lattice_parameter(
    crystal_system: System, lattice_type: str
) -> Tuple[float, float, float]:
    """
    Calculate the lattice constants for a given
    crystal structure and lattice type.

    Parameters
    ----------
    crystal_system : atomrdf.System
        The crystal structure as a atomrdf.System object.
    lattice_type : str
        The type of lattice ('fcc', 'bcc', 'hcp').

    Returns
    -------
    Tuple[float, float, float]
        The lattice constants of the crystal structure.

    Raises
    ------
    ValueError
        If the lattice type is not supported.
    """

    val, dist = crystal_system.calculate.radial_distribution_function(bins=500)
    peaks, _ = find_peaks(val, height=0)

    lattice_calculation = {
        "fcc": lambda d: (d * sqrt(2),) * 3,
        "bcc": lambda d: (d * 2 / sqrt(3),) * 3,
        "hcp": lambda d: (d, d, d * 2 * sqrt(6) / 3),
        "cubic diamond": lambda d: (d * 4 / sqrt(3),) * 3,
    }

    if lattice_type not in lattice_calculation:
        raise ValueError(f"Lattice type '{lattice_type}' not supported")

    lattice_constants = lattice_calculation[lattice_type](dist[peaks[0]])
    # Round to 3 decimal places
    lattice_constants = tuple(round(lc, 3) for lc in lattice_constants)
    return lattice_constants


def get_lattice_angle(lattice_type: str) -> Tuple[float, float, float]:
    """
    Get the lattice angles of the crystal structure

    Parameters
    ----------
    crystal_structure : ase.Atoms
        The crystal structure

    Returns
    -------
    Tuple[float, float, float]
        The lattice angles of the crystal structure
    """
    lattice_angle = {
        "fcc": (90, 90, 90),
        "bcc": (90, 90, 90),
        "hcp": (90, 90, 120),
        "cubic diamond": (90, 90, 90),
    }

    return lattice_angle[lattice_type]


def add_lattice_parameter(
    graph: KnowledgeGraph, system_name: str, values: Tuple[float, float, float]
) -> None:
    properties = ["hasLength_x", "hasLength_y", "hasLength_z"]
    for prop, value in zip(properties, values):
        graph.add(
            (
                URIRef(f"{system_name}_LatticeParameter"),
                CMSO[prop],
                Literal(value, datatype=XSD.float),
            )
        )


def add_lattice_angle(
    graph: KnowledgeGraph, system_name: str, angles: Tuple[float, float, float]
) -> None:
    properties = ["hasAngle_alpha", "hasAngle_beta", "hasAngle_gamma"]
    for prop, angle in zip(properties, angles):
        graph.add(
            (
                URIRef(f"{system_name}_LatticeAngle"),
                CMSO[prop],
                Literal(angle, datatype=XSD.float),
            )
        )


def get_bravis_lattice_type(crystal_structure_type: str) -> str:
    """
    Get the Bravais lattice type of the crystal structure

    Parameters
    ----------
    crystal_structure : ase.Atoms
        The crystal structure

    Returns
    -------
    lattice_type : str
        The Bravais lattice type of the crystal structure
    """
    bravais_lattice = {
        "fcc": "https://www.wikidata.org/wiki/Q3006714",
        "bcc": "https://www.wikidata.org/wiki/Q851536",
        "hcp": "https://www.wikidata.org/wiki/Q663314",
        "cubic diamond": "https://www.wikidata.org/wiki/Q3006714",
    }
    return bravais_lattice[crystal_structure_type]

def convert_plural_to_singular(form: str) -> str:
    """
    Convert if the form is plural and convert it to singular

    Parameters
    ----------
    form : str
        The form to check

    Returns
    -------
    str
        The singular form
    """
    mapping_plural_singular = {
        "Vacancies": "Vacancy",
        "Interstitials": "Interstitial",
        "Substitutionals": "Substitutional",
    }

    if form not in mapping_plural_singular:
        logging.warning(f"Form '{form}' not found in mapping.")
    return mapping_plural_singular.get(form, form)

def get_space_group(crystal_number_type: str) -> Tuple[int, str]:
    """
    Get the space group number of the crystal structure

    Parameters
    ----------
    crystal_structure_type : str
        The crystal structure type

    Returns
    -------
    lattice_type : str
        The space group number of the crystal structure
    """
    space_group = {
        "fcc": (225, "Fm-3m"),
        "bcc": (229, "Im-3m"),
        "hcp": (194, "P63/mmc"),
        "cubic diamond": (227, "Fd-3m"),
    }
    return space_group[crystal_number_type]



def add_defects_to_graph(kg: KnowledgeGraph, system_name: str, defects: dict) -> None:
    """
    Add defect data to the knowledge graph.

    Parameters:
    kg (KnowledgeGraph): The knowledge graph instance.
    system_name (str): The name of the system being modelled.
    defects (dict): A dictionary with defect types as keys and dictionaries with 'count' and 'fraction' as values.
    """
    for defect_type, defect_info in defects.items():
        if defect_info["count"] > 0:
            add_defect_relations(kg, system_name, defect_type, defect_info)

def add_defect_relations(kg: KnowledgeGraph, system_name: str, defect_type: str, defect_info: dict):
    """Helper function to add defect relations to the knowledge graph."""
    singular_defect_type = convert_plural_to_singular(defect_type)
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
        (URIRef(f"{system_name}_{singular_defect_type}"), RDF.type, PODO[singular_defect_type])
    )




def annotate_defects(
    input_turtle_file: str, reference_data_file: str, output_file: str
) -> None:
    """Annotates defects in the material sample described in the input turtle file using the reference data file."""
    kg = KnowledgeGraph(graph_file=input_turtle_file)
    sample = kg.samples[0]
    system_name = str(sample)
    system = kg.get_system_from_sample(sample)

    actual_positions = system.atoms.positions
    ref_system, _ = read_crystal_structure_file(reference_data_file, format="cif")
    ref_positions = ref_system.atoms.positions

    defects = analyze_defects(
        reference_positions=ref_positions, actual_positions=actual_positions
    )
    add_defects_to_graph(kg, system_name, defects)

    kg.write(output_file, format="ttl")


def annotate_crystal_structure(data_file: str, format: str, output_file: str) -> None:
    """Annotate the crystal structure using pyscal and save the results to a knowledge graph

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
    None"""
    system, kg = read_crystal_structure_file(data_file, format="cif")
    system.to_graph()

    crystal_type = get_crystal_structure_using_cna(system)

    kg.graph.add(
        (
            URIRef(f"{system._name}_CrystalStructure"),
            CMSO["hasAltName"],
            Literal(crystal_type, datatype=XSD.string),
        )
    )

    if crystal_type != "others":
        space_group_number, space_group_symbol = get_space_group(crystal_type)
        kg.graph.add(
            (
                URIRef(f"{system._name}_SpaceGroup"),
                CMSO["hasSpaceGroupNumber"],
                Literal(space_group_number),
            )
        )
        kg.graph.add(
            (
                URIRef(f"{system._name}_SpaceGroup"),
                CMSO["hasSpaceGroupSymbol"],
                Literal(space_group_symbol, datatype=XSD.string),
            )
        )
        bravice_lattice = get_bravis_lattice_type(crystal_type)
        kg.graph.add(
            (
                URIRef(f"{system._name}_UnitCell"),
                CMSO["hasBravaisLattice"],
                URIRef(bravice_lattice),
            )
        )
        lattice_constants = find_lattice_parameter(system, crystal_type)
        lattice_angles = get_lattice_angle(crystal_type)

        add_lattice_parameter(kg.graph, system._name, lattice_constants)
        add_lattice_angle(kg.graph, system._name, lattice_angles)

    kg.write(output_file, format="ttl")
