from ase.io import read as ase_read
from pyscal_rdf import System, KnowledgeGraph
from typing import Tuple, Dict, Optional
import ase
from scipy.signal import find_peaks
from math import sqrt
from rdflib import Literal, XSD, URIRef, Namespace

CMSO = Namespace("http://purls.helmholtz-metadaten.de/cmso/")


def read_crystal_structure_file(
    filename: str, format: Optional[str] = None
) -> Tuple[ase.Atoms, System, KnowledgeGraph]:
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

    return crystal_structure, system, kg


def get_crystal_structure_using_cna(pyscal_system: System) -> Dict[str, int]:
    """
    Get the crystal structure using adaptive common neighbor analysis

    Parameters
    ----------
    pyscal_system : pyscal.System
        The pyscal system object

    Returns
    -------
    cna_results : dict
        The results of the common neighbor analysis
    """

    cna_results = pyscal_system.analyze.common_neighbor_analysis()

    return dict(cna_results)


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
    crystal_system : pyscal_rdf.System
        The crystal structure as a pyscal_rdf.System object.
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
    }

    if lattice_type not in lattice_calculation:
        raise ValueError(f"Lattice type '{lattice_type}' not supported")

    lattice_constants = lattice_calculation[lattice_type](dist[peaks[0]])

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
    }

    return lattice_angle[lattice_type]


if __name__ == "__main__":
    crystal_structure, system, kg = read_crystal_structure_file(
        "/home/n.bhat/development/data/fcc/no_defects/Al/Al.cif", format="cif"
    )
    system.to_graph()
    cna_results = get_crystal_structure_using_cna(system)
    # Pick key with maximum value from cn_results
    crystal_type = max(cna_results, key=cna_results.__getitem__)

    volume = calculate_volume(crystal_structure)
    if crystal_type != "other":
        lattice_constants = find_lattice_parameter(system, crystal_type)
        lattice_angles = get_lattice_angle(crystal_type)

        kg.graph.add(
            (
                URIRef(f"{system._name}_LatticeParameter"),
                CMSO.hasLength_x,
                Literal(lattice_constants[0], datatype=XSD.float),
            )
        )
        kg.graph.add(
            (
                URIRef(f"{system._name}_LatticeParameter"),
                CMSO.hasLength_y,
                Literal(lattice_constants[1], datatype=XSD.float),
            )
        )
        kg.graph.add(
            (
                URIRef(f"{system._name}_LatticeParameter"),
                CMSO.hasLength_z,
                Literal(lattice_constants[2], datatype=XSD.float),
            )
        )

        kg.graph.add(
            (
                URIRef(f"{system._name}_LatticeAngle"),
                CMSO.hasAngle_alpha,
                Literal(lattice_angles[0], datatype=XSD.float),
            )
        )
        kg.graph.add(
            (
                URIRef(f"{system._name}_LatticeAngle"),
                CMSO.hasAngle_beta,
                Literal(lattice_angles[1], datatype=XSD.float),
            )
        )
        kg.graph.add(
            (
                URIRef(f"{system._name}_LatticeAngle"),
                CMSO.hasAngle_gamma,
                Literal(lattice_angles[2], datatype=XSD.float),
            )
        )

    kg.write("crystal_structure_test.ttl", format="ttl")
