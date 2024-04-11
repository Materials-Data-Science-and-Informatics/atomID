from ase.io import read as ase_read
from pyscal_rdf import System, KnowledgeGraph
from typing import Tuple, Dict, Optional
import ase
from scipy.signal import find_peaks
from math import sqrt


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

    cna_results = pyscal_system.analysis.common_neighbor_analysis()

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

    val, dist = crystal_system.get_calculator().radial_distribution_function(bins=500)
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
