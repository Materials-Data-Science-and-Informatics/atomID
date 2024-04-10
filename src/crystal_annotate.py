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


def find_lattice_parameter(crystal_structure: System, lattice_type: str) -> float:
    """
    Find the lattice parameter of the crystal structure

    Parameters
    ----------
    crystal_structure : ase.Atoms
        The crystal structure

    Returns
    -------
    lattice_parameter : float
        The lattice parameter of the crystal structure
    """

    val, dist = crystal_structure.calculate.radial_distribution_function()
    peaks, _ = find_peaks(val, height=0)
    if lattice_type == "fcc":
        lattice_parameter = dist[peaks[0]] * sqrt(2)
    elif lattice_type == "bcc":
        lattice_parameter = dist[peaks[0]] * sqrt(3)
    elif lattice_type == "hcp":
        lattice_parameter = dist[peaks[0]] * sqrt(2)
    else:
        # return the error stating that the lattice type is not supported
        raise ValueError("Lattice type not supported")

    return float(lattice_parameter)
