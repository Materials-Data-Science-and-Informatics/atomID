from ase.io import read as ase_read
from pyscal3.operations.input import read_inputfile
from pyscal3 import System
from typing import Tuple, Dict, Optional


def read_crystal_structure_file(filename: str, format: Optional[str] = None) -> Tuple:
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

    """

    crystal_structure = ase_read(filename, format=format)

    # Convert to pyscal atoms object
    system = System()
    read_inputfile(system=system, filename=crystal_structure, format="ase")

    return crystal_structure, system


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
