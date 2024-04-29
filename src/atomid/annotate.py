"""Annotate crystal class."""

import atomrdf as ardf
from ase.io import read as ase_read

from atomid.crystal.structure_identification import (
    find_lattice_parameter,
    get_crystal_structure_using_cna,
)
from atomid.point_defect_analysis.wigner_seitz_method import analyze_defects


class AnnotateCrystal:
    """Annotate crystal object."""

    def __init__(self) -> None:
        self.system = ardf.System()
        self.kg = ardf.KnowledgeGraph()

    def read_crystal_structure_file(self, data_file: str, format: str) -> None:
        """Read the crystal structure file.

        Parameters
        ----------
        data_file : str
            The name of the file to read
        format : str
            The format of the file. If None, the format is guessed from the file extension
        """
        crystal_data = ase_read(data_file, format=format)
        kg = ardf.KnowledgeGraph()
        crystal_structure = ardf.System.read.file(
            filename=crystal_data, format="ase", graph=kg
        )

        self.kg = kg
        self.system = crystal_structure

    def annotate_crystal_structure(self) -> None:
        """Identify and annotate the crystal structure.

        This method identifies the crystal structure using Common Neighbour Analysis
        and lattice constant using radial distribution function.
        """
        crystal_structure = get_crystal_structure_using_cna(self.system)

        if crystal_structure != "others":
            lattice_constants = find_lattice_parameter(self.system, crystal_structure)
            self.system = ardf.System.read.file(
                crystal_structure,
                format="ase",
                graph=self.kg,
                lattice=crystal_structure,
                lattice_constant=lattice_constants,
            )

    def identify_defects(self, reference_data_file: str, ref_format: str) -> dict:
        """Identify defects in the crystal structure using the reference data file.

        Parameters
        ----------
        reference_data_file : str
            The name of the file to read
        ref_format : str
            The format of the file. If None, the format is guessed from the file extension

        Returns
        -------
        defects : dict
            A dictionary containing the vacancy and interstitial defects.
        """
        actual_positions = self.system.atoms.positions
        ref_ase = ase_read(reference_data_file, format=ref_format)
        ref_positions = ref_ase.positions

        defects: dict[str, dict[str, float]] = analyze_defects(
            reference_positions_list=ref_positions,
            actual_positions_list=actual_positions,
        )
        return defects

    def annotate_defects(self, reference_data_file: str, ref_format: str) -> None:
        """Annotate defects in the crystal structure using the reference data file.

        Parameters
        ----------
        reference_data_file : str
            The name of the file to read
        ref_format : str
            The format of the file. If None, the format is guessed from the file extension

        """
        defects = self.identify_defects(reference_data_file, ref_format)

        vacancies = defects.get("Vacancies", {"count": 0, "fraction": 0})

        if vacancies["count"] > 0:
            self.system.add_vacancy(
                concentration=vacancies["fraction"], number=vacancies["count"]
            )

    def write_to_file(self, filename: str, format: str = "ttl") -> None:
        """Write the annotated system to a file.

        Parameters
        ----------
        filename : str
            The name of the file to write
        format : str
            The format of the file. If None, the format is guessed from the file extension
        """
        self.kg.write(filename, format=format)
