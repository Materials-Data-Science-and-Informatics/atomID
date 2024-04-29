"""Annotate system class."""

import atomrdf as ardf
from ase.io import read as ase_read
from crystal.structure_identification import (
    find_lattice_parameter,
    get_crystal_structure_using_cna,
)
from point_defect_analysis.wigner_seitz_method import analyze_defects


class System:
    """Annotate sytem obect."""

    def __init__(self) -> None:
        self.system = ardf.System()
        self.graph = ardf.KnowledgeGraph()

    def read_crystal_structure_file(self, data_file: str, format: str) -> None:
        """Read the crystal structure file."""
        crystal_data = ase_read(data_file, format=format)
        kg = ardf.KnowledgeGraph()
        crystal_structure = ardf.System.read.file(
            filename=crystal_data, format="ase", graph=kg
        )

        self.graph = kg
        self.system = crystal_structure

    def annotate_crystal_structure(self) -> None:
        """Identify and annotate the crystal structure."""
        crystal_structure = get_crystal_structure_using_cna(self.system)

        if crystal_structure != "others":
            lattice_constants = find_lattice_parameter(self.system, crystal_structure)
            self.system = ardf.System.read.file(
                crystal_structure,
                format="ase",
                graph=self.graph,
                lattice=crystal_structure,
                lattice_constant=lattice_constants,
            )

    def identify_defects(self, reference_data_file: str, ref_format: str) -> dict:
        actual_positions = self.system.atoms.positions
        ref_ase = ase_read(reference_data_file, format=ref_format)
        ref_positions = ref_ase.positions

        defects: dict[str, dict[str, float]] = analyze_defects(
            reference_positions=ref_positions, actual_positions=actual_positions
        )
        return defects

    def annotate_defects(self, reference_data_file: str, ref_format: str) -> None:
        defects = self.identify_defects(reference_data_file, ref_format)

        vacancies = defects.get("Vacancies", {"count": 0, "fraction": 0})

        if vacancies["count"] > 0:
            self.system.add_vacancy(
                concentration=vacancies["fraction"], number=vacancies["count"]
            )
