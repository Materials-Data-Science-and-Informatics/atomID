"""Annotate crystal class."""

from typing import Optional

import atomrdf as ardf
from ase.io import read as ase_read
from ovito.data import DataCollection
from ovito.io import import_file
from ovito.modifiers import PolyhedralTemplateMatchingModifier

from atomid.crystal.structure_identification import (
    analyse_polyhedral_template_matching_data,
    find_lattice_parameter_2,
    get_crystal_structure_using_cna,
)
from atomid.point_defect_analysis.wigner_seitz_method import analyze_defects


class AnnotateCrystal:
    """Annotate crystal object."""

    def __init__(self) -> None:
        self.system = ardf.System()
        self.kg = ardf.KnowledgeGraph()
        self.ase_crystal = None

    def read_crystal_structure_file(self, data_file: str, format: str) -> None:
        """Read the crystal structure file.

        Parameters
        ----------
        data_file : str
            The name of the file to read
        format : str
            The format of the file. If None, the format is guessed from the file extension
        """
        self.ase_crystal = ase_read(data_file, format=format)
        kg = ardf.KnowledgeGraph()

        crystal_structure = ardf.System.read.file(
            filename=self.ase_crystal, format="ase", graph=kg
        )

        self.ovito_pipeline = import_file(data_file)

        self.kg = kg
        self.system = crystal_structure

    def get_polyhedral_template_matching_data(self) -> DataCollection:
        """Get the polyhedral template matching data from the ovito pipeline.

        Parameters
        ----------
        ovito_pipeline : OvitoPipeline
            The ovito pipeline object.

        Returns
        -------
        dict
            The polyhedral template matching data.
        """
        self.ovito_pipeline.modifiers.append(
            PolyhedralTemplateMatchingModifier(output_interatomic_distance=True)
        )
        data = self.ovito_pipeline.compute()
        return data

    def annotate_crystal_structure(self) -> None:
        """Identify and annotate the crystal structure.

        This method identifies the crystal structure using Common Neighbour Analysis
        and lattice constant using radial distribution function.
        """
        crystal_type = get_crystal_structure_using_cna(self.system)
        # get crystal structure from polyhedral template matching
        structure_data = self.get_polyhedral_template_matching_data()
        structure_type_atoms = structure_data.particles["Structure Type"][...]  # noqa
        structure_id, crystal_type = analyse_polyhedral_template_matching_data(
            structure_type_atoms
        )

        if crystal_type != "others":
            interatomic_distance = structure_data.particles["Interatomic Distance"][...]  # noqa
            lattice_constants = find_lattice_parameter_2(
                interatomic_distance, structure_type_atoms, int(structure_id)
            )
            # lattice_constants = find_lattice_parameter(self.system, crystal_type)
            self.system = None
            self.kg = ardf.KnowledgeGraph()
            self.system = ardf.System.read.file(
                self.ase_crystal,
                format="ase",
                graph=self.kg,
                lattice=crystal_type,
                lattice_constant=lattice_constants,
            )

    def identify_defects(
        self, reference_data_file: str, ref_format: str, method: Optional[str] = None
    ) -> dict:
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
        species_reference = ref_ase.get_chemical_symbols()
        species_actual = self.system.atoms["species"]
        defects: dict[str, dict[str, float]] = analyze_defects(
            reference_positions=ref_positions,
            actual_positions=actual_positions,
            method=method,
            species_ref=species_reference,
            species_actual=species_actual,
        )
        return defects

    def annotate_defects(
        self, reference_data_file: str, ref_format: str, method: Optional[str] = None
    ) -> None:
        """Annotate defects in the crystal structure using the reference data file.

        Parameters
        ----------
        reference_data_file : str
            The name of the file to read
        ref_format : str
            The format of the file. If None, the format is guessed from the file extension

        """
        defects = self.identify_defects(reference_data_file, ref_format, method)

        vacancies = defects.get("Vacancies", {"count": 0, "fraction": 0})
        interstitials = defects.get("Interstitials", {"count": 0, "fraction": 0})
        substitutions = defects.get("Substitutions", {"count": 0, "fraction": 0})
        if vacancies["count"] > 0:
            self.system.add_vacancy(
                concentration=vacancies["fraction"], number=vacancies["count"]
            )

        if interstitials["count"] > 0:
            self.system.add_triples_for_interstitial_impurities(
                conc_of_impurities=interstitials["fraction"],
                no_of_impurities=interstitials["count"],
            )

        if substitutions["count"] > 0:
            self.system.add_triples_for_substitutional_impurities(
                conc_of_impurities=substitutions["fraction"],
                no_of_impurities=substitutions["count"],
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
