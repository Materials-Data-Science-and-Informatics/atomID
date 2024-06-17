"""Annotate crystal class."""

from typing import Optional

import atomrdf as ardf
from ase.io import read as ase_read
from ovito.data import DataCollection
from ovito.io import import_file
from ovito.modifiers import (
    PolyhedralTemplateMatchingModifier,
)

from atomid.crystal.structure_identification import (
    analyse_polyhedral_template_matching_data,
    find_lattice_parameter,
)
from atomid.line_defect_analysis.dislocation_extraction import identify_dislocations
from atomid.plane_defect_analysis.grain_segmentation import identify_grain_orientations
from atomid.point_defect_analysis.wigner_seitz_method import analyze_defects


class AnnotateCrystal:
    """Annotate crystal object."""

    def __init__(
        self,
        data_file: Optional[str] = None,
        format: Optional[str] = None,
        **kwargs: dict[str, str],
    ) -> None:
        """Initialize the AnnotateCrystal object.

        Parameters
        ----------
        data_file : str
            The name of the file to read
        format : str
            The format of the file. If None, the format is guessed from the file extension

        Returns
        -------
        None
        """
        if data_file is not None and format is not None:
            self.read_crystal_structure_file(data_file, format, **kwargs)

    def read_crystal_structure_file(
        self, data_file: str, format: str, **kwargs: dict[str, str]
    ) -> None:
        """Read the crystal structure file.

        Parameters
        ----------
        data_file : str
            The name of the file to read
        format : str
            The format of the file. If None, the format is guessed from the file extension

        Returns
        -------
        None
        """
        self.ase_crystal = ase_read(data_file, format=format, **kwargs)
        kg = ardf.KnowledgeGraph()

        crystal_structure = ardf.System.read.file(
            filename=self.ase_crystal, format="ase", graph=kg
        )

        self.ovito_pipeline = import_file(data_file)

        self.kg = kg
        self.system = crystal_structure

    def validate_parameters_for_crystal_annotation(self) -> None:
        if hasattr(self, "lattice_constant") is False or self.lattice_constant is None:
            self._raise_error(
                "Lattice constants have not been set.",
                "Please run the 'identify_crystal_structure' method first or"
                "set the lattice constants manually.",
            )

        if hasattr(self, "crystal_type") is False or self.crystal_type is None:
            self._raise_error(
                "Crystal structure has not been set.",
                "Please run the 'identify_crystal_structure' method first or"
                "set the crystal structure manually.",
            )

    def _raise_error(self, error_message: str, suggestion: str) -> None:
        full_message = f"\033[91mError: {error_message} {suggestion}\033[0m"
        raise ValueError(full_message)

    def get_polyhedral_template_matching_data(self) -> DataCollection:
        """Get the polyhedral template matching data from the ovito pipeline.

        Parameters
        ----------
        ovito_pipeline : OvitoPipeline
            The ovito pipeline object.

        Returns
        -------
        ovito.data.DataCollection
            The polyhedral template matching data.
        """
        polyhedral_modifier = PolyhedralTemplateMatchingModifier(
            output_interatomic_distance=True, output_orientation=True
        )

        polyhedral_modifier.structures[
            PolyhedralTemplateMatchingModifier.Type.CUBIC_DIAMOND
        ].enabled = True
        polyhedral_modifier.structures[
            PolyhedralTemplateMatchingModifier.Type.HEX_DIAMOND
        ].enabled = True
        polyhedral_modifier.structures[
            PolyhedralTemplateMatchingModifier.Type.SC
        ].enabled = True
        polyhedral_modifier.structures[
            PolyhedralTemplateMatchingModifier.Type.GRAPHENE
        ].enabled = True

        self.ovito_pipeline.modifiers.append(polyhedral_modifier)
        data = self.ovito_pipeline.compute()
        return data

    def identify_crystal_structure(self) -> None:
        """Identify and annotate the crystal structure.

        This method identifies the crystal structure using polyhedral template matching
        and lattice constant using radial distribution function.

        Returns
        -------
        None
        """
        # get crystal structure from polyhedral template matching
        structure_data = self.get_polyhedral_template_matching_data()
        structure_type_atoms = structure_data.particles["Structure Type"][...]  # noqa
        structure_id, crystal_type = analyse_polyhedral_template_matching_data(
            structure_type_atoms
        )
        self.crystal_type = crystal_type

        if crystal_type != "other":
            interatomic_distance = structure_data.particles["Interatomic Distance"][...]  # noqa
            self.lattice_constant = find_lattice_parameter(
                interatomic_distance, structure_type_atoms, int(structure_id)
            )
        else:
            # Warn user that the crystal structure could not be identified
            # and set the lattice constant can not be determined
            print("\033[91mCrystal structure could not be identified.\033[0m")
            print("\033[91mLattice constant can not be determined.\033[0m")
            self.lattice_constant = None

    def set_lattice_constant(self, lattice_constant: float) -> None:
        """Set the lattice constant.

        Parameters
        ----------
        lattice_constant : float
            The lattice constant.

        Returns
        -------
        None
        """
        self.lattice_constant = lattice_constant

    def set_crystal_structure(self, crystal_type: str) -> None:
        """Set the crystal structure.

        Parameters
        ----------
        crystal_type : str
            The crystal structure.

        Returns
        -------
        None
        """
        self.crystal_type = crystal_type

    def annotate_crystal_structure(self) -> None:
        """Annotate the crystal structure.

        Returns
        -------
        None
        """
        self.validate_parameters_for_crystal_annotation()

        self.kg = ardf.KnowledgeGraph()
        self.system = ardf.System.read.file(
            self.ase_crystal,
            format="ase",
            graph=self.kg,
            lattice=self.crystal_type,
            lattice_constant=self.lattice_constant,
        )

    def identify_point_defects(
        self,
        reference_data_file: str,
        ref_format: str,
        method: Optional[str] = None,
        **kwargs: dict[str, str],
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

        ref_ase = ase_read(reference_data_file, format=ref_format, **kwargs)
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

    def identify_line_defects(self) -> tuple[list, list]:
        """Identify line defects in the crystal structure.

        Returns
        -------
        tuple[list, list]
            A tuple containing the burgers vectors and lengths of the dislocations.
        """
        (burgers_vectors, lengths) = identify_dislocations(self.ovito_pipeline)

        if len(burgers_vectors) == 0:
            print("\033[91mNo dislocations found.\033[0m")
            return None, None

        return burgers_vectors, lengths

    def identify_grains(self) -> tuple[list, list]:
        """Identify grains in the crystal structure.

        Returns
        -------
        tuple[list, list]
            A tuple containing the orientations and angles of the grains.
        """
        orientations, angles = identify_grain_orientations(self.ovito_pipeline)

        return orientations, angles

    def annotate_point_defects(
        self, reference_data_file: str, ref_format: str, method: Optional[str] = None
    ) -> None:
        """Annotate defects in the crystal structure using the reference data file.

        Parameters
        ----------
        reference_data_file : str
            The name of the file to read
        ref_format : str
            The format of the file. If None, the format is guessed from the file extension
        method : str
            The method to use for defect identification

        Returns
        -------
        None
        """
        defects = self.identify_point_defects(reference_data_file, ref_format, method)

        vacancies = defects.get("vacancies", {"count": 0, "fraction": 0})
        interstitials = defects.get("interstitials", {"count": 0, "fraction": 0})
        substitutions = defects.get("substitutions", {"count": 0, "fraction": 0})
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

        Returns
        -------
        None
        """
        self.kg.write(filename, format=format)
