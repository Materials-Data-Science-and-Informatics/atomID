"""Annotate crystal class."""

from typing import Any, Optional

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
        kg: Optional[ardf.KnowledgeGraph] = None,
        **kwargs: Any,
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
        if kg is None:
            self.kg = ardf.KnowledgeGraph()
        else:
            self.kg = kg

        if data_file is not None:
            if format is not None:
                self.read_crystal_structure_file(data_file, format, **kwargs)
            else:
                self.read_crystal_structure_file(data_file, **kwargs)

    def read_crystal_structure_file(
        self, data_file: str, format: Optional[str] = None, **kwargs: dict[str, str]
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

        self.ovito_pipeline = import_file(data_file)

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

    def identify_crystal_structure(self, log: bool = True) -> None:
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
            if log:
                print(f"\033[92mCrystal structure: {crystal_type}\033[0m")
                print(f"\033[92mLattice constant: {self.lattice_constant}\033[0m")
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
        log: bool = True,
        **kwargs: dict[str, str],
    ) -> None:
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

        if hasattr(self, "defects") is False:
            self.defects = dict()

        self.defects.update(defects)
        # Print identified defects
        if log:
            print("\033[92mIdentified defects:\033[0m")
        for defect, defect_info in defects.items():
            if defect_info["count"] != 0:
                if log:
                    print(
                        f"\033[92m{defect}:\033[0m Count: {defect_info['count']} "
                        f"Concentration: {defect_info['concentration']:.2f}"
                    )

    def add_vacancy_information(self, concentration: float, number: int) -> None:
        """Add vacancy information to the system.

        Parameters
        ----------
        concentration : float
            The concentration of vacancies.
        number : int
            The number of vacancies.

        Returns
        -------
        None
        """
        if hasattr(self, "defects") is False:
            self.defects = {}

        self.defects.update(
            {"vacancies": {"concentration": concentration, "count": number}}
        )

    def add_interstitial_information(self, concentration: float, number: int) -> None:
        """Add interstitial information to the system.

        Parameters
        ----------
        concentration : float
            The concentration of interstitials.
        number : int
            The number of interstitials.

        Returns
        -------
        None
        """
        if hasattr(self, "defects") is False:
            self.defects = {}

        self.defects.update(
            {"interstitials": {"concentration": concentration, "count": number}}
        )

    def add_substitution_information(self, concentration: float, number: int) -> None:
        """Add substitution information to the system.

        Parameters
        ----------
        concentration : float
            The concentration of substitutions.
        number : int
            The number of substitutions.

        Returns
        -------
        None
        """
        if hasattr(self, "defects") is False:
            self.defects = {}

        self.defects.update(
            {"substitutions": {"concentration": concentration, "count": number}}
        )

    def annotate_point_defects(self) -> None:
        """Annotate defects in the crystal structure using the reference data file.

        Returns
        -------
        None
        """
        if hasattr(self, "defects") is False:
            self.defects = {}

        defects = self.defects

        vacancies = defects.get("vacancies", {"count": 0, "concentration": 0})
        interstitials = defects.get("interstitials", {"count": 0, "concentration": 0})
        substitutions = defects.get("substitutions", {"count": 0, "concentration": 0})
        if vacancies["count"] > 0:
            self.system.add_vacancy(
                concentration=vacancies["concentration"], number=vacancies["count"]
            )

        if interstitials["count"] > 0:
            self.system.add_triples_for_interstitial_impurities(
                conc_of_impurities=interstitials["concentration"],
                no_of_impurities=interstitials["count"],
            )

        if substitutions["count"] > 0:
            self.system.add_triples_for_substitutional_impurities(
                conc_of_impurities=substitutions["concentration"],
                no_of_impurities=substitutions["count"],
            )

    def identify_line_defects(self) -> None:
        """Identify line defects in the crystal structure.

        Returns
        -------
        tuple[list, list]
            A tuple containing the burgers vectors and lengths of the dislocations.
        """
        (burgers_vectors, lengths) = identify_dislocations(self.ovito_pipeline)

        if len(burgers_vectors) == 0:
            print("\033[91mNo dislocations found.\033[0m")
            return None
        # Lengths to 4 decimal places
        print(
            "\033[92mDislocations found:\033[0m"
            f"\nBurgers vectors: {burgers_vectors}"
            f"\nLengths: {[round(length, 4) for length in lengths]}"
        )
        self.burgers_vectors = burgers_vectors
        # round the lengths to 4 decimal places
        self.dislocation_lengths = [round(length, 4) for length in lengths]

    def add_dislocation_information(self, burgers_vectors: list, lengths: list) -> None:
        """Add dislocation information to the system.

        Parameters
        ----------
        burgers_vectors : list
            A list of burgers vectors.
        lengths : list
            A list of dislocation lengths.

        Returns
        -------
        None
        """
        # add to burgers vectors if not already present
        if hasattr(self, "burgers_vectors") is False:
            self.burgers_vectors = []
        self.burgers_vectors.extend(burgers_vectors)
        # add to dislocation lengths if not already present
        if hasattr(self, "dislocation_lengths") is False:
            self.dislocation_lengths = []
        self.dislocation_lengths.extend(lengths)

    def identify_grains(self) -> None:
        """Identify grains in the crystal structure.

        Returns
        -------
        tuple[list, list]
            A tuple containing the orientations and angles of the grains.
        """
        orientations, angles = identify_grain_orientations(self.ovito_pipeline)
        self.orientations = orientations
        self.angles = angles
        # enumerate the orientations and angles
        for i, (orientation, angle) in enumerate(zip(orientations, angles)):
            print(
                f"\033[92mGrain {i + 1}:\033[0m Orientation: {orientation}"
                f"Angle: {round(angle, 4)}"
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
