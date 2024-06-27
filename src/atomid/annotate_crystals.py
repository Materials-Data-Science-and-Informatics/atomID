"""AnnotateCrystals class."""

from typing import Dict, List, Optional, TypedDict

import atomrdf as ardf

from atomid.annotate import AnnotateCrystal


class CrystalEntry(TypedDict):
    """Type definition for a crystal entry."""

    data_file: str
    system: AnnotateCrystal


class AnnotateCrystals:
    """Annotate multiple crystal structures."""

    def __init__(self, data_files: Optional[List[str]] = None):
        self._kg = ardf.KnowledgeGraph()
        self.crystals_dict: Dict[int, CrystalEntry] = {}

        if data_files:
            for idx, data_file in enumerate(data_files, start=1):
                self.crystals_dict[idx] = {
                    "data_file": data_file,
                    "system": AnnotateCrystal(data_file, kg=self._kg),
                }

    @property
    def num_samples(self) -> int:
        """Return the number of samples."""
        return len(self.crystals_dict)

    @property
    def kg(self) -> ardf.KnowledgeGraph:
        """Return the knowledge graph."""
        return self._kg

    def add_data_file(self, data_file: str) -> None:
        """Add a data file to the list of samples.

        Parameters
        ----------
        data_file : str
            The path to the data file.

        Returns
        -------
        None
        """
        idx = self.num_samples + 1
        self.crystals_dict[idx] = {
            "data_file": data_file,
            "system": AnnotateCrystal(data_file, kg=self._kg),
        }

    def get_crystal(self, idx: int) -> AnnotateCrystal:
        """Return the AnnotateCrystal object for a given sample.

        Parameters
        ----------
        idx : int
            The index of the sample.

        Returns
        -------
        AnnotateCrystal
            The AnnotateCrystal object for the sample.
        """
        return self.crystals_dict[idx]["system"]

    def get_data_file(self, idx: int) -> str:
        """Return the data file for a given sample.

        Parameters
        ----------
        idx : int
            The index of the sample.

        Returns
        -------
        str
            The path to the data file.
        """
        return self.crystals_dict[idx]["data_file"]

    def get_sample(self, idx: int) -> CrystalEntry:
        """Return the dictionary for a given sample.

        Parameters
        ----------
        idx : int
            The index of the sample.

        Returns
        -------
        dict
            The dictionary for the sample.
        """
        return self.crystals_dict[idx]

    def annotate_all_crystal_structures(self) -> None:
        """Annotate the crystal structures for all samples.

        This method identifies the crystal structure for each sample and annotates it.

        Returns
        -------
        None
        """
        for idx in self.crystals_dict:
            self.crystals_dict[idx]["system"].identify_crystal_structure()
            if self.crystals_dict[idx]["system"].crystal_type != "other":
                self.crystals_dict[idx]["system"].annotate_crystal_structure()
            else:
                print(
                    f"Crystal structure not identified for sample {idx}"
                    f"{self.crystals_dict[idx]['data_file']}"
                )
                print(
                    "Manually annotate the crystal structure using the "
                    "set_crystal_structure method."
                )

    def annotate_crystal_structure(self, idx: int, log: bool = False) -> None:
        """Annotate the crystal structure for a given sample.

        This method identifies the crystal structure for a given sample and annotates it.

        Parameters
        ----------
        idx : int
            The index of the sample.
        log : bool, optional
            Whether to log the output. Defaults to False.

        Returns
        -------
        None
        """
        self.crystals_dict[idx]["system"].identify_crystal_structure(log=log)
        if self.crystals_dict[idx]["system"].crystal_type is not None:
            self.crystals_dict[idx]["system"].annotate_crystal_structure()
        else:
            print(
                f"Crystal structure not identified for sample {idx}"
                f"{self.crystals_dict[idx]['data_file']}"
            )
            print(
                "Manually annotate the crystal structure using the "
                "set_crystal_structure method."
            )

    def identify_point_defects_all_samples(
        self,
        reference_files: List[str],
        reference_format_list: Optional[List[str]] = None,
    ) -> None:
        """Identify point defects in all samples.

        Parameters
        ----------
        reference_files : List[str]
            The list of reference files.
        reference_format_list : List[str], optional
            The list of reference file formats.


        Returns
        -------
        None
        """
        print(self.crystals_dict)
        for idx in self.crystals_dict:
            print(idx)
            self.identify_point_defects(
                idx,
                reference_files[idx - 1],
                reference_format_list[idx - 1] if reference_format_list else None,
            )

    def identify_point_defects(
        self, idx: int, reference_file: str, reference_format: Optional[str] = None
    ) -> None:
        """Identify point defects in a given sample.

        Parameters
        ----------
        idx : int
            The index of the sample.
        reference_file : str
            The reference file.
        reference_format : str, optional
            The format of the reference file.

        Returns
        -------
        None
        """
        self.crystals_dict[idx]["system"].identify_point_defects(
            reference_file, ref_format=reference_format
        )

    def annotate_point_defects_all_samples(self) -> None:
        """Annotate point defects in all samples.

        Returns
        -------
        None
        """
        for idx in self.crystals_dict:
            self.annotate_point_defects(idx)

    def annotate_point_defects(self, idx: int) -> None:
        """Annotate point defects in a given sample.

        Parameters
        ----------
        idx : int
            The index of the sample.

        Returns
        -------
        None
        """
        self.crystals_dict[idx]["system"].annotate_point_defects()
