"""General utility functions for the project."""

import logging


def convert_plural_to_singular(form: str) -> str:
    """
    Convert if the form is plural and convert it to singular.

    Parameters
    ----------
    form : str
        The form to check

    Returns
    -------
    str
        The singular form
    """
    mapping_plural_singular = {
        "Vacancies": "Vacancy",
        "Interstitials": "Interstitial",
        "Substitutionals": "Substitutional",
    }

    if form not in mapping_plural_singular:
        logging.warning(f"Form '{form}' not found in mapping.")
    return mapping_plural_singular.get(form, form)
