"""Grain segmentation analysis module."""

import numpy as np
from ovito.modifiers import GrainSegmentationModifier
from ovito.pipeline import Pipeline


def identify_grain_orientations(input_pipeline: Pipeline) -> tuple:
    """
    Identify grains in the crystal structure.

    Parameters
    ----------
    input_pipeline : Pipeline
        The pipeline object that handles the input data.

    Returns
    -------
    tuple
        A tuple containing the grain orientation axes and angles.
    """
    grain_modifier = GrainSegmentationModifier()
    input_pipeline.modifiers.append(grain_modifier)
    grain_data = input_pipeline.compute()

    orientations = grain_data.tables["grains"]["Orientation"]
    axes, angles = get_grain_orientations(orientations)

    return axes, angles


def get_grain_orientations(orientations: np.ndarray) -> tuple:
    """
    Get grain orientations from quaternion representations.

    Parameters
    ----------
    orientations : list of list of float
        List of quaternions representing the grain orientations.

    Returns
    -------
    tuple
        A tuple containing arrays of axes and angles for each grain.
    """
    quaternions = np.array(orientations)
    qx, qy, qz, qw = quaternions.T

    scale = np.sqrt(qx**2 + qy**2 + qz**2)
    axes = np.zeros((len(quaternions), 3))
    valid_scale = scale > 1e-12

    axes[valid_scale, 0] = qx[valid_scale] / scale[valid_scale]
    axes[valid_scale, 1] = qy[valid_scale] / scale[valid_scale]
    axes[valid_scale, 2] = qz[valid_scale] / scale[valid_scale]

    axes = np.round(axes, decimals=4)

    angles = np.zeros(len(quaternions))
    angles[valid_scale] = 2.0 * np.arccos(qw[valid_scale])

    return axes, angles
