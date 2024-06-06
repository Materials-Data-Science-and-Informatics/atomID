from ovito.modifiers import DislocationAnalysisModifier


def identify_dislocations(data_pipeline) -> None:
    """Identify dislocations in the crystal structure.

    Parameters
    ----------
    data_file : str
        The name of the file to read
    """
    dislocation_modifier = DislocationAnalysisModifier()
    data_pipeline.modifiers.append(dislocation_modifier)
    discloation_data = data_pipeline.compute()

    burgers_vectors = []
    lengths = []
    for lines in discloation_data.dislocations.lines:
        burgers_vectors.append(lines.true_burgers_vector)
        lengths.append(lines.length)

    return burgers_vectors, lengths