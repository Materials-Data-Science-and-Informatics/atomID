"""Utility functions for comparing RDF graphs."""

from typing import Any

from rdflib import Graph, Namespace, URIRef
from rdflib.compare import graph_diff, to_isomorphic

NS1 = Namespace("http://purls.helmholtz-metadaten.de/cmso/")


def find_unique_identifier(graph: Graph) -> Any:
    """Find the unique identifier in the graph."""
    query = """
        PREFIX ns1: <http://purls.helmholtz-metadaten.de/cmso/>
        SELECT ?subject WHERE {
            ?subject ?predicate ns1:AtomicScaleSample .
        }
    """
    results = list(graph.query(query))
    if len(results) > 1:
        raise ValueError("Multiple unique identifiers are not supported")
    return results[0].subject.split(":")[1] if results else None


def substitute_uuid_in_file(
    filepath: str, current_uuid: str, replacement_uuid: str
) -> str:
    """Substitute the current UUID with the replacement UUID in the file."""
    with open(filepath, "r") as file:
        filedata = file.read().replace(current_uuid, replacement_uuid)
    return filedata


def get_substitute_path(filepath: str) -> str:
    """Get the path for the substituted file."""
    return "".join(filepath.split(".")[:-1]) + "_sub.ttl"


def load_and_parse_graph(filepath: str) -> Graph:
    """Load and parse the graph from the given file."""
    graph = Graph()
    graph.parse(filepath, format="turtle")
    return graph


def graph_difference(substitute_graph: Graph, reference_graph: Graph) -> Any:
    """Compute the differences between two graphs."""
    iso_substitute_graph = to_isomorphic(substitute_graph)
    iso_reference_graph = to_isomorphic(reference_graph)
    return graph_diff(iso_substitute_graph, iso_reference_graph)


def compare_graphs(current_filepath: str, reference_filepath: str) -> tuple[bool, list]:
    """Compare two RDF graphs given their file paths and return True if they are the same, otherwise False."""
    current_graph = load_and_parse_graph(current_filepath)
    reference_graph = load_and_parse_graph(reference_filepath)

    current_uuid = find_unique_identifier(current_graph)
    reference_uuid = find_unique_identifier(reference_graph)

    filedata = substitute_uuid_in_file(current_filepath, current_uuid, reference_uuid)

    substituted_graph = Graph()
    substituted_graph.parse(data=filedata, format="turtle")

    in_both, in_first, in_second = graph_difference(substituted_graph, reference_graph)

    differences = handle_graph_differences(in_first, in_second)

    return len(differences) == 0, differences


def handle_graph_differences(in_first: Graph, in_second: Graph) -> list:
    """
    Handle the differences found in graph comparison.

    This function compares two RDF graphs and identifies differences.

    Args:
        in_first (Graph): The first RDF graph.
        in_second (Graph): The second RDF graph.

    Returns
        list: A list of differences found between the two graphs.
    """
    compare_results = []

    for subject, predicate, obj in in_first:
        if "hasPath" in str(predicate) or "hasIdentifier" in str(predicate):
            continue  # Skip known irrelevant differences

        query = f"SELECT ?object WHERE {{ <{subject}> <{predicate}> ?object .}}"
        results = in_second.query(query)
        if not results:
            compare_results.append(handle_no_matches(subject, predicate, in_second))
        else:
            compare_results.append(
                handle_result_differences(results, obj, subject, predicate)
            )

    return [result for result in compare_results if result]


def handle_no_matches(s: URIRef, p: URIRef, in_second: Graph) -> str:
    """Handle cases where no matching subjects were found in the reference graph."""
    query = f"SELECT ?object WHERE {{ <{s}> ?predicate ?object .}}"
    results_sub = in_second.query(query)
    if not results_sub:
        return f"No match found in reference graph for: {s}"
    return f"Match not found for predicate {p} with subject {s}"


def handle_result_differences(
    results: list, original_object: URIRef, s: URIRef, p: URIRef
) -> str:
    """Evaluate and print differences in query results."""

    for row in results:
        try:
            if (
                abs(float(row.object) - float(original_object)) / float(row.object)
                > 0.1
            ):
                return f"Values differ significantly for: {s}, {p}"
        except ValueError:
            return f"Non-numeric comparison for: {s}, {p}\nOriginal: {original_object}, Found: {row.object}"

    return ""
