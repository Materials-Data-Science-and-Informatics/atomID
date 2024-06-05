"""Utility functions for comparing RDF graphs."""

from rdflib import Graph, Namespace
from rdflib.compare import graph_diff, to_isomorphic

NS1 = Namespace("http://purls.helmholtz-metadaten.de/cmso/")


def find_unique_identifier(graph):
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


def substitute_uuid_in_file(filepath, current_uuid, replacement_uuid):
    """Substitute the current UUID with the replacement UUID in the file."""
    with open(filepath, "r") as file:
        filedata = file.read().replace(current_uuid, replacement_uuid)
    return filedata


def get_substitute_path(filepath):
    """Get the path for the substituted file."""
    return "".join(filepath.split(".")[:-1]) + "_sub.ttl"


def load_and_parse_graph(filepath):
    """Load and parse the graph from the given file."""
    graph = Graph()
    graph.parse(filepath, format="turtle")
    return graph


def graph_difference(substitute_graph, reference_graph):
    """Compute the differences between two graphs."""
    iso_substitute_graph = to_isomorphic(substitute_graph)
    iso_reference_graph = to_isomorphic(reference_graph)
    return graph_diff(iso_substitute_graph, iso_reference_graph)


def compare_graphs(current_filepath, reference_filepath):
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


def handle_graph_differences(in_first, in_second):
    """Handle the differences found in graph comparison."""
    compare_results = []
    for s, p, o in in_first:
        if "hasPath" in str(p) or "hasIdentifier" in str(p):
            continue  # Skip known irrelevant differences

        query = f"SELECT ?object WHERE {{ <{s}> <{p}> ?object .}}"
        results = in_second.query(query)
        if not results:
            compare_results.append(handle_no_matches(s, p, in_second))
        else:
            compare_results.append(handle_result_differences(results, o, s, p))
    return [res for res in compare_results if res]  # Filter out None results


def handle_no_matches(s, p, in_second):
    """Handle cases where no matching subjects were found in the reference graph."""
    query = f"SELECT ?object WHERE {{ <{s}> ?predicate ?object .}}"
    results_sub = in_second.query(query)
    if not results_sub:
        return f"No match found in reference graph for: {s}"
    return f"Match not found for predicate {p} with subject {s}"


def handle_result_differences(results, original_object, s, p):
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
    return f"No significant differences for: {s}, {p}"
