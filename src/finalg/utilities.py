# ==========
# Utilities
# ==========

import numpy as np
import networkx as nx
import itertools as it
from collections import Counter
import re

def np_arr_to_tuple(arr: np.ndarray) -> tuple:
    """Convert a 2d numpy array into a tuple of tuples to use as a dictionary key.

    EXAMPLE: This function is the 'make_key' input to 'generate_algebra_from_element_dict'.
    """
    return tuple([tuple(row) for row in arr.tolist()])


def delete_row_col(np_arr, row, col):
    """Removes the specified row and col from a Numpy array.
    A new np array is returned, so this does not affect the input array."""
    return np.delete(np.delete(np_arr, row, 0), col, 1)


def get_duplicates(lst):
    """Return a list of the duplicate items in the input list."""
    return [item for item, count in Counter(lst).items() if count > 1]


def all_strings(iterable_object):
    """Checks if all elements in an iterable are strings."""
    return all(isinstance(item, str) for item in iterable_object)


def replace_item(data, old, new):
    """Within a list or tuple, 'data', replace all occurrences
    of 'old' with 'new' while maintaining the original positions."""
    new_data = [new if item == old else item for item in data]
    if isinstance(data, tuple):
        return tuple(new_data)
    else:
        return new_data


def yes_or_no(true_or_false):
    """A convenience function for turning True or False into Yes or No, respectively."""
    if true_or_false:
        return "Yes"
    else:
        return "No"


def symm_diff_of_two_lists_of_lists (list1, list2):
    """Return the symmetric difference between two lists of lists"""
    # Turn the inner lists into tuples, because tuples
    # are hashable and lists are not. Then turn the outer
    # lists into sets, before applying the symm diff operator.
    return set(map(tuple, list1)) ^ set(map(tuple, list2))


def same_lists_of_lists(list_of_lists1, list_of_lists2):
    """Handy for determining, for example, whether a list of left cosets
    is the same as a list of right cosets."""
    return not symm_diff_of_two_lists_of_lists(list_of_lists1, list_of_lists2)


# See https://docs.python.org/3/library/itertools.html#itertools-recipes
def powerset(iterable):
    """Returns the powerset of the input iterable.
    e.g., powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)
    """
    s = list(iterable)
    return it.chain.from_iterable(it.combinations(s, r) for r in range(len(s)+1))


# def remove_singleton_tuple_commas(text):
#     """Python regex code that takes a string containing tuples,
#     and finds and replaces all tuples that consist of a single
#     digit followed by a comma inside parentheses. So, something
#     like '(4, 3), (7,), (1, 2), (5, 6)' would become
#     '(4, 3), (7), (1, 2), (5)'
#     """
#     return re.sub(r'\((\d),\)', r'(\1)', text)


# def print_nx_multidigraph(nx_multidigraph):
#     """Prints nodes and edges of a NetworkX directed graph.
#     Just here for testing code related to Cayley Diagrams."""
#     print("Nodes:", nx_multidigraph.nodes(data=True))
#     print("\nEdges:")
#     for u, v, data in nx_multidigraph.edges(data=True):
#         print(f"  {u} -> {v} : {data}")


def cayley_graph_to_json(cayley_graph):
    """Uses NetworkX to return graph in node-link format that is suitable
    for JSON serialization and use in JavaScript documents."""
    return nx.node_link_data(cayley_graph)


def make_table_from_xml(table_string):
    """This function helps turn the XML-based tables at https://groupprops.subwiki.org/wiki/Main_Page
    into a list of lists for use here.

    Instructions for use:
    1. Copy the table from there and paste it here;
    2. Find & Replace the strings, "<row>" and "</row>", with nothing;
    3. Place triple quotes around the result and give it a variable name;
    4. Then run make_table on the variable.

    Parameters
    ----------
    table_string : str
      XML-based table at https://groupprops.subwiki.org/wiki/Main_Page

    Returns
    -------
    list
      A list of lists of ints, representing a group's multiplication table.
    """
    return [[int(n) for n in row.strip().split(" ")]
            for row in table_string.splitlines()]


def compress_runs(s: str) -> str:
    """
    Replace every run of n identical characters with char^n (n > 1),
    or leave it as a single character when n == 1.
    For example:
        'ffrrrf' -> 'f^2r^3f'
        'aaabbc' -> 'a^3b^2c'
    """
    def replace_run(m: re.Match) -> str:
        char, n = m.group()[0], len(m.group())
        # return char if n == 1 else f'{char}^{n}'
        if n == 1:
            return char
        else:
            return f'{char}^{n}'

    return re.sub(r'(.)\1*', replace_run, s)


# ── Demo ──────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    tests = [
        'ffrrrf',
        'aaabbc',
        'aabbcc',
        'abcd',
        'aaaa',
        'a',
        '',
        'ffrrrfrrf',
    ]
    for t in tests:
        print(f"  {t!r:20s} -> {compress_runs(t)!r}")
# End of File
