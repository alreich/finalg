# ==========
# Utilities
# ==========

import numpy as np
import networkx as nx
import itertools as it
import json
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


def is_power_of_two(n):
    return n > 0 and (n & (n - 1)) == 0


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


# def cayley_graph_to_json(cayley_graph):
#     """Uses NetworkX to return graph in node-link format that is suitable
#     for JSON serialization and use in JavaScript documents."""
#     return nx.node_link_data(cayley_graph)


def cayley_graph_to_json(cayley_graph):
    """Uses NetworkX to return graph in node-link format that is suitable
    for JSON serialization and use in JavaScript documents."""
    return nx.node_link_data(cayley_graph, edges="edges")


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


# ── Helpers for building groups and writing them as pretty-printed JSON ────────
# (Used, e.g., by the notebooks in notebooks/small_groups.)

def mult_action(modulus, k):
    """The action of Z_m on Z_modulus (elements are strings of integers) in which h acts by
    x -> k^h * x (mod modulus). Use with ``Group.semidirect_product``."""
    return lambda h, x: str((int(x) * pow(k, int(h), modulus)) % modulus)


def _extend_hom(N, gen_images):
    """Extend {generator: image} to a map on all of N by f(a*g) = f(a)*f(g)."""
    if not gen_images:
        return {n: n for n in N.elements}                 # no images given: the identity automorphism
    f = {N.identity: N.identity}
    f.update(gen_images)
    queue = list(f)
    while queue:
        a = queue.pop()
        for g, img in gen_images.items():
            p = N.op(a, g)
            if p not in f:
                f[p] = N.op(f[a], img)
                queue.append(p)
    return f


def action_from_generators(N, H, images):
    """Build an action of the group H on the group N by automorphisms, from generators.

    `images` maps each of a set of generators of H to a dictionary {n_gen: image, ...} that says where the
    corresponding automorphism of N sends a set of generators of N (an empty dictionary means that the
    H-generator acts trivially). The result is a dictionary {h: {n: phi_h(n)}} that covers every element
    h of H and n of N, ready to be passed to ``N.semidirect_product(H, action)``, which verifies that it
    really is a homomorphism from H into the automorphisms of N."""
    identity_map = {n: n for n in N.elements}
    gen_autos = {h: _extend_hom(N, imgs) for h, imgs in images.items()}
    phi = {H.identity: identity_map}
    queue = [H.identity]
    while queue:
        h = queue.pop()
        for g, auto in gen_autos.items():
            hg = H.op(h, g)
            if hg not in phi:
                phi[hg] = {n: phi[h][auto[n]] for n in N.elements}     # phi_(hg) = phi_h o phi_g
                queue.append(hg)
    if len(phi) != H.order:
        raise ValueError(f"The generators given for {H.name} generate only {len(phi)} of its {H.order} elements.")
    return phi


# ----------------------------------------------------------------------------------------
#  Pretty JSON
# ----------------------------------------------------------------------------------------

def _wrapped_list(items, indent, width):
    """JSON-encode a list of scalars on as few lines as fit within `width` columns."""
    pieces = [json.dumps(x, ensure_ascii=False) for x in items]
    one_line = '[' + ', '.join(pieces) + ']'
    if indent + len(one_line) <= width:
        return one_line
    lines, current = [], ''
    for p in pieces:
        candidate = (current + ', ' + p) if current else p
        if current and indent + 2 + len(candidate) + 1 > width:
            lines.append(current + ',')
            current = p
        else:
            current = candidate
    lines.append(current)
    pad = ' ' * (indent + 2)
    return '[\n' + '\n'.join(pad + ln for ln in lines) + '\n' + ' ' * indent + ']'


def pretty_json(d, width=100):
    """Return the JSON text for an algebra dictionary (from ``alg.to_dict()``), laid out like the example
    files in finalg/data/algebras: one key per line, the elements wrapped to `width` columns, and one
    row of each Cayley table per line."""
    parts = []
    for key, value in d.items():
        k = json.dumps(key)
        if key in ('table', 'table2'):
            rows = ',\n'.join('    ' + json.dumps(list(r)) for r in value)
            parts.append(f'  {k}: [\n{rows}\n  ]')
        elif key == 'elements':
            parts.append(f'  {k}: ' + _wrapped_list(list(value), 2, width))
        else:
            parts.append(f'  {k}: ' + json.dumps(value, ensure_ascii=False))
    return '{\n' + ',\n'.join(parts) + '\n}\n'


def dump_pretty(alg, path, width=100):
    """Write an algebra to `path` as pretty-printed JSON that ``make_finite_algebra(path)`` can read back."""
    with open(path, 'w', encoding='utf-8') as fout:
        fout.write(pretty_json(alg.to_dict(), width))


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
