__all__ = ["make_finite_algebra",
           "generate_cyclic_group",
           "generate_symmetric_group",
           "generate_powerset_group",
           "generate_commutative_monoid",
           "generate_relative_primes_group",
           "generate_powerset_ring",
           "generate_algebra_mod_n",
           "generate_nxn_matrix_algebra",
           "generate_dihedral_group",
           "generate_algebra_from_element_dict",
           "make_cayley_table",
           "about_tables",
           "InfixNotation",
           "Perm",
           "examples",
           "about_subalgebras",
           "find_isomorphic_subalgebra",
           "np_arr_to_tuple",
           "cayley_graph_to_json",
           "compress_runs",
           "from_sympy_permutation_group"
]

from finalg.make_finite_algebra import make_finite_algebra
from finalg.cayley_table import make_cayley_table, about_tables
from finalg.examples import Examples
from finalg.infix_notation import InfixNotation
from finalg.about_subalgebras import about_subalgebras, find_isomorphic_subalgebra
from finalg.algebra_generators import generate_cyclic_group, generate_symmetric_group, generate_powerset_group,\
    generate_commutative_monoid, generate_relative_primes_group, generate_powerset_ring, generate_algebra_mod_n,\
    generate_nxn_matrix_algebra, generate_algebra_from_element_dict, generate_dihedral_group
from finalg.utilities import np_arr_to_tuple, cayley_graph_to_json, compress_runs
from finalg.permutation import Perm
from finalg.sympy_interop import from_sympy_permutation_group

from importlib import resources

# Determine alg_dir, the directory that contains example algebras.
# By default, alg_dir contains a file, examples.json, that lists
# file names of example algebras.
#
# This uses importlib.resources rather than a path relative to __file__,
# so it works correctly whether finalg is used from a source checkout,
# an editable install, or a wheel installed from PyPI (where there is no
# "sibling" algebras/ directory next to the installed package).
alg_dir = resources.files("finalg.data.algebras")
example_file_names = "examples.json"
file_names_path = alg_dir / example_file_names

if file_names_path.is_file():
    if __name__ != "__main__":
        print("Finite Algebras:")
        print(f"  To see a listing of built-in example algebras, run {__name__}.examples.about()")
        print(f"  To retrieve an example algebra, use {__name__}.examples[INDEX]\n")
        examples = Examples(alg_dir, example_file_names)
else:
    print("\nExamples not found.\n")
