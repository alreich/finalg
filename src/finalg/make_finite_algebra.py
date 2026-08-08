# ===========================
#    make_finite_algebra
# ===========================

import json

from finalg.magma import Magma
from finalg.quasigroup_and_loop import Quasigroup, Loop
from finalg.semigroup import Semigroup
from finalg.monoid import Monoid
from finalg.group import Group
from finalg.ring import Ring
from finalg.field import Field, is_field
from finalg.module import Module
from finalg.vector_space import VectorSpace
from finalg.cayley_table import make_cayley_table
from finalg.utilities import get_duplicates, all_strings


def make_finite_algebra(*args):
    """This is the recommended function to use to create any finite algebra.
    It analyzes the input and returns the appropriate finite algebra:
    Group, Ring, Field, VectorSpace, Module, Monoid, Semigroup, or Magma.

    If only 1 input argument, then it must either be a string or a Python
    dictionary.  If it's a string, then it must be a path to a JSON file
    that defines a FiniteAlgebra (i.e., Magma, Semigroup, Monoid,
    Group, Ring, or Field), as described below for the first five arguments.
    If it's a Python dictionary, then it must be the dictionary version of
    such a JSON file. (No JSON or dictionary formats are defined for
    FiniteCompositeAlgebras.)

    Otherwise, the first argument should always be the name (str) of the
    algebra and the second argument should be a description (str) of the
    algebra.

    The remaining arguments depend on whether the algebra being constructed
    is a FiniteAlgebra (i.e, Magma, Semigroup, Monoid, Group,
    Ring, or Field) or a FiniteCompositeAlgebra (i.e., Module or Vector
    Space).

    If constructing a FiniteAlgebra:

    The third argument should be a list of element names (str).

    The fourth argument should be a list of lists of either all integers
    or all strings that represent a finite binary operation.  That is, a
    2-dimensional, square "table" (Cayley table).  The meaning of a table
    entry C corresponding to row A and column B, is that A * B = C, where
    * is the binary operator. If the items in the table are all integers,
    then they must all represent the positions of elements in the element
    list given by the third argument, above. If they are all strings, then
    they must all be members of the list of strings given by the third
    argument.

    A fifth argument is required only if a Ring or Field is being
    constructed, and it should also be a table with structure similar to
    the fourth argument.

    If constructing a FiniteCompositeAlgebra:

    The third argument should be a Ring or Field (the "scalars").

    The fourth argument should be a Group (the "vectors").

    And the fifth argument should be a function that implements the binary
    operation for "scaling vectors".

    See the definitions and examples at https://finalg.readthedocs.io
    """

    if len(args) == 1:

        # Create from a JSON file
        if isinstance(args[0], str):
            with open(args[0], 'r') as fin:
                finalg_dict = json.load(fin)

        # Create from a dictionary
        elif isinstance(args[0], dict):
            finalg_dict = args[0]

        else:
            raise ValueError("If there's a single input, then it must be a string or a dictionary.")

    elif len(args) == 4:

        # The inputs define a Group, Monoid, Semigroup, or Magma.
        # More checks to come farther below.
        finalg_dict = {'name': args[0],
                       'description': args[1],
                       'elements': args[2],
                       'table': args[3]
                       }

    elif len(args) == 5:

        # The inputs define a VectorSpace or Module.
        # It gets created & returned immediately, right here.
        if isinstance(args[3], Group):
            if isinstance(args[2], Field):
                return VectorSpace(args[0], args[1], args[2], args[3], args[4])
            elif isinstance(args[2], Ring):
                return Module(args[0], args[1], args[2], args[3], args[4])
            else:
                raise ValueError(f"{args[2]} must be a Ring or a Field")

        # The inputs define a Field or Ring.
        # More checks to come farther below.
        else:
            finalg_dict = {'name': args[0],
                           'description': args[1],
                           'elements': args[2],
                           'table': args[3],
                           'table2': args[4]
                           }

    elif len(args) == 6:  # Only happens when we're constructing a Cayley-Dickson ring or field

        finalg_dict = {'name': args[0],
                       'description': args[1],
                       'elements': args[2],
                       'table': args[3],
                       'table2': args[4],
                       'conj_map': args[5]  # Lookup table for conjugates
                       }

    else:
        raise ValueError("Incorrect number of input arguments.")

    name = finalg_dict['name']
    desc = finalg_dict['description']
    elems = finalg_dict['elements']
    tbl = finalg_dict['table']

    # Check for duplicate element names
    dups = get_duplicates(elems)
    if len(dups) == 0 and all_strings(elems):
        pass
    else:
        raise ValueError(f"All elements must be unique strings.")

    table = make_cayley_table(tbl, elems)

    # If a second table was input, turn it into a CayleyTable
    # and determine if it supports associativity
    table2 = None
    is_assoc2 = False
    if 'table2' in finalg_dict:
        table2 = make_cayley_table(finalg_dict['table2'], elems)
        is_assoc2 = table2.is_associative()

    # Conjugate Mapping: A lookup table for conjugates in rings or fields that are Cayley-Dickson algebras
    if 'conj_map' in finalg_dict:
        conj_map = finalg_dict['conj_map']
    else:
        conj_map = None

    is_assoc = table.is_associative()
    identity = table.identity()  # this is the integer index of the identity, not the name str
    if identity is not None:
        inverses = table.has_inverses()
    else:
        inverses = None

    # Based on the properties of the inputs, create & return the appropriate algebraic structure
    if is_assoc:
        if identity is not None:
            if inverses:
                if table2 is not None and is_assoc2:
                    # is_field will either build the abelian Group, mentioned in the Field definition,
                    # or it will return False.  In the latter case, this becomes a Ring, instead of a Field.
                    # NOTE: Additional, required checks for multiplicative associativity and distributivity
                    # of multiplication over addition are done within the Field & Ring constructors themselves.
                    abelian_group = is_field(elems[identity], elems, table2.table)
                    if abelian_group:
                        return Field(name, desc, elems, table, table2, check_inputs=False,
                                     mult_sub_grp=abelian_group, conjugate_mapping=conj_map)
                    else:
                        return Ring(name, desc, elems, table, table2, check_inputs=False, conjugate_mapping=conj_map)
                else:
                    return Group(name, desc, elems, table, check_inputs=False)
            else:
                return Monoid(name, desc, elems, table, check_inputs=False)
        else:
            return Semigroup(name, desc, elems, table, check_inputs=False)
    else:
        if table.has_cancellation():
            if table.identity() is not None:
                return Loop(name, desc, elems, table)
            return Quasigroup(name, desc, elems, table)
        else:
            return Magma(name, desc, elems, table)


