# ================================================================
#   sympy_interop
# ================================================================
#
# A converter that pulls a group definition out of SymPy's
# combinatorics library (sympy.combinatorics.PermutationGroup) and
# turns it into a finalg Group, so that the very large catalog of
# groups SymPy already knows how to build or compute (symmetric,
# alternating, dihedral, cyclic, direct products, derived/commutator
# subgroups, Sylow subgroups, point/setwise stabilizers, and any
# custom group built from a handful of permutation generators) can be
# brought into finalg without hand-authoring a Cayley table.
#
# Three things make this fast rather than merely correct:
#
# 1. The group's elements are enumerated by SymPy itself
#    (`PermutationGroup.generate()`), which uses SymPy's own (much
#    more efficient) group-theoretic algorithms, rather than being
#    rediscovered with finalg's generic closure algorithm (repeatedly
#    multiplying known elements together until nothing new appears --
#    fine for small algebras, but quadratic-or-worse per pass, and
#    unaware of any of the group structure SymPy already knows).
#
# 2. Once the element list is known, building the n x n Cayley table
#    only needs *some* correct, associative composition -- it doesn't
#    need SymPy's own Permutation objects, which carry a fair amount
#    of per-call overhead for validity checks, cycle bookkeeping, etc.
#    SymPy defines (p*q)(i) = q(p(i)); precomputing each permutation's
#    array form once and then composing with plain tuple indexing
#    reproduces exactly the same table using ordinary Python integers.
#
# 3. finalg's usual entry point, `make_finite_algebra`, always runs an
#    O(n^3) pure-Python associativity check to decide what class of
#    algebra to build. A `PermutationGroup` is associative and fully
#    invertible by construction (that's what makes it a group), so
#    this builds a `Group` directly instead and skips that redundant
#    check (see the `verify` parameter to opt back into it).
#
# Together, these mean a several-hundred-element group -- e.g. all of
# S6, order 720 -- converts in well under a second, where routing
# through the usual construction path would take minutes.
#
# finalg's own `Perm` class is used only to produce readable
# cycle-notation element names; it plays no part in the actual group
# operation.

from sympy.combinatorics import Permutation, PermutationGroup

from finalg.group import Group
from finalg.permutation import Perm as FinalgPerm


def from_sympy_permutation_group(sympy_group, name=None, description=None, verify=False):
    """Convert a sympy.combinatorics.PermutationGroup into a finalg Group.

    Works with anything SymPy can hand you as a PermutationGroup: the named
    constructors (SymmetricGroup, AlternatingGroup, DihedralGroup, CyclicGroup,
    AbelianGroup, ...), groups built directly from a handful of Permutation
    generators, or groups SymPy derives from another group (e.g.
    `.derived_subgroup()`, `.sylow_subgroup(p)`, a point stabilizer, a direct
    product via `DirectProduct(...)`).

    The group's elements are enumerated by SymPy itself, and the O(n^2) Cayley
    table is filled in using plain-tuple permutation composition rather than
    SymPy's own (comparatively heavyweight) Permutation objects. The Group is
    then built directly (bypassing finalg's usual `make_finite_algebra`
    dispatcher, which would otherwise re-verify associativity with a pure-Python
    O(n^3) triple loop over the whole table) -- a `PermutationGroup` is
    associative and invertible by construction, so that check is redundant
    here. Net effect: even a few-hundred-element group, e.g. all of S6 (order
    720), converts in a fraction of a second, where going through the usual
    O(n^3) verification path would take minutes. It's still an O(n^2) table,
    though, so this isn't the tool for something with tens of thousands of
    elements or more.

    Parameters
    ----------
    sympy_group : sympy.combinatorics.PermutationGroup
        The group to convert.
    name : str, optional
        Name for the resulting finalg algebra. Defaults to a name built from
        the group's order and degree.
    description : str, optional
        Description for the resulting algebra. Defaults to a note that it was
        converted from a SymPy PermutationGroup, together with the group's
        order and degree.
    verify : bool, default False
        If True, explicitly re-checks associativity and inverses on the
        resulting table (the same O(n^3) work `make_finite_algebra` would
        otherwise do) and raises ValueError if either check fails. Off by
        default for speed; useful as a one-off sanity check on smaller
        groups, or while testing changes to this function itself.

    Returns
    -------
    (finalg_group, elem_dict)
        finalg_group : the resulting finalg Group.
        elem_dict : dict mapping each finalg element name (str, in cycle
            notation) to the corresponding sympy Permutation (aligned to the
            group's degree), in case the caller wants to work with the
            original SymPy objects.

    Raises
    ------
    TypeError
        If `sympy_group` isn't a sympy.combinatorics.PermutationGroup.
    ValueError
        If the number of elements SymPy enumerates doesn't match the order it
        reports for the group, if two distinct elements were assigned the
        same name (either would indicate a bug rather than a normal failure
        mode), or if `verify=True` and the resulting table fails an
        associativity or inverses check.
    """
    if not isinstance(sympy_group, PermutationGroup):
        raise TypeError(
            f"Expected a sympy.combinatorics.PermutationGroup, got {type(sympy_group).__name__}."
        )

    degree = sympy_group.degree
    expected_order = sympy_group.order()

    # Enumerate elements via SymPy's own algorithms, then drop down to plain
    # tuples of ints (each permutation's array form, aligned to `degree`) for
    # the O(n^2) work below -- SymPy's Permutation objects are correct but
    # comparatively expensive to construct and multiply at this scale.
    sympy_perms = [Permutation(p.array_form, size=degree) for p in sympy_group.generate()]
    if len(sympy_perms) != expected_order:
        raise ValueError(
            f"SymPy enumerated {len(sympy_perms)} element(s) via generate(), but reports "
            f"this group's order as {expected_order}."
        )

    arrays = [tuple(p.array_form) for p in sympy_perms]
    names = [str(FinalgPerm(list(a))) for a in arrays]
    if len(set(names)) != len(names):
        raise ValueError("Two distinct group elements were assigned the same name; this is a bug.")

    index_of = {a: i for i, a in enumerate(arrays)}

    # SymPy defines (p*q)(i) = q(p(i)), i.e. (p*q).array_form[i] == q[p[i]].
    table = [[index_of[tuple(q[x] for x in p)] for q in arrays] for p in arrays]

    alg_name = name or f"SymPy_Group_order{expected_order}_deg{degree}"
    alg_description = description or (
        f"Converted from a sympy.combinatorics.PermutationGroup "
        f"of order {expected_order} acting on {degree} points."
    )

    # A PermutationGroup is associative and fully invertible by construction,
    # so build the Group directly rather than routing through
    # make_finite_algebra (which would otherwise re-derive that fact with an
    # O(n^3) table scan).
    algebra = Group(alg_name, alg_description, names, table, check_inputs=False)

    if verify:
        if not algebra.is_associative():
            raise ValueError("Converted table failed an associativity check (verify=True).")
        if not algebra.has_inverses():
            raise ValueError("Converted table failed an inverses check (verify=True).")

    elem_dict = dict(zip(names, sympy_perms))

    return algebra, elem_dict
