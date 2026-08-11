# ================================================================
#   sympy_interop
# ================================================================
#
# Two converters between finalg and sympy.combinatorics, in opposite
# directions:
#
#   from_sympy_permutation_group -- pulls a group definition out of
#     SymPy (sympy.combinatorics.PermutationGroup) and turns it into a
#     finalg Group, so that the very large catalog of groups SymPy
#     already knows how to build or compute (symmetric, alternating,
#     dihedral, cyclic, direct products, derived/commutator subgroups,
#     Sylow subgroups, point/setwise stabilizers, and any custom group
#     built from a handful of permutation generators) can be brought
#     into finalg without hand-authoring a Cayley table.
#
#   to_sympy_permutation_group -- goes the other way: embeds a finalg
#     Group into SymPy as a PermutationGroup via its regular
#     representation (Cayley's theorem), so that a group built or
#     derived inside finalg (by hand, from JSON, as a quotient or
#     subalgebra, etc.) can make use of SymPy's own group-theoretic
#     toolkit. This direction only makes sense for Group -- see that
#     function's docstring for why it doesn't extend to finalg's other
#     algebra types.
#
# Three things make from_sympy_permutation_group fast rather than
# merely correct:
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
from finalg.ring import Ring
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


def to_sympy_permutation_group(finalg_group, verify_order=True):
    """Convert a finalg Group into a sympy.combinatorics.PermutationGroup, via its
    regular representation: Cayley's theorem guarantees that every finite group of
    order n embeds faithfully into Sym(n) this way.

    Concretely, each element g is mapped to the permutation R_g(x) = x * g ("right
    multiplication by g") of the group's own n elements. Associativity of finalg's
    operation makes g -> R_g a homomorphism; because a Group has inverses, it's also
    injective, so the resulting map is a faithful embedding, and the SymPy group's
    own multiplication reproduces this group's table exactly (not just up to
    isomorphism -- literally: for any a, b in the group, R_a * R_b (SymPy
    composition) equals R_(a*b)).

    Note it's *right*, not left, multiplication that works here: SymPy composes
    permutations as (p*q)(i) = q(p(i)), and under that convention the left-regular
    representation turns out to be an anti-homomorphism (it would silently
    reproduce the opposite group's table instead).

    This only makes sense for a genuine finalg Group -- the argument only being a
    single operation with inverses is exactly what makes the embedding a
    homomorphism in the first place. It does not extend to Quasigroup/Loop
    (cancellative, but not associative, so the regular-representation map isn't a
    homomorphism at all and composing the resulting permutations would not
    reproduce the loop's actual product), to Magma/Semigroup in general (without
    cancellation, "multiply by g" need not even be a bijection), or to Ring/Field
    (two operations; SymPy's combinatorics module models single-operation
    permutation groups).

    Every element's permutation is computed (not just a generating set -- SymPy's
    PermutationGroup handles redundant generators fine, and finding a genuinely
    minimal generating set is its own, unrelated combinatorial search), so the
    returned dict covers every element, and the caller can look up the permutation
    for any specific one.

    Parameters
    ----------
    finalg_group : finalg.group.Group
        The group to convert.
    verify_order : bool, default True
        If True (the default -- this check is cheap, unlike the O(n^3) associativity
        check `from_sympy_permutation_group`'s `verify` guards), confirms that the
        SymPy group SymPy computes from the chosen generators has the same order as
        `finalg_group`, and raises ValueError if not.

    Returns
    -------
    (sympy_group, elem_dict)
        sympy_group : the resulting sympy.combinatorics.PermutationGroup.
        elem_dict : dict mapping each of finalg_group's element names to its
            corresponding sympy Permutation (i.e. the full regular representation,
            not just the generators).

    Raises
    ------
    TypeError
        If `finalg_group` isn't a finalg Group.
    ValueError
        If `verify_order=True` and the SymPy group's order doesn't match
        finalg_group's order -- this would indicate a bug rather than a normal
        failure mode, since this embedding is always faithful for a genuine group.
    """
    if not isinstance(finalg_group, Group):
        raise TypeError(
            f"Expected a finalg Group, got {type(finalg_group).__name__}. The regular-"
            f"representation embedding used here relies on the operation being associative "
            f"with inverses, which is exactly what makes a finalg algebra a Group; it isn't "
            f"meaningful for non-associative structures (Magma, Quasigroup, Loop) or "
            f"algebras with two operations (Ring, Field)."
        )
    if isinstance(finalg_group, Ring):  # Ring (and its subclass Field) extend Group via +
        raise TypeError(
            f"{finalg_group.name} is a {type(finalg_group).__name__}, which has two operations "
            f"(addition and multiplication). Converting it as though it were a plain Group would "
            f"silently embed only its additive structure -- self.op -- and throw away multiplication "
            f"entirely, which would be misleading rather than merely partial. If you want the "
            f"additive group on its own, that's unambiguous, so build it directly instead, e.g. "
            f"to_sympy_permutation_group(make_finite_algebra(finalg_group.name, finalg_group.description, "
            f"finalg_group.elements, finalg_group.table.tolist())). If you want the multiplicative "
            f"structure, use finalg_group.units_subgroup() first (its elements form a genuine Group "
            f"under multiplication) and convert that."
        )

    elements = list(finalg_group.elements)
    n = len(elements)
    index_of = {e: i for i, e in enumerate(elements)}

    # Full right-regular representation: R_g(x) = x * g, for every element g (not
    # just generators), so callers can look up any specific element's permutation.
    perm_of = {
        g: Permutation([index_of[finalg_group.op(x, g)] for x in elements])
        for g in elements
    }

    # SymPy's PermutationGroup doesn't need a *minimal* generating set to compute
    # correctly or efficiently (Schreier-Sims handles redundant generators just
    # fine), and finding a genuinely minimal one via finalg's own
    # `_smallest_generating_set` can itself be slow for larger groups -- it's a
    # combinatorial search unrelated to anything this function otherwise needs.
    # Since every element's permutation is already on hand, just hand SymPy all
    # of them.
    sympy_group = PermutationGroup(list(perm_of.values()))

    if verify_order and sympy_group.order() != n:
        raise ValueError(
            f"The SymPy group generated from {finalg_group.name}'s generators has order "
            f"{sympy_group.order()}, but {finalg_group.name} has order {n}. For a genuine "
            f"Group this should never happen and would indicate a bug."
        )

    return sympy_group, perm_of
