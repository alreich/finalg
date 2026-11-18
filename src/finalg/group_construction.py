# ==============================================================
#   Group construction: automorphisms and cyclic extensions
# ==============================================================
"""Build groups from structure, instead of searching through all possible Cayley tables.

Every solvable group G has a normal subgroup N of prime index p. If t is any element of G outside N,
then conjugation by t is an automorphism alpha of N, and c = t^p lies in N, with

    alpha(c) = c    and    alpha^p = conjugation by c.

Conversely, every pair (alpha, c) that satisfies those two conditions defines a group of order p * |N|,
whose elements are x t^i (x in N, 0 <= i < p), with the product

    (x t^i)(y t^j) = x alpha^i(y) t^(i+j),    where t^p = c.

If c is the identity, the result is the semidirect product of N by Z_p (see ``Group.semidirect_product``).
Otherwise, it is a non-split extension (e.g., Q8, Z8, and the dicyclic groups). The same holds with any
m >= 2 in place of p, and then the quotient G/N is cyclic of order m.

The functions are:

* ``automorphisms(G)``, ``inner_automorphisms(G)`` and ``automorphism_group(G)`` compute the automorphisms
  of a Group. An automorphism is represented by a tuple ``f``, where ``f[i]`` is the index, in
  ``G.elements``, of the image of the element with index ``i``.
* ``cyclic_extensions(N, m)`` returns all groups, up to isomorphism, that contain N as a normal subgroup
  with a cyclic quotient of order m. Its helpers are ``valid_extension_data``, ``reduce_by_conjugacy`` and
  ``build_cyclic_extension``.
* ``groups_of_order(n)`` returns all *solvable* groups of order n, up to isomorphism, by taking cyclic
  extensions by Z_p of all the groups of order n/p, for every prime p that divides n. Every group of order
  less than 60 is solvable, so, below 60, this is all of the groups of that order.

Isomorphism is always tested with ``fast_isomorphic``, after first sorting the candidates into buckets by a
cheap isomorphism invariant.
"""

import itertools as it
import warnings
from collections import defaultdict

import numpy as np
from sympy import factorint

from finalg.group import Group
from finalg.make_finite_algebra import make_finite_algebra

__all__ = ["automorphisms",
           "inner_automorphisms",
           "automorphism_group",
           "valid_extension_data",
           "reduce_by_conjugacy",
           "build_cyclic_extension",
           "cyclic_extensions",
           "relabel_group",
           "groups_of_order"]


# ----------------------------------------------------------------------------------------
#  Integer views of a Group
# ----------------------------------------------------------------------------------------

def _int_view(G):
    """Return (T, e), where T is the Cayley table of G as a NumPy array of int, and e is the index of
    the identity element."""
    T = np.asarray(G.table.table, dtype=np.int64)
    e = list(G.elements).index(G.identity)
    return T, e


def _inverse_array(T, e):
    """Return the array whose i-th entry is the index of the inverse of element i."""
    inv = np.empty(len(T), dtype=np.int64)
    rows, cols = np.where(T == e)
    inv[rows] = cols
    return inv


def _element_orders(T, e):
    """Return the array whose i-th entry is the order of element i."""
    n = len(T)
    x = np.arange(n)
    power = x.copy()
    order = np.ones(n, dtype=np.int64)
    todo = power != e
    while todo.any():
        power[todo] = T[power[todo], x[todo]]
        order[todo] += 1
        todo = power != e
    return order


def _invariant(T, e):
    """An isomorphism invariant that is cheap to compute: the sorted list, over the elements x, of
    (order of x, size of the centralizer of x, order of x^2, number of y such that y^2 = x)."""
    n = len(T)
    order = _element_orders(T, e)
    centralizer_size = (T == T.T).sum(axis=1)
    squares = T[np.arange(n), np.arange(n)]
    num_roots = np.bincount(squares, minlength=n)
    return tuple(sorted(zip(order.tolist(), centralizer_size.tolist(),
                            order[squares].tolist(), num_roots.tolist())))


# ----------------------------------------------------------------------------------------
#  Automorphisms
# ----------------------------------------------------------------------------------------

def _extend(table_rows, e, n, gens, images):
    """Extend the assignment gens[i] -> images[i] to a map f such that f(a * g) = f(a) * f(g), for every
    element a and every generator g. Return f as a list of indices, or None if the assignment is
    inconsistent or the result is not a bijection.

    Every edge a -> a * g is checked, so an assignment that survives is a homomorphism (by induction on
    the length of a product of generators), and, being a bijection, an automorphism."""
    f = [-1] * n
    f[e] = e
    for g, image in zip(gens, images):
        f[g] = image
    pairs = list(zip(gens, images))
    queue, seen = [e], {e}
    for a in queue:
        row, f_row = table_rows[a], table_rows[f[a]]
        for g, image in pairs:
            p, value = row[g], f_row[image]
            if f[p] == -1:
                f[p] = value
            elif f[p] != value:
                return None
            if p not in seen:
                seen.add(p)
                queue.append(p)
    if len(seen) != n or len(set(f)) != n:
        return None
    return f


def automorphisms(G):
    """Return a list of all of the automorphisms of the Group G, each one a tuple f, where f[i] is
    the index, in G.elements, of the image of the element with index i.

    An automorphism is determined by its values on a generating set, so only the assignments of
    images to the elements of a small generating set are tried, and, for each generator, only
    the elements that have the same invariants (order, row and column profiles, size of centralizer) as it
    does. This is fast if the number of automorphisms is up to the tens of thousands (GL(4,2) has 20160),
    but it is not meant for groups like Z2^5, which has almost 10 million."""
    T, e = _int_view(G)
    n = len(T)
    if n == 1:
        return [(0,)]
    names = list(G.elements)
    gens = [names.index(g) for g in G._smallest_generating_set()]
    invariants = G._element_invariants()
    with_invariant = defaultdict(list)
    for i, name in enumerate(names):
        with_invariant[invariants[name]].append(i)
    candidates = [with_invariant[invariants[names[g]]] for g in gens]
    table_rows = T.tolist()
    result = []
    for images in it.product(*candidates):
        if len(set(images)) < len(images):
            continue
        f = _extend(table_rows, e, n, gens, images)
        if f is not None:
            result.append(tuple(f))
    return result


def _inner_array(T, e):
    """Row x of the result is the inner automorphism a -> x * a * x^-1 (as an array of indices)."""
    n = len(T)
    inv, indices = _inverse_array(T, e), np.arange(n)
    return np.array([T[T[x, indices], inv[x]] for x in range(n)], dtype=np.int64)


def inner_automorphisms(G):
    """Return a sorted list of the distinct inner automorphisms of the Group G, as tuples of indices."""
    T, e = _int_view(G)
    return sorted({tuple(row) for row in _inner_array(T, e).tolist()})


def automorphism_group(G, max_order=2000, return_maps=False):
    """Return the automorphism group of the Group G, as a Group whose elements are named 'f0' (the
    identity), 'f1', 'f2', ..., and whose product is composition: f_i * f_j = f_i o f_j.

    If return_maps is True, then the list of automorphisms (tuples of indices, as returned by
    `automorphisms`) that the names stand for is returned as well, as a second value, so that 'f<k>' is
    the k-th item in the list.

    A ValueError is raised if the automorphism group has more than max_order elements, because its
    Cayley table would be too big."""
    auts = automorphisms(G)
    k, n = len(auts), G.order
    if k > max_order:
        raise ValueError(f"|Aut({G.name})| = {k} > max_order = {max_order}")
    identity = tuple(range(n))
    auts.sort(key=lambda f: f != identity)  # The identity goes first
    A = np.array(auts, dtype=np.int64)
    index = {A[i].tobytes(): i for i in range(k)}
    table = np.empty((k, k), dtype=np.int64)
    for i in range(k):  # table[i][j] is the index of f_i o f_j
        table[i] = [index[row.tobytes()] for row in A[i][A]]
    aut_group = make_finite_algebra(f"Aut({G.name})", f"Automorphism group of {G.name}",
                                    [f"f{i}" for i in range(k)], table.tolist())
    return (aut_group, auts) if return_maps else aut_group


# ----------------------------------------------------------------------------------------
#  Cyclic extensions
# ----------------------------------------------------------------------------------------

def valid_extension_data(N, m, auts=None):
    """Return all pairs (alpha, c), where alpha is an automorphism of the Group N (a tuple of indices)
    and c is the index of an element of N, such that alpha(c) = c and alpha^m is conjugation by c.

    These are exactly the data that define a group G = <N, t>, which contains N as a normal subgroup, with
    G/N cyclic of order m, where t^m = c and t x t^-1 = alpha(x) for every x in N. If the automorphisms
    of N have already been computed, they can be passed in as `auts`."""
    T, e = _int_view(N)
    n = len(T)
    auts = automorphisms(N) if auts is None else auts
    inner = _inner_array(T, e)
    c_for_inner = defaultdict(list)  # conjugation by c, as bytes -> all of the c that give it
    for c in range(n):
        c_for_inner[inner[c].tobytes()].append(c)
    data = []
    for f in auts:
        F = np.array(f, dtype=np.int64)
        f_to_the_m = np.arange(n)
        for _ in range(m):
            f_to_the_m = F[f_to_the_m]
        for c in c_for_inner.get(f_to_the_m.tobytes(), ()):
            if f[c] == c:
                data.append((f, c))
    return data


def reduce_by_conjugacy(data, auts):
    """Return one pair (alpha, c) from each orbit of the automorphism group of N, whose automorphisms
    are `auts`, acting on the pairs `data` by (alpha, c) -> (b alpha b^-1, b(c)).

    The pairs in an orbit define isomorphic groups, via x t^i -> b(x) t^i, so this loses nothing."""
    A = np.array(auts, dtype=np.int64)  # Row b is the automorphism b
    A_inv = np.argsort(A, axis=1)       # Row b is the inverse of b
    rows = np.arange(len(auts))[:, None]
    seen, representatives = set(), []
    for f, c in data:
        F = np.array(f, dtype=np.int64)
        if (F.tobytes(), c) in seen:
            continue
        representatives.append((f, c))
        conjugates = A[rows, F[A_inv]]  # Row b is b o alpha o b^-1
        seen.update(zip((row.tobytes() for row in conjugates), A[:, c].tolist()))
    return representatives


def build_cyclic_extension(N, m, alpha, c, name=None, description=None):
    """Return the group G = <N, t>, with t^m = c and t x t^-1 = alpha(x), where N is a Group, alpha is an
    automorphism of N (a tuple of indices), and c is the index of an element of N that satisfy the
    conditions of `valid_extension_data`.

    The elements of G are named 'x:i', for x in N.elements and 0 <= i < m (the delimiter is the one
    returned by N.direct_product_delimiter()), and 'x:i' stands for x t^i. They are in the same order
    as the elements of the direct product N * Z_m, and the product is
    (x t^i)(y t^j) = x alpha^i(y) t^(i+j), with t^m replaced by c."""
    T, _ = _int_view(N)
    n = len(T)
    indices = np.arange(n)
    F = np.array(alpha, dtype=np.int64)
    alpha_to_the = [indices]
    for _ in range(1, m):
        alpha_to_the.append(F[alpha_to_the[-1]])
    table = np.empty((n * m, n * m), dtype=np.int64)
    for i in range(m):
        X = T[indices[:, None], alpha_to_the[i][None, :]]  # x * alpha^i(y), for all x and y
        Xc = T[X, c]                                       # ... times c, for when i + j >= m
        for j in range(m):
            table[i::m, j::m] = (Xc if i + j >= m else X) * m + (i + j) % m
    delimiter = N.direct_product_delimiter()
    elements = [f"{x}{delimiter}{i}" for x in N.elements for i in range(m)]
    if name is None:
        name = f"{N.name}_ext{m}"
    if description is None:
        moved = {N.elements[a]: N.elements[b] for a, b in enumerate(alpha) if a != b}
        description = (f"Cyclic extension of {N.name} by Z{m}: t^{m} = {N.elements[c]}, "
                       f"t x t^-1 = alpha(x), alpha = {moved if moved else 'identity'}")
    G = make_finite_algebra(name, description, elements, table.tolist())
    if not isinstance(G, Group):
        raise ValueError("The extension data do not define a group.")
    return G


def cyclic_extensions(N, m, return_data=False):
    """Return a list of all of the groups G, up to isomorphism, that contain the Group N as a normal subgroup,
    with G/N cyclic of order m. There is one for each pair (alpha, c) of `valid_extension_data`, up to
    conjugation by automorphisms of N, and up to isomorphism of the resulting groups.

    If `return_data` is True, then a list of triples (G, alpha, c) is returned instead, where alpha is a tuple
    of indices and c is an index into N.elements. If c is the identity of N, then the extension splits (G
    is a semidirect product of N by Z_m); if c is not the identity, it may still split."""
    auts = automorphisms(N)
    representatives = reduce_by_conjugacy(valid_extension_data(N, m, auts), auts)
    found, buckets = [], defaultdict(list)
    for alpha, c in representatives:
        G = build_cyclic_extension(N, m, alpha, c, name=f"{N.name}_ext{m}_{len(found) + 1}")
        bucket = buckets[_invariant(*_int_view(G))]
        if any(G.fast_isomorphic(H) for H, _, _ in bucket):
            continue
        bucket.append((G, alpha, c))
        found.append((G, alpha, c))
    return found if return_data else [G for G, _, _ in found]


# ----------------------------------------------------------------------------------------
#  All of the solvable groups of a given order
# ----------------------------------------------------------------------------------------

def relabel_group(G, name=None):
    """Return a copy of the Group G whose elements are named '0', '1', ..., with the identity named
    '0', and, optionally, with a new name."""
    T, e = _int_view(G)
    n = len(T)
    order = [e] + [i for i in range(n) if i != e]
    new_index = np.empty(n, dtype=np.int64)  # old index -> new index
    new_index[order] = np.arange(n)
    table = new_index[T[np.ix_(order, order)]]
    return make_finite_algebra(name or G.name, G.description, [str(i) for i in range(n)], table.tolist())


_GROUPS_CACHE = {}


def groups_of_order(n, cache=None):
    """Return a list of all of the solvable groups of order n, up to isomorphism. Their elements are named
    '0', ..., 'n-1', with the identity named '0', and the groups are named 'G<n>_1', 'G<n>_2', ....

    The groups are built recursively, as the cyclic extensions by Z_p of the groups of order n/p, for each
    prime p that divides n, and the groups of each order are computed only once. Every solvable group has a
    normal subgroup of prime index, so this finds all of them. Every group of order less than 60 is
    solvable, so, up to order 59, this is all of the groups. A group that is not solvable (the smallest is
    A5, of order 60) is not found, so a warning is issued for n >= 60.

    The groups of each order that have been computed are remembered between calls, unless a dictionary is given
    as `cache`, and then that one is used and updated instead."""
    if isinstance(n, bool) or not isinstance(n, int) or n < 1:
        raise ValueError(f"n = {n!r} is invalid; it must be a positive integer.")
    if n >= 60:
        warnings.warn(f"Groups of order {n} that are not solvable are not found.", stacklevel=2)
    return list(_groups_of_order(n, _GROUPS_CACHE if cache is None else cache))


def _groups_of_order(n, cache):
    if n in cache:
        return cache[n]
    if n == 1:
        result = [make_finite_algebra("G1_1", "Trivial group", ["0"], [[0]])]
    else:
        result, buckets = [], defaultdict(list)
        for p in sorted(factorint(n)):
            for N in _groups_of_order(n // p, cache):
                for G in cyclic_extensions(N, p):
                    bucket = buckets[_invariant(*_int_view(G))]
                    if any(G.fast_isomorphic(H) for H in bucket):
                        continue
                    bucket.append(G)
                    result.append(G)
        result = [relabel_group(G, name=f"G{n}_{k}") for k, G in enumerate(result, start=1)]
    cache[n] = result
    return result
