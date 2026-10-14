# ===============================================
#   Small algebras: a catalogue of small groups
# ===============================================
"""A catalogue of all groups of order 1 through 20, one for each isomorphism class (54 groups in all).

The catalogue is built, the first time it is needed, from simpler groups: cyclic groups, direct
products (``*`` and ``**``), the dihedral, dicyclic, and alternating group generators,
``Group.semidirect_product``, and ``Group.central_product``. It is then cached, so that every
call returns the same Group objects.

Two functions are provided:

* ``small_groups(order=None)`` returns the dictionary of all small groups, keyed by order,
  or, if an order is given, the list of all groups of that order.
* ``identify_small_group(alg)`` returns the small group that a given group is isomorphic to,
  or None.

Group names follow the conventions of the generators: ``Z6`` is cyclic of order 6, ``D4`` is the
dihedral group of order 8 (the symmetries of a square), ``Q8`` and ``Q16`` are the quaternion and
generalized quaternion groups, ``Dic3`` and ``Dic5`` are the dicyclic groups of orders 12 and 20,
``Z4_x_Z2`` is a direct product, and so on.
"""

from functools import lru_cache

from finalg.group import Group

MAX_SMALL_ORDER = 20

# The known number of groups of order 1, 2, ..., 20 (OEIS A000001), used as a sanity check.
_GROUP_COUNTS = (1, 1, 1, 2, 1, 2, 1, 5, 2, 2, 1, 5, 1, 2, 1, 14, 1, 5, 1, 5)


@lru_cache(maxsize=None)
def _catalogue():
    """Build, once, the dictionary that maps each order from 1 to 20 to the list of all groups of
    that order, up to isomorphism."""
    from finalg.algebra_generators import (generate_cyclic_group, generate_dihedral_group,
                                           generate_dicyclic_group, generate_alternating_group)

    def z(n):
        """The cyclic group of order n."""
        return generate_cyclic_group(n)

    def dihedral(n):
        """The dihedral group D_n, of order 2n."""
        return generate_dihedral_group(n)[0]

    def dicyclic(n):
        """The dicyclic group Dic_n, of order 4n."""
        return generate_dicyclic_group(n)[0]

    def mult_action(modulus, k):
        """Action of Z_m on Z_modulus (elements are strings of integers), where h acts by x -> k^h * x."""
        return lambda h, x: str((int(x) * pow(k, int(h), modulus)) % modulus)

    # Groups of order 4, 6, 8, 12, 16, 18, 20 that are named or need more than direct products.
    v4 = z(2) * z(2)
    v4.name = 'V4'
    v4.description += ' (Klein-4 group)'
    d3, d4, d5, d6, d7, d8, d9, d10 = [dihedral(n) for n in range(3, 11)]
    a4 = generate_alternating_group(4)
    q8, dic3, q16, dic5 = dicyclic(2), dicyclic(3), dicyclic(4), dicyclic(5)

    # Order 16
    z4_sdp_z4 = z(4).semidirect_product(z(4), mult_action(4, 3), name='Z4_sdp_Z4')   # inversion
    m16 = z(8).semidirect_product(z(2), mult_action(8, 5), name='M16')               # x -> x^5
    sd16 = z(8).semidirect_product(z(2), mult_action(8, 3), name='SD16')             # x -> x^3
    swap = lambda h, x: x if int(h) % 2 == 0 else ':'.join(reversed(x.split(':')))
    z2z2_sdp_z4 = (z(2) * z(2)).semidirect_product(z(4), swap, name='Z2Z2_sdp_Z4')   # swap the two Z2's
    d4_center = [c for c in d4.center() if c != d4.identity]                         # one involution
    pauli = z(4).central_product(d4, ('2', d4_center[0]), name='Pauli')              # Z4 o D4

    # Order 18
    invert = lambda h, x: x if int(h) % 2 == 0 else ':'.join(str((-int(t)) % 3) for t in x.split(':'))
    dih_z3z3 = (z(3) * z(3)).semidirect_product(z(2), invert, name='Dih_Z3xZ3')      # (Z3 x Z3) X| Z2

    # Order 20
    f20 = z(5).semidirect_product(z(4), mult_action(5, 2), name='F20')               # Z5 X| Z4, faithful

    groups = {
        1: [z(1)],
        2: [z(2)],
        3: [z(3)],
        4: [z(4), v4],
        5: [z(5)],
        6: [z(6), d3],
        7: [z(7)],
        8: [z(8), z(4) * z(2), z(2) ** 3, d4, q8],
        9: [z(9), z(3) * z(3)],
        10: [z(10), d5],
        11: [z(11)],
        12: [z(12), z(6) * z(2), a4, d6, dic3],
        13: [z(13)],
        14: [z(14), d7],
        15: [z(15)],
        16: [z(16), z(4) * z(4), z2z2_sdp_z4, z4_sdp_z4, z(8) * z(2), m16, d8, sd16, q16,
             z(4) * z(2) ** 2, d4 * z(2), q8 * z(2), pauli, z(2) ** 4],
        17: [z(17)],
        18: [z(18), d9, d3 * z(3), dih_z3z3, z(6) * z(3)],
        19: [z(19)],
        20: [z(20), z(10) * z(2), d10, dic5, f20],
    }

    for n, grps in groups.items():
        if len(grps) != _GROUP_COUNTS[n - 1] or any(g.order != n for g in grps):
            raise RuntimeError(f"The catalogue of small groups of order {n} is wrong: "
                               f"{len(grps)} groups, but there should be {_GROUP_COUNTS[n - 1]}.")
    return groups


def small_groups(order=None):
    """Return the small groups, which are all the groups of order 1 through 20, up to isomorphism.

    If `order` is None (the default), a dictionary is returned whose keys are the orders 1, ..., 20 and
    whose values are the lists of all groups of that order, e.g., small_groups()[4] is [Z4, V4].
    If `order` is an integer from 1 to 20, then just the list of groups of that order is returned,
    e.g., small_groups(4) is [Z4, V4].

    The Group objects are created once and then shared, so a group returned by one call is the very
    same object returned by another. The dictionary and lists are new each time, so changing them does
    not change the catalogue."""
    catalogue = _catalogue()
    if order is None:
        return {n: list(grps) for n, grps in catalogue.items()}
    if isinstance(order, bool) or not isinstance(order, int) or not 1 <= order <= MAX_SMALL_ORDER:
        raise ValueError(f"order = {order!r} is invalid; it must be an integer from 1 to {MAX_SMALL_ORDER}, "
                         f"or None for all of the small groups.")
    return list(catalogue[order])


def identify_small_group(alg):
    """Return the small group that the given algebra is isomorphic to, or None if there isn't one.

    There is one only if `alg` is a Group of order 20 or less, and then there is exactly one, because
    the small groups contain every isomorphism class once. Isomorphism is tested with `fast_isomorphic`."""
    if not isinstance(alg, Group) or alg.order > MAX_SMALL_ORDER:
        return None
    for small in _catalogue()[alg.order]:
        if alg.fast_isomorphic(small):
            return small
    return None
