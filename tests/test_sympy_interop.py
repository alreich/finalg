"""
Unit tests for finalg.sympy_interop: converting between a
sympy.combinatorics.PermutationGroup and a finalg Group, in both directions.
"""

from unittest import TestCase

from sympy.combinatorics import (
    Permutation,
    PermutationGroup,
    SymmetricGroup,
    AlternatingGroup,
    CyclicGroup,
    DihedralGroup,
)

from finalg.sympy_interop import from_sympy_permutation_group, to_sympy_permutation_group
from finalg.group import Group
from finalg import (
    generate_symmetric_group,
    generate_dihedral_group,
    generate_cyclic_group,
    make_finite_algebra,
    examples,
)


class TestFromSympyPermutationGroup(TestCase):

    def test_rejects_non_permutation_group(self):
        with self.assertRaises(TypeError):
            from_sympy_permutation_group("not a group")

    def test_returns_a_finalg_group(self):
        alg, _ = from_sympy_permutation_group(SymmetricGroup(3))
        self.assertIsInstance(alg, Group)

    def test_order_matches_sympy(self):
        for sympy_grp in [SymmetricGroup(4), AlternatingGroup(4), CyclicGroup(6), DihedralGroup(5)]:
            with self.subTest(sympy_grp=sympy_grp):
                alg, _ = from_sympy_permutation_group(sympy_grp)
                self.assertEqual(alg.order, sympy_grp.order())

    def test_custom_name_and_description(self):
        alg, _ = from_sympy_permutation_group(CyclicGroup(4), name='MyZ4', description='custom desc')
        self.assertEqual(alg.name, 'MyZ4')
        self.assertEqual(alg.description, 'custom desc')

    def test_elem_dict_maps_names_to_sympy_permutations(self):
        alg, elem_dict = from_sympy_permutation_group(SymmetricGroup(3))
        self.assertEqual(set(elem_dict.keys()), set(alg.elements))
        for perm in elem_dict.values():
            self.assertIsInstance(perm, Permutation)

    def test_trivial_group(self):
        triv = PermutationGroup([Permutation([0, 1, 2])])
        alg, _ = from_sympy_permutation_group(triv)
        self.assertEqual(alg.order, 1)

    def test_isomorphic_to_finalg_symmetric_group(self):
        alg, _ = from_sympy_permutation_group(SymmetricGroup(3))
        s3 = generate_symmetric_group(3)
        self.assertTrue(alg.fast_isomorphic(s3))

    def test_isomorphic_to_finalg_dihedral_group(self):
        alg, _ = from_sympy_permutation_group(DihedralGroup(4))
        d4, _, _ = generate_dihedral_group(4)
        self.assertTrue(alg.fast_isomorphic(d4))

    def test_isomorphic_to_finalg_cyclic_group(self):
        alg, _ = from_sympy_permutation_group(CyclicGroup(6))
        z6 = generate_cyclic_group(6)
        self.assertTrue(alg.fast_isomorphic(z6))

    def test_custom_subgroup_generates_full_symmetric_group(self):
        # A 5-cycle and a transposition together generate all of S5.
        p1 = Permutation([1, 2, 3, 4, 0])
        p2 = Permutation([1, 0, 2, 3, 4])
        g = PermutationGroup([p1, p2])
        alg, _ = from_sympy_permutation_group(g)
        self.assertEqual(alg.order, 120)
        self.assertFalse(alg.is_commutative())

    def test_derived_subgroup_of_s4_is_a4(self):
        s4 = SymmetricGroup(4)
        der = s4.derived_subgroup()
        alg, _ = from_sympy_permutation_group(der)
        self.assertEqual(alg.order, 12)

    def test_sylow_subgroup_of_s4(self):
        s4 = SymmetricGroup(4)
        syl = s4.sylow_subgroup(2)
        alg, _ = from_sympy_permutation_group(syl)
        self.assertEqual(alg.order, 8)

    def test_element_names_are_readable_cycle_notation(self):
        alg, _ = from_sympy_permutation_group(CyclicGroup(3))
        # Every element name should look like finalg's own cycle notation,
        # e.g. '(0)(1)(2)' or '(0 1 2)', not a raw SymPy repr.
        for elem in alg.elements:
            self.assertNotIn('Permutation', elem)
            self.assertTrue(elem.startswith('('))

    def test_verify_true_passes_for_a_genuine_group(self):
        # Should not raise: a real PermutationGroup is always associative
        # and fully invertible, so verify=True's extra checks should pass.
        alg, _ = from_sympy_permutation_group(SymmetricGroup(4), verify=True)
        self.assertTrue(alg.is_associative())
        self.assertTrue(alg.has_inverses())

    def test_table_is_associative_even_without_verify(self):
        # verify=False (the default) skips the eager check, but the resulting
        # table should still genuinely be associative -- confirm by asking
        # for it explicitly after the fact.
        alg, _ = from_sympy_permutation_group(DihedralGroup(4))
        self.assertTrue(alg.is_associative())

    def test_moderately_large_group_converts(self):
        # S6 has 720 elements; this is mainly a regression guard against
        # reintroducing an O(n^3) construction path.
        alg, _ = from_sympy_permutation_group(SymmetricGroup(6))
        self.assertEqual(alg.order, 720)
        self.assertTrue(alg.has_identity())


def _exact_homomorphism(finalg_group, perm_of):
    """True iff perm_of[a] * perm_of[b] (SymPy composition) == perm_of[a*b]
    (finalg's own op) for every pair -- i.e. a genuine, order-preserving
    embedding, not merely an isomorphism up to relabeling."""
    elements = list(finalg_group.elements)
    return all(
        perm_of[a] * perm_of[b] == perm_of[finalg_group.op(a, b)]
        for a in elements
        for b in elements
    )


class TestToSympyPermutationGroup(TestCase):

    def test_rejects_non_group(self):
        rps = make_finite_algebra('RPS', 'Rock Paper Scissors',
                                   ['r', 'p', 's'], [[0, 1, 0], [1, 1, 2], [0, 2, 2]])
        with self.assertRaises(TypeError):
            to_sympy_permutation_group(rps)

    def test_rejects_quasigroup(self):
        latin = examples[17]
        with self.assertRaises(TypeError):
            to_sympy_permutation_group(latin)

    def test_rejects_ring_and_field(self):
        f4 = examples[9]  # Field with 4 elements -- a subclass of Group via +
        with self.assertRaises(TypeError):
            to_sympy_permutation_group(f4)

    def test_returns_a_sympy_permutation_group(self):
        z4 = generate_cyclic_group(4)
        sympy_group, _ = to_sympy_permutation_group(z4)
        self.assertIsInstance(sympy_group, PermutationGroup)

    def test_order_matches(self):
        for grp in [generate_symmetric_group(3), generate_dihedral_group(4)[0],
                    generate_cyclic_group(6), examples[13]]:  # Q8
            with self.subTest(grp=grp.name):
                sympy_group, _ = to_sympy_permutation_group(grp)
                self.assertEqual(sympy_group.order(), grp.order)

    def test_elem_dict_covers_every_element_not_just_generators(self):
        s3 = generate_symmetric_group(3)
        _, perm_of = to_sympy_permutation_group(s3)
        self.assertEqual(set(perm_of.keys()), set(s3.elements))
        for perm in perm_of.values():
            self.assertIsInstance(perm, Permutation)

    def test_exact_homomorphism_abelian(self):
        z6 = generate_cyclic_group(6)
        _, perm_of = to_sympy_permutation_group(z6)
        self.assertTrue(_exact_homomorphism(z6, perm_of))

    def test_exact_homomorphism_nonabelian(self):
        # This is the case that actually distinguishes left- from
        # right-regular representation under SymPy's composition convention.
        s3 = generate_symmetric_group(3)
        _, perm_of = to_sympy_permutation_group(s3)
        self.assertTrue(_exact_homomorphism(s3, perm_of))

    def test_exact_homomorphism_quaternion(self):
        q8 = examples[13]
        _, perm_of = to_sympy_permutation_group(q8)
        self.assertTrue(_exact_homomorphism(q8, perm_of))

    def test_round_trip_is_isomorphic(self):
        s4 = generate_symmetric_group(4)
        sympy_group, _ = to_sympy_permutation_group(s4)
        back, _ = from_sympy_permutation_group(sympy_group)
        self.assertTrue(back.fast_isomorphic(s4))

    def test_order_60_group_converts_quickly(self):
        # Mainly a regression guard against reintroducing a dependency on
        # finding a *minimal* generating set, which is its own (slow) search.
        a5 = examples[15]
        sympy_group, _ = to_sympy_permutation_group(a5)
        self.assertEqual(sympy_group.order(), 60)

