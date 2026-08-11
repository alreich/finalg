"""
Unit tests for finalg.sympy_interop: converting a
sympy.combinatorics.PermutationGroup into a finalg Group.
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

from finalg.sympy_interop import from_sympy_permutation_group
from finalg.group import Group
from finalg import generate_symmetric_group, generate_dihedral_group, generate_cyclic_group


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
