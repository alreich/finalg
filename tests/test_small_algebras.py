"""
Unit tests for finalg.small_algebras: the catalogue of all groups of order 1 through 20.
"""

import itertools as it
from unittest import TestCase

from finalg import (small_groups, identify_small_group, generate_cyclic_group, generate_symmetric_group,
                    generate_dihedral_group, generate_dicyclic_group, generate_alternating_group,
                    generate_algebra_mod_n)
from finalg.group import Group

# Numbers of groups of order 1, 2, ..., 20, and numbers of abelian groups of those orders
GROUP_COUNTS = [1, 1, 1, 2, 1, 2, 1, 5, 2, 2, 1, 5, 1, 2, 1, 14, 1, 5, 1, 5]
ABELIAN_COUNTS = [1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 2, 1, 1, 1, 5, 1, 2, 1, 2]


class TestSmallGroups(TestCase):

    def test_keys_are_orders_1_through_20(self):
        self.assertEqual(sorted(small_groups()), list(range(1, 21)))

    def test_number_of_groups_of_each_order(self):
        catalogue = small_groups()
        self.assertEqual([len(catalogue[n]) for n in range(1, 21)], GROUP_COUNTS)
        self.assertEqual(sum(len(v) for v in catalogue.values()), 54)

    def test_every_algebra_is_a_group_of_the_right_order(self):
        for n, grps in small_groups().items():
            for g in grps:
                self.assertIsInstance(g, Group)
                self.assertEqual(g.order, n)

    def test_number_of_abelian_groups_of_each_order(self):
        catalogue = small_groups()
        self.assertEqual([sum(g.is_commutative() for g in catalogue[n]) for n in range(1, 21)],
                         ABELIAN_COUNTS)

    def test_no_two_groups_of_the_same_order_are_isomorphic(self):
        for n, grps in small_groups().items():
            for g1, g2 in it.combinations(grps, 2):
                self.assertFalse(g1.fast_isomorphic(g2), f"{g1.name} and {g2.name} are isomorphic")

    def test_names_within_an_order_are_unique(self):
        for n, grps in small_groups().items():
            names = [g.name for g in grps]
            self.assertEqual(len(names), len(set(names)), f"duplicate names for order {n}: {names}")

    def test_groups_of_order_4(self):
        z4, v4 = small_groups(4)
        self.assertEqual((z4.name, v4.name), ('Z4', 'V4'))
        self.assertFalse(z4.fast_isomorphic(v4))

    def test_order_argument_returns_the_list_for_that_order(self):
        self.assertEqual(small_groups(16), small_groups()[16])
        self.assertEqual(len(small_groups(16)), 14)
        self.assertEqual([g.name for g in small_groups(1)], ['Z1'])

    def test_order_16_contains_the_pauli_group_and_the_order_16_dihedral_and_quaternion_groups(self):
        names = [g.name for g in small_groups(16)]
        for name in ('Pauli', 'D8', 'Q16', 'SD16', 'M16', 'Z4_sdp_Z4', 'Z2Z2_sdp_Z4'):
            self.assertIn(name, names)

    def test_standard_groups_are_in_the_catalogue(self):
        catalogue = small_groups()
        standard = [(6, generate_symmetric_group(3)),
                    (8, generate_dihedral_group(4)[0]),
                    (8, generate_dicyclic_group(2)[0]),
                    (12, generate_alternating_group(4)),
                    (12, generate_dicyclic_group(3)[0]),
                    (16, generate_dihedral_group(8)[0]),
                    (20, generate_dihedral_group(10)[0]),
                    (20, generate_dicyclic_group(5)[0])]
        for n, g in standard:
            self.assertTrue(any(g.fast_isomorphic(s) for s in catalogue[n]), f"{g.name} is missing")

    def test_groups_are_created_once_and_shared(self):
        self.assertIs(small_groups(4)[0], small_groups()[4][0])
        self.assertIs(small_groups()[16][5], small_groups(16)[5])

    def test_changing_a_returned_container_does_not_change_the_catalogue(self):
        catalogue = small_groups()
        catalogue[4].clear()
        catalogue.pop(5)
        small_groups(6).append('junk')
        self.assertEqual(len(small_groups(4)), 2)
        self.assertIn(5, small_groups())
        self.assertEqual(len(small_groups(6)), 2)

    def test_invalid_order_raises(self):
        for bad in (0, -3, 21, 2.0, True, '4'):
            with self.assertRaises(ValueError, msg=repr(bad)):
                small_groups(bad)


class TestIdentifySmallGroup(TestCase):

    def test_every_small_group_identifies_as_itself(self):
        for grps in small_groups().values():
            for g in grps:
                self.assertIs(identify_small_group(g), g)

    def test_relabeled_copies_are_identified(self):
        for grps in small_groups().values():
            for g in grps:
                relabeled = g.reorder_elements(list(reversed(g.elements)))
                self.assertIs(identify_small_group(relabeled), g, g.name)

    def test_symmetric_group_on_3_letters_is_d3(self):
        self.assertEqual(identify_small_group(generate_symmetric_group(3)).name, 'D3')

    def test_cyclic_groups_are_identified(self):
        for n in (1, 2, 6, 12, 16, 18, 20):
            self.assertEqual(identify_small_group(generate_cyclic_group(n)).name, f"Z{n}")

    def test_direct_product_of_coprime_cyclic_groups_is_cyclic(self):
        z2_x_z3 = generate_cyclic_group(2) * generate_cyclic_group(3)
        self.assertEqual(identify_small_group(z2_x_z3).name, 'Z6')

    def test_klein_four_group(self):
        v4 = generate_cyclic_group(2) * generate_cyclic_group(2)
        self.assertEqual(identify_small_group(v4).name, 'V4')

    def test_dihedral_group_of_order_20_is_not_the_dicyclic_group_of_order_20(self):
        d10 = identify_small_group(generate_dihedral_group(10)[0])
        dic5 = identify_small_group(generate_dicyclic_group(5)[0])
        self.assertEqual((d10.name, dic5.name), ('D10', 'Dic5'))

    def test_groups_of_order_greater_than_20_return_none(self):
        self.assertIsNone(identify_small_group(generate_cyclic_group(21)))
        self.assertIsNone(identify_small_group(generate_alternating_group(5)))

    def test_algebras_that_are_not_groups_return_none(self):
        self.assertIsNone(identify_small_group(generate_algebra_mod_n(4)))
        self.assertIsNone(identify_small_group('not an algebra'))
