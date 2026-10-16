"""
Unit tests for finalg.small_algebras: the catalogue of all groups of order 1 through 40, which is
loaded from the JSON files in the package finalg.data.small_groups. (The loading functions and the
data files themselves are tested in test_small_groups_loader.py.)
"""

import itertools as it
from unittest import TestCase

import finalg.small_algebras as sa
from finalg import (small_groups, identify_small_group, generate_cyclic_group, generate_symmetric_group,
                    generate_dihedral_group, generate_dicyclic_group, generate_alternating_group,
                    generate_algebra_mod_n)
from finalg.group import Group
from finalg.utilities import action_from_generators

# Numbers of groups of order 1, 2, ..., 40, and numbers of abelian groups of those orders
GROUP_COUNTS = [1, 1, 1, 2, 1, 2, 1, 5, 2, 2, 1, 5, 1, 2, 1, 14, 1, 5, 1, 5,
                2, 2, 1, 15, 2, 2, 5, 4, 1, 4, 1, 51, 1, 2, 1, 14, 1, 2, 2, 14]
ABELIAN_COUNTS = [1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 2, 1, 1, 1, 5, 1, 2, 1, 2,
                  1, 1, 1, 3, 2, 1, 3, 2, 1, 1, 1, 7, 1, 1, 1, 4, 1, 1, 1, 3]


class TestSmallGroups(TestCase):

    def test_max_small_order_is_40(self):
        self.assertEqual(sa.MAX_SMALL_ORDER, 40)

    def test_keys_are_orders_1_through_40(self):
        self.assertEqual(sorted(small_groups()), list(range(1, 41)))

    def test_number_of_groups_of_each_order(self):
        catalogue = small_groups()
        self.assertEqual([len(catalogue[n]) for n in range(1, 41)], GROUP_COUNTS)
        self.assertEqual(sum(len(v) for v in catalogue.values()), 181)

    def test_every_algebra_is_a_group_of_the_right_order(self):
        for n, grps in small_groups().items():
            for g in grps:
                self.assertIsInstance(g, Group)
                self.assertEqual(g.order, n)

    def test_number_of_abelian_groups_of_each_order(self):
        catalogue = small_groups()
        self.assertEqual([sum(g.is_commutative() for g in catalogue[n]) for n in range(1, 41)],
                         ABELIAN_COUNTS)

    def test_no_two_groups_of_the_same_order_are_isomorphic(self):
        for n, grps in small_groups().items():
            for g1, g2 in it.combinations(grps, 2):
                self.assertFalse(g1.fast_isomorphic(g2), f"{g1.name} and {g2.name} are isomorphic")

    def test_names_within_an_order_are_unique(self):
        for n, grps in small_groups().items():
            names = [g.name for g in grps]
            self.assertEqual(len(names), len(set(names)), f"duplicate names for order {n}: {names}")

    def test_the_groups_are_in_the_order_of_the_smallgroup_ids(self):
        index = sa.groups_index()
        flat = [g for n in range(1, 41) for g in small_groups(n)]
        self.assertEqual([(g.order, g.name) for g in flat], [(r['order'], r['name']) for r in index])

    def test_groups_of_order_4(self):
        z4, v4 = small_groups(4)
        self.assertEqual((z4.name, v4.name), ('Z4', 'V4'))
        self.assertFalse(z4.fast_isomorphic(v4))

    def test_order_argument_returns_the_list_for_that_order(self):
        self.assertEqual(small_groups(16), small_groups()[16])
        self.assertEqual(len(small_groups(16)), 14)
        self.assertEqual(len(small_groups(32)), 51)
        self.assertEqual([g.name for g in small_groups(1)], ['Z1'])

    def test_order_16_contains_the_pauli_group_and_the_order_16_dihedral_and_quaternion_groups(self):
        names = [g.name for g in small_groups(16)]
        for name in ('Pauli', 'D8', 'Q16', 'SD16', 'M16', 'Z4_sdp_Z4', 'Z2Z2_sdp_Z4'):
            self.assertIn(name, names)

    def test_the_groups_that_are_new_with_order_21_to_40_are_there(self):
        names = {g.name for g in small_groups(24)}
        for name in ('S4', 'SL2_3', 'Dic6'):
            self.assertIn(name, names)
        self.assertEqual(small_groups(21)[0].name, 'F21')
        self.assertIn('Q32', [g.name for g in small_groups(32)])
        self.assertEqual([g.name for g in small_groups(32) if g.name.startswith('SG32_')],
                         ['SG32_8', 'SG32_15', 'SG32_32'])

    def test_standard_groups_are_in_the_catalogue(self):
        standard = [generate_symmetric_group(3), generate_dihedral_group(4)[0], generate_dicyclic_group(2)[0],
                    generate_alternating_group(4), generate_dicyclic_group(3)[0], generate_dihedral_group(8)[0],
                    generate_dihedral_group(10)[0], generate_dicyclic_group(5)[0],
                    generate_symmetric_group(4), generate_dihedral_group(20)[0], generate_dicyclic_group(8)[0],
                    generate_dicyclic_group(10)[0]]
        for g in standard:
            small = identify_small_group(g)
            self.assertIsNotNone(small, g.name)
            self.assertEqual(small.order, g.order)

    def test_groups_are_created_once_and_shared(self):
        self.assertIs(small_groups(4)[0], small_groups()[4][0])
        self.assertIs(small_groups()[16][5], small_groups(16)[5])
        self.assertIs(small_groups(32)[19], small_groups()[32][19])

    def test_changing_a_returned_container_does_not_change_the_catalogue(self):
        catalogue = small_groups()
        catalogue[4].clear()
        catalogue.pop(5)
        small_groups(6).append('junk')
        self.assertEqual(len(small_groups(4)), 2)
        self.assertIn(5, small_groups())
        self.assertEqual(len(small_groups(6)), 2)

    def test_invalid_order_raises(self):
        for bad in (0, -3, 41, 2.0, True, '4'):
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

    def test_groups_of_order_32_and_40_are_identified_quickly_enough_to_be_useful(self):
        # Before the generating-set search was fixed, this took about 40 seconds for some of these groups
        for n in (32, 40):
            for g in small_groups(n):
                relabeled = g.reorder_elements([g.elements[0]] + list(reversed(g.elements[1:])))
                self.assertIs(identify_small_group(relabeled), g, g.name)

    def test_symmetric_group_on_3_letters_is_d3(self):
        self.assertEqual(identify_small_group(generate_symmetric_group(3)).name, 'D3')

    def test_symmetric_group_on_4_letters_is_s4(self):
        self.assertEqual(identify_small_group(generate_symmetric_group(4)).name, 'S4')

    def test_cyclic_groups_are_identified(self):
        for n in (1, 2, 6, 12, 16, 18, 20, 21, 24, 30, 32, 36, 40):
            self.assertEqual(identify_small_group(generate_cyclic_group(n)).name, f"Z{n}")

    def test_direct_product_of_coprime_cyclic_groups_is_cyclic(self):
        z2_x_z3 = generate_cyclic_group(2) * generate_cyclic_group(3)
        self.assertEqual(identify_small_group(z2_x_z3).name, 'Z6')
        z5_x_z8 = generate_cyclic_group(5) * generate_cyclic_group(8)
        self.assertEqual(identify_small_group(z5_x_z8).name, 'Z40')

    def test_klein_four_group(self):
        v4 = generate_cyclic_group(2) * generate_cyclic_group(2)
        self.assertEqual(identify_small_group(v4).name, 'V4')

    def test_dihedral_group_of_order_20_is_not_the_dicyclic_group_of_order_20(self):
        d10 = identify_small_group(generate_dihedral_group(10)[0])
        dic5 = identify_small_group(generate_dicyclic_group(5)[0])
        self.assertEqual((d10.name, dic5.name), ('D10', 'Dic5'))

    def test_the_generalized_quaternion_group_of_order_32(self):
        self.assertEqual(identify_small_group(generate_dicyclic_group(8)[0]).name, 'Q32')

    def test_a_group_built_as_a_semidirect_product_is_identified(self):
        q8, z3 = small_groups(8)[3], generate_cyclic_group(3)
        action = action_from_generators(q8, z3, {'1': {'a': 'ax', 'x': 'a'}})
        self.assertEqual(identify_small_group(q8.semidirect_product(z3, action)).name, 'SL2_3')

    def test_groups_of_order_greater_than_40_return_none(self):
        self.assertIsNone(identify_small_group(generate_cyclic_group(41)))
        self.assertIsNone(identify_small_group(generate_alternating_group(5)))

    def test_algebras_that_are_not_groups_return_none(self):
        self.assertIsNone(identify_small_group(generate_algebra_mod_n(4)))
        self.assertIsNone(identify_small_group('not an algebra'))
