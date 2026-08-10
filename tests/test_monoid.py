"""
Unit tests for finalg.monoid.Monoid.
"""

import numpy as np
from unittest import TestCase

from finalg import generate_commutative_monoid, generate_cyclic_group
from finalg.monoid import Monoid


class TestMonoidConstruction(TestCase):

    def test_generate_commutative_monoid_is_monoid(self):
        m = generate_commutative_monoid(4)
        self.assertEqual(type(m).__name__, 'Monoid')
        self.assertEqual(m.elements, ('a0', 'a1', 'a2', 'a3'))
        self.assertEqual(m.identity, 'a1')

    def test_direct_construction_raises_without_identity(self):
        ex141_tbl = [[0, 3, 0, 3, 0, 3],
                     [1, 4, 1, 4, 1, 4],
                     [2, 5, 2, 5, 2, 5],
                     [3, 0, 3, 0, 3, 0],
                     [4, 1, 4, 1, 4, 1],
                     [5, 2, 5, 2, 5, 2]]
        with self.assertRaises(ValueError):
            Monoid('bad', 'no identity', ['a', 'b', 'c', 'd', 'e', 'f'], ex141_tbl)

    def test_direct_construction_raises_when_not_associative(self):
        with self.assertRaises(ValueError):
            Monoid('bad', 'not associative', ['r', 'p', 's'],
                   [[0, 1, 0], [1, 1, 2], [0, 2, 2]])


class TestElementOrder(TestCase):

    def setUp(self):
        self.m4 = generate_commutative_monoid(4)  # mult mod 4
        self.z4 = generate_cyclic_group(4)

    def test_identity_has_order_1(self):
        self.assertEqual(self.m4.element_order('a1'), 1)

    def test_other_orders_in_z4(self):
        self.assertEqual(self.z4.element_order('0'), 1)
        self.assertEqual(self.z4.element_order('1'), 4)
        self.assertEqual(self.z4.element_order('2'), 2)
        self.assertEqual(self.z4.element_order('3'), 4)

    def test_element_order_is_cached(self):
        # Call twice; the cached value should be returned and remain correct.
        first = self.z4.element_order('2')
        second = self.z4.element_order('2')
        self.assertEqual(first, second)
        self.assertEqual(second, 2)


class TestUnits(TestCase):

    def setUp(self):
        self.m4 = generate_commutative_monoid(4)  # mult mod 4: 0,1,2,3

    def test_units_returns_names_by_default(self):
        self.assertEqual(self.m4.units(), ['a1', 'a3'])

    def test_units_returns_indices_when_requested(self):
        self.assertEqual(self.m4.units(return_names=False), [1, 3])

    def test_units_subgroup(self):
        sub = self.m4.units_subgroup()
        self.assertEqual(sub.elements, ('a1', 'a3'))
        self.assertEqual(sub.op('a3', 'a3'), 'a1')  # 3*3 mod 4 = 1


class TestRegularRepresentation(TestCase):

    def setUp(self):
        self.z4 = generate_cyclic_group(4)

    def test_dense_regular_representation_shapes(self):
        mapping, inv_mapping, elem_to_arr, arr_to_elem = self.z4.regular_representation()
        self.assertEqual(set(mapping.keys()), set(self.z4.elements))
        for arr in mapping.values():
            self.assertEqual(arr.shape, (4, 4))

    def test_regular_representation_identity_matrix_maps_to_identity_elem(self):
        mapping, inv_mapping, elem_to_arr, arr_to_elem = self.z4.regular_representation()
        self.assertTrue(np.array_equal(elem_to_arr(self.z4.identity), np.eye(4, dtype=int)))

    def test_verify_regular_representation_true(self):
        _, _, elem_to_arr, arr_to_elem = self.z4.regular_representation()
        self.assertTrue(self.z4.verify_regular_representation(elem_to_arr, arr_to_elem))

    def test_sparse_regular_representation(self):
        import scipy.sparse as sp
        mapping, _, _, _ = self.z4.regular_representation(sparse='CSR')
        for arr in mapping.values():
            self.assertIsInstance(arr, sp.csr_array)

    def test_array_to_element_round_trip(self):
        _, _, elem_to_arr, arr_to_elem = self.z4.regular_representation()
        for elem in self.z4.elements:
            arr = elem_to_arr(elem)
            self.assertEqual(arr_to_elem(arr), elem)

    def test_matrix_product_matches_group_operation(self):
        _, _, elem_to_arr, arr_to_elem = self.z4.regular_representation()
        a, b = '1', '2'
        prod_matrix = np.dot(elem_to_arr(a), elem_to_arr(b))
        self.assertEqual(arr_to_elem(prod_matrix), self.z4.op(a, b))


class TestMonoidIsomorphismMappings(TestCase):

    def test_make_element_mappings_fixes_identity(self):
        m1 = generate_commutative_monoid(4)
        m2 = generate_commutative_monoid(4)
        mappings = m1.make_element_mappings(m2)
        # Every mapping must send the identity to the identity.
        self.assertTrue(all(mapping[m1.identity] == m2.identity for mapping in mappings))
        # The identity is now genuinely excluded before permuting, so there are
        # (n-1)! mappings over the remaining elements, with no duplicates.
        self.assertEqual(len(mappings), 6)  # 3! for the other 3 elements
        self.assertEqual(len(mappings), len({tuple(sorted(m.items())) for m in mappings}))

    def test_make_element_mappings_different_orders_raises(self):
        m1 = generate_commutative_monoid(4)
        m2 = generate_commutative_monoid(5)
        with self.assertRaises(ValueError):
            m1.make_element_mappings(m2)
