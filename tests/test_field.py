"""
Unit tests for finalg.field.Field and is_field.
"""

from unittest import TestCase

from finalg import generate_algebra_mod_n, make_finite_algebra
from finalg.field import Field, is_field
from finalg.group import Group


class TestIsField(TestCase):

    def test_prime_modulus_yields_field(self):
        # Multiplication table for Z/5Z (excluding zero-row/col considerations
        # handled internally by is_field).
        elements = ('0', '1', '2', '3', '4')
        mult_tbl = [[(a * b) % 5 for b in range(5)] for a in range(5)]
        result = is_field('0', elements, mult_tbl)
        self.assertIsInstance(result, Group)
        self.assertTrue(result.is_commutative())

    def test_composite_modulus_is_not_a_field(self):
        elements = ('0', '1', '2', '3')
        mult_tbl = [[(a * b) % 4 for b in range(4)] for a in range(4)]
        result = is_field('0', elements, mult_tbl)
        self.assertFalse(result)

    def test_trivial_single_element_is_not_a_field(self):
        self.assertFalse(is_field('0', ('0',), [[0]]))


class TestFieldConstruction(TestCase):

    def test_generate_algebra_mod_n_prime_is_field(self):
        f5 = generate_algebra_mod_n(5)
        self.assertEqual(type(f5).__name__, 'Field')

    def test_direct_construction_raises_for_non_field_ring(self):
        add_tbl = [[0, 1, 2, 3], [1, 2, 3, 0], [2, 3, 0, 1], [3, 0, 1, 2]]
        mult_tbl = [[0, 0, 0, 0], [0, 1, 2, 3], [0, 2, 0, 2], [0, 3, 2, 1]]
        with self.assertRaises(ValueError):
            Field('bad', 'not a field', ['0', '1', '2', '3'], add_tbl, mult_tbl)


class TestFieldOperations(TestCase):

    def setUp(self):
        self.f5 = generate_algebra_mod_n(5)

    def test_mult_abelian_subgroup(self):
        sub = self.f5.mult_abelian_subgroup()
        self.assertEqual(set(sub.elements), {'1', '2', '3', '4'})
        self.assertTrue(sub.is_commutative())

    def test_mult_inv(self):
        self.assertEqual(self.f5.mult_inv('2'), '3')  # 2*3 = 6 = 1 mod 5
        self.assertEqual(self.f5.mult(self.f5.mult_inv('2'), '2'), self.f5.one)

    def test_mult_inv_of_zero_is_none(self):
        self.assertIsNone(self.f5.mult_inv('0'))

    def test_div(self):
        self.assertEqual(self.f5.div('2', '3'), self.f5.mult('2', self.f5.mult_inv('3')))

    def test_div_by_zero_is_none(self):
        self.assertIsNone(self.f5.div('2', '0'))

    def test_element_to_power(self):
        # 2^3 mod 5 = 3
        self.assertEqual(self.f5.element_to_power('2', 3), '3')

    def test_field_has_no_zero_divisors(self):
        self.assertEqual(self.f5.zero_divisors(), [])

    def test_field_elements_all_nonzero_are_units(self):
        self.assertEqual(sorted(self.f5.units()), ['1', '2', '3', '4'])


class TestFieldWithSymbolicElements(TestCase):
    """F4, built with symbolic (non-numeric) element names '0','1','a','1+a',
    exercises Field over element names that aren't simple integers-as-strings."""

    def setUp(self):
        self.f4_elems = ['0', '1', 'a', '1+a']
        self.f4_add_table = [['0', '1', 'a', '1+a'],
                             ['1', '0', '1+a', 'a'],
                             ['a', '1+a', '0', '1'],
                             ['1+a', 'a', '1', '0']]
        self.f4_mult_table = [['0', '0', '0', '0'],
                              ['0', '1', 'a', '1+a'],
                              ['0', 'a', '1+a', '1'],
                              ['0', '1+a', '1', 'a']]
        self.f4 = make_finite_algebra('F4', 'Field with 4 elements',
                                      self.f4_elems, self.f4_add_table, self.f4_mult_table)

    def test_is_a_field(self):
        self.assertEqual(type(self.f4).__name__, 'Field')

    def test_equality_with_independently_built_copy(self):
        f4x = make_finite_algebra('F4X', 'Field with 4 elements EXTRA COPY',
                                  self.f4_elems, self.f4_add_table, self.f4_mult_table)
        self.assertEqual(self.f4, f4x)

    def test_mult_inv(self):
        self.assertEqual(self.f4.mult_inv('a'), '1+a')

    def test_div(self):
        self.assertEqual(self.f4.div('1', '1+a'), 'a')

    def test_mult_abelian_subgroup(self):
        subgrp = self.f4.mult_abelian_subgroup()
        self.assertEqual(subgrp.elements, ('1', 'a', '1+a'))
        self.assertEqual(subgrp.table.tolist(), [[0, 1, 2], [1, 2, 0], [2, 0, 1]])
