"""
Unit tests for finalg.finite_operator.FiniteOperator.
"""

from unittest import TestCase

from finalg.cayley_table import CayleyTable
from finalg.finite_operator import FiniteOperator


class TestFiniteOperator(TestCase):

    def setUp(self):
        # Z4 cyclic group
        self.elements = ('e', 'a', 'a2', 'a3')
        self.table = CayleyTable([[0, 1, 2, 3],
                                   [1, 2, 3, 0],
                                   [2, 3, 0, 1],
                                   [3, 0, 1, 2]])
        self.op = FiniteOperator(self.elements, 'e', self.table)

        # RPS Magma (non-associative), no identity
        self.rps_elements = ('r', 'p', 's')
        self.rps_table = CayleyTable([[0, 1, 0], [1, 1, 2], [0, 2, 2]])
        self.rps_op = FiniteOperator(self.rps_elements, None, self.rps_table)

    def test_zero_args_returns_identity(self):
        self.assertEqual(self.op(), 'e')

    def test_zero_args_returns_none_when_no_identity(self):
        self.assertIsNone(self.rps_op())

    def test_one_arg_valid_element_echoes_it(self):
        self.assertEqual(self.op('a2'), 'a2')

    def test_one_arg_invalid_element_raises(self):
        with self.assertRaises(ValueError):
            self.op('not_an_element')

    def test_two_args_returns_product(self):
        self.assertEqual(self.op('a', 'a2'), 'a3')
        self.assertEqual(self.op('a3', 'a'), 'e')

    def test_two_args_non_associative_magma(self):
        self.assertEqual(self.rps_op('r', 'p'), 'p')
        self.assertEqual(self.rps_op('p', 's'), 's')

    def test_more_than_two_args_left_associates(self):
        # a * a * a * a = a2 * a * a = a3 * a = e
        self.assertEqual(self.op('a', 'a', 'a', 'a'), 'e')

    def test_more_than_two_args_matches_manual_left_association(self):
        manual = self.rps_op(self.rps_op(self.rps_op('r', 'p'), 's'), 'p')
        self.assertEqual(self.rps_op('r', 'p', 's', 'p'), manual)

    def test_callable_directly(self):
        # __call__ delegates to _op
        self.assertTrue(callable(self.op))
        self.assertEqual(self.op('a', 'a'), self.op._op('a', 'a'))
