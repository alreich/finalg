"""
Unit tests for finalg.infix_notation.InfixNotation.
"""

from unittest import TestCase

from finalg import generate_symmetric_group, generate_cyclic_group
from finalg.infix_notation import InfixNotation
from finalg.element import Element


class TestInfixNotation(TestCase):

    def setUp(self):
        self.s3 = generate_symmetric_group(3)
        self.z4 = generate_cyclic_group(4)

    def test_enter_returns_element_map(self):
        with InfixNotation(self.s3) as f:
            self.assertEqual(set(f.keys()), set(self.s3.elements))
            for name, elem in f.items():
                self.assertIsInstance(elem, Element)
                self.assertEqual(elem.name, name)

    def test_infix_arithmetic_matches_op(self):
        elems = self.s3.elements
        with InfixNotation(self.s3) as f:
            result = (f[elems[1]] + f[elems[2]]).name
        self.assertEqual(result, self.s3.op(elems[1], elems[2]))

    def test_context_manager_exit_does_not_suppress_exceptions(self):
        with self.assertRaises(ZeroDivisionError):
            with InfixNotation(self.z4):
                raise ZeroDivisionError("boom")

    def test_element_map_is_freshly_built_per_instantiation(self):
        with InfixNotation(self.z4) as f1:
            e1 = f1['1']
        with InfixNotation(self.z4) as f2:
            e2 = f2['1']
        self.assertEqual(e1, e2)
        self.assertIsNot(e1, e2)
