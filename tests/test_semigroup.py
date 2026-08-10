"""
Unit tests for finalg.semigroup.Semigroup.
"""

from unittest import TestCase

from finalg import make_finite_algebra
from finalg.semigroup import Semigroup


class TestSemigroup(TestCase):

    def setUp(self) -> None:
        self.ex141_tbl = [[0, 3, 0, 3, 0, 3],
                          [1, 4, 1, 4, 1, 4],
                          [2, 5, 2, 5, 2, 5],
                          [3, 0, 3, 0, 3, 0],
                          [4, 1, 4, 1, 4, 1],
                          [5, 2, 5, 2, 5, 2]]
        self.ex141_sg = make_finite_algebra('ex141', 'foobar',
                                            ['a', 'b', 'c', 'd', 'e', 'f'],
                                            self.ex141_tbl)

    def test_make_finite_algebra_returns_semigroup(self):
        self.assertIsInstance(self.ex141_sg, Semigroup)
        self.assertEqual(type(self.ex141_sg).__name__, 'Semigroup')

    def test_is_associative(self):
        self.assertTrue(self.ex141_sg.is_associative())

    def test_is_not_commutative(self):
        self.assertFalse(self.ex141_sg.is_commutative())

    def test_has_no_identity(self):
        self.assertIsNone(self.ex141_sg.identity)
        self.assertFalse(self.ex141_sg.has_identity())

    def test_associativity_directly_via_op(self):
        ab = self.ex141_sg.op('a', 'b')
        bc = self.ex141_sg.op('b', 'c')
        ab_c = self.ex141_sg.op(ab, 'c')
        a_bc = self.ex141_sg.op('a', bc)
        self.assertEqual(ab_c, a_bc)

    def test_equality(self):
        sg2 = make_finite_algebra(
            'ex141', 'foobar', ['a', 'b', 'c', 'd', 'e', 'f'], self.ex141_tbl)
        self.assertEqual(self.ex141_sg, sg2)

    def test_direct_construction_raises_when_not_associative(self):
        # RPS is a Magma, not a Semigroup.
        with self.assertRaises(ValueError):
            Semigroup('bad', 'not associative', ['r', 'p', 's'],
                      [[0, 1, 0], [1, 1, 2], [0, 2, 2]])

    def test_direct_construction_skips_check_when_disabled(self):
        # With check_inputs=False, even a non-associative table is accepted
        # (this is used internally by make_finite_algebra's dispatch logic
        # is NOT exercised this way in practice, but the flag itself must work).
        sg = Semigroup('bad', 'not associative', ['r', 'p', 's'],
                        [[0, 1, 0], [1, 1, 2], [0, 2, 2]], check_inputs=False)
        self.assertFalse(sg.is_associative())

    def test_is_regular(self):
        self.assertTrue(self.ex141_sg.is_regular())

    def test_weak_inverses(self):
        result = self.ex141_sg.weak_inverses()
        self.assertEqual(result,
                          {'a': ['a', 'c', 'e'], 'b': ['b', 'd', 'f'],
                           'c': ['a', 'c', 'e'], 'd': ['b', 'd', 'f'],
                           'e': ['a', 'c', 'e'], 'f': ['b', 'd', 'f']})

    def test_weak_inverses_keys_match_elements(self):
        result = self.ex141_sg.weak_inverses()
        self.assertEqual(set(result.keys()), set(self.ex141_sg.elements))

    def test_is_regular_for_a_group_is_always_true(self):
        # Every element of a group is its own weak inverse's partner (via true
        # inverse), so groups are trivially regular.
        from finalg import generate_cyclic_group
        z4 = generate_cyclic_group(4)
        self.assertTrue(z4.is_regular())
