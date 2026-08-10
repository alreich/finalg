"""
Unit tests for finalg.quasigroup_and_loop.Quasigroup and Loop.
"""

from unittest import TestCase

from finalg import make_finite_algebra
from finalg.quasigroup_and_loop import Quasigroup, Loop


class TestQuasigroup(TestCase):

    def setUp(self):
        self.tbl9 = [[0, 4, 8, 2, 3, 9, 6, 7, 1, 5],
                     [3, 6, 2, 8, 7, 1, 9, 5, 0, 4],
                     [8, 9, 3, 1, 0, 6, 4, 2, 5, 7],
                     [1, 7, 6, 5, 4, 8, 0, 3, 2, 9],
                     [2, 1, 9, 0, 6, 7, 5, 8, 4, 3],
                     [5, 2, 7, 4, 9, 3, 1, 0, 8, 6],
                     [4, 3, 0, 6, 1, 5, 2, 9, 7, 8],
                     [9, 8, 5, 7, 2, 0, 3, 4, 6, 1],
                     [7, 0, 1, 9, 5, 4, 8, 6, 3, 2],
                     [6, 5, 4, 3, 8, 2, 7, 1, 9, 0]]
        self.elements = [str(i) for i in range(10)]
        self.qg = make_finite_algebra('QG9', 'A quasigroup', self.elements, self.tbl9)

    def test_make_finite_algebra_returns_quasigroup(self):
        self.assertEqual(type(self.qg).__name__, 'Quasigroup')
        self.assertIsInstance(self.qg, Quasigroup)

    def test_has_cancellation(self):
        self.assertTrue(self.qg.has_cancellation())

    def test_is_not_associative(self):
        self.assertFalse(self.qg.is_associative())

    def test_has_no_identity(self):
        self.assertFalse(self.qg.has_identity())

    def test_order(self):
        self.assertEqual(self.qg.order, 10)

    def test_direct_construction(self):
        qg2 = Quasigroup('QG9', 'A quasigroup', self.elements, self.tbl9)
        self.assertEqual(qg2, self.qg)

    def test_op_matches_table(self):
        self.assertEqual(self.qg.op('0', '1'), '4')
        self.assertEqual(self.qg.op('9', '9'), '0')


class TestLoop(TestCase):

    def setUp(self):
        self.loop_tbl = [[0, 1, 2, 3, 4, 5, 6],
                         [1, 2, 0, 5, 6, 4, 3],
                         [2, 0, 1, 6, 5, 3, 4],
                         [3, 6, 5, 4, 0, 1, 2],
                         [4, 5, 6, 0, 3, 2, 1],
                         [5, 3, 4, 2, 1, 6, 0],
                         [6, 4, 3, 1, 2, 0, 5]]
        self.elements = [str(i) for i in range(7)]
        self.loop = make_finite_algebra('L7', 'A loop', self.elements, self.loop_tbl)

    def test_make_finite_algebra_returns_loop(self):
        self.assertEqual(type(self.loop).__name__, 'Loop')
        self.assertIsInstance(self.loop, Loop)
        self.assertIsInstance(self.loop, Quasigroup)

    def test_has_identity(self):
        self.assertTrue(self.loop.has_identity())
        self.assertEqual(self.loop.identity, '0')

    def test_is_not_associative(self):
        self.assertFalse(self.loop.is_associative())

    def test_has_cancellation(self):
        self.assertTrue(self.loop.has_cancellation())

    def test_direct_construction(self):
        loop2 = Loop('L7', 'A loop', self.elements, self.loop_tbl)
        self.assertEqual(loop2, self.loop)

    def test_identity_behaves_as_identity_for_op(self):
        for elem in self.loop.elements:
            self.assertEqual(self.loop.op('0', elem), elem)
            self.assertEqual(self.loop.op(elem, '0'), elem)
