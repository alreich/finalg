"""
@author: Alfred J. Reich

"""

from unittest import TestCase
from finalg.permutation import Perm


class TestPermutations(TestCase):

    def setUp(self) -> None:

        # The examples here are from [Pinter, 1990], pages 70-71.

        self.s3_1 = {'epsilon': Perm((1, 2, 3)),
                     'alpha'  : Perm((1, 3, 2)),
                     'beta'   : Perm((3, 1, 2)),
                     'gamma'  : Perm((2, 1, 3)),
                     'delta'  : Perm((2, 3, 1)),
                     'kappa'  : Perm((3, 2, 1))}

        self.s3_0 = {'epsilon': Perm((0, 1, 2)),
                     'alpha'  : Perm((0, 2, 1)),
                     'beta'   : Perm((2, 0, 1)),
                     'gamma'  : Perm((1, 0, 2)),
                     'delta'  : Perm((1, 2, 0)),
                     'kappa'  : Perm((2, 1, 0))}

        # A reverse lookup dictionary, so that names can be looked up by permutation.
        # This is used when creating the multiplication table.
        self.s3_1_rev = {val: key for key, val in self.s3_1.items()}
        self.s3_0_rev = {val: key for key, val in self.s3_0.items()}

        self.a1 = self.s3_1['alpha']
        self.b1 = self.s3_1['beta']
        self.g1 = self.s3_1['gamma']

        self.a0 = self.s3_0['alpha']
        self.b0 = self.s3_0['beta']
        self.g0 = self.s3_0['gamma']

        self.p1 = Perm((0, 1, 2, 4, 3))  # [(3, 4)]
        self.p2 = Perm((1, 0, 3, 2, 4))  # [(0, 1), (2, 3), (4,)]
        self.p3 = Perm((4, 2, 0, 1, 3))  # [(0, 4, 3, 1, 2)]
        self.p4 = Perm((0, 2, 1, 4, 3))  # [(1, 2), (3, 4)]

    def test_to_cycles(self):
        self.assertEqual(self.p1.to_cycles(), [(0,), (3, 4)])
        self.assertEqual(self.p2.to_cycles(), [(0, 1), (2, 3), (4,)])
        self.assertEqual(self.p3.to_cycles(), [(0, 4, 3, 1, 2)])
        self.assertEqual(self.p4.to_cycles(), [(0,), (1, 2), (3, 4)])

    def test_from_cycles(self):
        self.assertEqual(self.p1, Perm.from_cycles([(3, 4)]))
        self.assertEqual(self.p2, Perm.from_cycles([(0, 1), (2, 3), (4,)]))
        self.assertEqual(self.p3, Perm.from_cycles([(0, 4, 3, 1, 2)]))
        self.assertEqual(self.p4, Perm.from_cycles([(1, 2), (3, 4)]))

    # Perm((0, 2, 1)) o Perm((2, 0, 1)) = Perm((1, 0, 2))
    # (i.e., alpha o beta = gamma)
    def test_multiplication(self):
        self.assertEqual(self.a1 * self.b1, self.g1)
        self.assertEqual(self.a0 * self.b0, self.g0)

    def test_mult_table_for_s3_1(self):
        expect = [['epsilon', 'alpha', 'beta', 'gamma', 'delta', 'kappa'],
                  ['alpha', 'epsilon', 'gamma', 'beta', 'kappa', 'delta'],
                  ['beta', 'kappa', 'delta', 'alpha', 'epsilon', 'gamma'],
                  ['gamma', 'delta', 'kappa', 'epsilon', 'alpha', 'beta'],
                  ['delta', 'gamma', 'epsilon', 'kappa', 'beta', 'alpha'],
                  ['kappa', 'beta', 'alpha', 'delta', 'gamma', 'epsilon']]
        s3_1_mul_tbl = [[self.s3_1_rev[self.s3_1[a] * self.s3_1[b]] for b in self.s3_1] for a in self.s3_1]
        self.assertEqual(s3_1_mul_tbl, expect)

    def test_mult_table_for_s3_0(self):
        expect = [['epsilon', 'alpha', 'beta', 'gamma', 'delta', 'kappa'],
                  ['alpha', 'epsilon', 'gamma', 'beta', 'kappa', 'delta'],
                  ['beta', 'kappa', 'delta', 'alpha', 'epsilon', 'gamma'],
                  ['gamma', 'delta', 'kappa', 'epsilon', 'alpha', 'beta'],
                  ['delta', 'gamma', 'epsilon', 'kappa', 'beta', 'alpha'],
                  ['kappa', 'beta', 'alpha', 'delta', 'gamma', 'epsilon']]
        s3_0_mul_tbl = [[self.s3_0_rev[self.s3_0[a] * self.s3_0[b]] for b in self.s3_0] for a in self.s3_0]
        self.assertEqual(s3_0_mul_tbl, expect)

