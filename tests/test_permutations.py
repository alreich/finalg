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


class TestPermConstruction(TestCase):

    def test_base_0_valid(self):
        p = Perm((0, 1, 2))
        self.assertEqual(p.base, 0)

    def test_base_1_valid(self):
        p = Perm((1, 2, 3))
        self.assertEqual(p.base, 1)

    def test_gap_in_values_raises(self):
        with self.assertRaises(ValueError):
            Perm((1, 2, 4))

    def test_duplicate_values_raises(self):
        with self.assertRaises(ValueError):
            Perm((0, 0, 1))

    def test_base_other_than_0_or_1_raises(self):
        with self.assertRaises(ValueError):
            Perm((2, 3, 4))  # min value 2 is neither 0 nor 1


class TestPermCall(TestCase):

    def setUp(self):
        self.p = Perm((4, 2, 1, 5, 3))

    def test_call_on_string(self):
        self.assertEqual(self.p("ABCDE"), "DBAEC")

    def test_call_on_tuple(self):
        self.assertEqual(self.p((1, 2, 3, 4, 5)), (4, 2, 1, 5, 3))

    def test_call_on_list(self):
        self.assertEqual(self.p([1, 2, 3, 4, 5]), [4, 2, 1, 5, 3])

    def test_call_on_range_returns_list(self):
        self.assertEqual(self.p(range(1, 6)), [4, 2, 1, 5, 3])

    def test_call_on_perm_returns_perm(self):
        result = self.p(Perm((1, 2, 3, 4, 5)))
        self.assertIsInstance(result, Perm)
        self.assertEqual(result.values, (4, 2, 1, 5, 3))

    def test_call_wrong_length_raises(self):
        with self.assertRaises(ValueError):
            self.p("ABC")


class TestPermProperties(TestCase):

    def setUp(self):
        self.p = Perm((4, 2, 1, 5, 3))

    def test_values(self):
        self.assertEqual(self.p.values, (4, 2, 1, 5, 3))

    def test_mapping_is_zero_based(self):
        self.assertEqual(self.p.mapping, [3, 1, 0, 4, 2])

    def test_size_equals_len(self):
        self.assertEqual(self.p.size, len(self.p))
        self.assertEqual(self.p.size, 5)

    def test_getitem(self):
        self.assertEqual(self.p[0], 4)
        self.assertEqual(self.p[4], 3)

    def test_inverse(self):
        inv = self.p.inverse()
        # Composing a permutation with its inverse should give the identity.
        self.assertEqual((self.p * inv).values, self.p.id().values)
        self.assertEqual((inv * self.p).values, self.p.id().values)

    def test_id_matches_size_and_base(self):
        ident = self.p.id()
        self.assertEqual(ident.values, (1, 2, 3, 4, 5))
        self.assertEqual(ident.base, self.p.base)

    def test_str_and_repr(self):
        self.assertEqual(str(self.p), "(1 4 5 3)")
        self.assertEqual(repr(self.p), "Perm((4, 2, 1, 5, 3))")

    def test_hash_consistent_with_equality(self):
        p2 = Perm((4, 2, 1, 5, 3))
        self.assertEqual(self.p, p2)
        self.assertEqual(hash(self.p), hash(p2))

    def test_not_equal_operator(self):
        p2 = Perm((1, 2, 3, 4, 5))
        self.assertTrue(self.p != p2)
        self.assertFalse(self.p != Perm((4, 2, 1, 5, 3)))


class TestPermStaticMethods(TestCase):

    def test_identity_base_0(self):
        ident = Perm.identity(4)
        self.assertEqual(ident.values, (0, 1, 2, 3))

    def test_identity_base_1(self):
        ident = Perm.identity(4, base=1)
        self.assertEqual(ident.values, (1, 2, 3, 4))

    def test_identity_zero_size_raises(self):
        with self.assertRaises(ValueError):
            Perm.identity(0)

    def test_random_default_base(self):
        p = Perm.random(6)
        self.assertEqual(sorted(p.values), list(range(6)))

    def test_random_base_1(self):
        p = Perm.random(6, base=1)
        self.assertEqual(sorted(p.values), list(range(1, 7)))

    def test_random_invalid_base_raises(self):
        with self.assertRaises(ValueError):
            Perm.random(4, base=5)

    def test_from_cycles_base_0(self):
        p = Perm.from_cycles([[0, 2, 1], [3, 4]])
        self.assertEqual(p.values, (2, 0, 1, 4, 3))

    def test_from_cycles_base_1(self):
        p = Perm.from_cycles([[1, 3, 2], [4, 5]], base=1)
        self.assertEqual(p.values, (3, 1, 2, 5, 4))


class TestPermMultiplicationErrors(TestCase):

    def test_different_sizes_raises(self):
        p1 = Perm((0, 1, 2))
        p2 = Perm((0, 1, 2, 3))
        with self.assertRaises(ValueError):
            p1 * p2

    def test_different_bases_raises(self):
        p1 = Perm((0, 1, 2))
        p2 = Perm((1, 2, 3))
        with self.assertRaises(ValueError):
            p1 * p2


class TestPermParity(TestCase):

    def test_identity_is_even(self):
        ident = Perm((0, 1, 2))
        self.assertTrue(ident.is_even)
        self.assertEqual(ident.sign, 1)
        self.assertEqual(ident.parity, "even")

    def test_single_transposition_is_odd(self):
        p = Perm((1, 0, 2))
        self.assertFalse(p.is_even)
        self.assertEqual(p.sign, -1)
        self.assertEqual(p.parity, "odd")

    def test_three_cycle_is_even(self):
        p = Perm((1, 2, 0))
        self.assertTrue(p.is_even)
        self.assertEqual(p.sign, 1)

    def test_is_even_is_memoized(self):
        p = Perm((1, 0, 2))
        first = p.is_even
        second = p.is_even  # exercises the cached branch
        self.assertEqual(first, second)
        self.assertFalse(second)

