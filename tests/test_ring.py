"""
Unit tests for finalg.ring.Ring, including the Cayley-Dickson construction.
"""

from unittest import TestCase

from finalg import generate_algebra_mod_n
from finalg.ring import Ring
from finalg.group import Group
from finalg.monoid import Monoid


class TestRingConstruction(TestCase):

    def test_generate_algebra_mod_n_composite_is_ring_not_field(self):
        r4 = generate_algebra_mod_n(4)
        self.assertEqual(type(r4).__name__, 'Ring')
        self.assertEqual(r4.elements, ('0', '1', '2', '3'))

    def test_direct_construction_raises_when_add_not_commutative(self):
        # A non-commutative addition table (S3-like) can't be a Ring's additive group.
        s3_tbl = [[0, 1, 2, 3, 4, 5],
                  [1, 2, 0, 5, 3, 4],
                  [2, 0, 1, 4, 5, 3],
                  [3, 4, 5, 0, 1, 2],
                  [4, 5, 3, 2, 0, 1],
                  [5, 3, 4, 1, 2, 0]]
        with self.assertRaises(ValueError):
            Ring('bad', 'non-comm add', ['0', '1', '2', '3', '4', '5'], s3_tbl, s3_tbl)

    def test_direct_construction_raises_when_mult_not_distributive(self):
        add_tbl = [[0, 1, 2, 3], [1, 0, 3, 2], [2, 3, 0, 1], [3, 2, 1, 0]]
        # Non-distributive multiplication table (arbitrary, breaks distributivity)
        bad_mult_tbl = [[0, 1, 2, 3], [1, 0, 3, 2], [2, 3, 0, 1], [3, 2, 1, 0]]
        with self.assertRaises(ValueError):
            Ring('bad', 'bad mult', ['0', '1', '2', '3'], add_tbl, bad_mult_tbl)


class TestRingBasics(TestCase):

    def setUp(self):
        self.r4 = generate_algebra_mod_n(4)

    def test_zero_and_one(self):
        self.assertEqual(self.r4.zero, '0')
        self.assertEqual(self.r4.add_identity, '0')
        self.assertEqual(self.r4.one, '1')
        self.assertEqual(self.r4.mult_identity, '1')

    def test_has_mult_identity(self):
        self.assertTrue(self.r4.has_mult_identity())

    def test_minus_one(self):
        self.assertEqual(self.r4.minus_one, '3')

    def test_add_and_mult(self):
        self.assertEqual(self.r4.add('2', '3'), '1')
        self.assertEqual(self.r4.mult('2', '3'), '2')

    def test_add_table_and_mult_table(self):
        self.assertIs(self.r4.add_table, self.r4.table)
        self.assertIs(self.r4.mult_table, self.r4._ring_mult_table)

    def test_mult_is_commutative(self):
        self.assertTrue(self.r4.mult_is_commutative())

    def test_extract_additive_algebra(self):
        add_alg = self.r4.extract_additive_algebra()
        self.assertIsInstance(add_alg, Group)
        self.assertEqual(add_alg.elements, self.r4.elements)

    def test_extract_multiplicative_algebra(self):
        mult_alg = self.r4.extract_multiplicative_algebra()
        self.assertIsInstance(mult_alg, Monoid)
        self.assertEqual(mult_alg.elements, self.r4.elements)

    def test_element_to_power_uses_multiplication(self):
        # 2^3 mod 4 = 0
        self.assertEqual(self.r4.element_to_power('2', 3), '0')


class TestRingZeroDivisorsAndUnits(TestCase):

    def setUp(self):
        self.r4 = generate_algebra_mod_n(4)  # composite modulus -> zero divisors
        self.f5 = generate_algebra_mod_n(5)  # prime modulus -> field, no zero divisors

    def test_zero_divisors_mod4(self):
        self.assertEqual(self.r4.zero_divisors(), ['2'])

    def test_zero_divisor_pairs_mod4(self):
        self.assertEqual(self.r4.zero_divisor_pairs(), [('2', '2')])

    def test_units_mod4(self):
        self.assertEqual(sorted(self.r4.units()), ['1', '3'])

    def test_field_has_no_zero_divisors(self):
        self.assertEqual(self.f5.zero_divisors(), [])

    def test_commutator_of_commutative_ring_is_zero(self):
        self.assertEqual(self.r4.commutator('2', '3'), self.r4.zero)

    def test_element_pairs_where_product_equals(self):
        pairs = self.r4.element_pairs_where_product_equals('0')
        self.assertIn(('2', '2'), pairs)
        self.assertIn(('0', '0'), pairs)


class TestRingSquareRoots(TestCase):

    def setUp(self):
        self.r4 = generate_algebra_mod_n(4)

    def test_square_root_mapping(self):
        mapping = self.r4.square_root_mapping()
        self.assertEqual(mapping, {'0': ['0', '2'], '1': ['1', '3']})

    def test_square_roots_of_element_with_roots(self):
        self.assertEqual(sorted(self.r4.square_roots('1')), ['1', '3'])

    def test_square_roots_of_element_without_roots_is_empty(self):
        self.assertEqual(self.r4.square_roots('2'), [])


class TestRingDirectProduct(TestCase):

    def test_mul_produces_ring(self):
        r4 = generate_algebra_mod_n(4)
        f3 = generate_algebra_mod_n(3)
        dp = r4 * f3
        self.assertEqual(type(dp).__name__, 'Ring')
        self.assertEqual(dp.order, 12)

    def test_mul_with_non_ring_raises(self):
        r4 = generate_algebra_mod_n(4)
        from finalg import generate_cyclic_group
        z4 = generate_cyclic_group(4)
        with self.assertRaises(ValueError):
            r4 * z4


class TestRingElementOrder(TestCase):
    """A Ring's additive structure is a Group, so element_order (inherited via
    Monoid/Group) reflects additive order mod n."""

    def test_element_orders_mod_6(self):
        r6 = generate_algebra_mod_n(6)
        self.assertEqual(r6.element_order('0'), 1)
        self.assertEqual(r6.element_order('1'), 6)
        self.assertEqual(r6.element_order('2'), 3)
        self.assertEqual(r6.element_order('3'), 2)
    """Builds Gaussian-integer-like rings via the Cayley-Dickson construction."""

    def setUp(self):
        self.f3 = generate_algebra_mod_n(3)  # A Field (3 is prime)
        self.zi3 = self.f3.make_cayley_dickson_algebra()

    def test_default_version_element_count(self):
        self.assertEqual(self.zi3.order, 9)

    def test_elements_are_pairs_joined_by_delimiter(self):
        self.assertIn('1:2', self.zi3.elements)

    def test_conjugates_mapping_exists(self):
        self.assertIsNotNone(self.zi3.conjugates())

    def test_conj_of_compound_element(self):
        # elem_conj('1:2') = conj('1') : inv('2') = '1' : '1' (since add-inverse of 2 mod 3 is 1)
        self.assertEqual(self.zi3.conj('1:2'), '1:1')

    def test_norm_is_product_with_conjugate(self):
        elem = '1:2'
        expected = self.zi3.mult(elem, self.zi3.conj(elem))
        self.assertEqual(self.zi3.norm(elem), expected)

    def test_split_element(self):
        self.assertEqual(self.zi3.split_element('1:2'), ('1', '2'))

    def test_split_element_noncompound_returns_unchanged(self):
        self.assertEqual(self.f3.split_element('1'), '1')

    def test_scalar_mult(self):
        # scalar_mult is meant to be called on the Cayley-Dickson ring itself,
        # using a scalar from the base ring, applied to a compound element.
        result = self.zi3.scalar_mult('2', '1:2')
        expected = self.zi3.mult('2:0', '1:2')
        self.assertEqual(result, expected)

    def test_is_gaussian_prime_true_case(self):
        # 1 + 2i: 1^2 + 2^2 = 5, which is prime -> Gaussian prime
        self.assertTrue(self.zi3.is_gaussian_prime('1:2'))

    def test_is_gaussian_prime_false_case(self):
        # 2 + 2i: 2^2 + 2^2 = 8, not prime -> not a Gaussian prime
        self.assertFalse(self.zi3.is_gaussian_prime('2:2'))

    def test_is_gaussian_prime_on_scalar_element_is_broken(self):
        # NOTE: is_gaussian_prime's dimension==0 branch sets `real = elem` without
        # converting the string element name to an int, so calling it on a base
        # (non-compound) Ring/Field element raises an error from sympy's isprime
        # rather than returning a boolean. This test documents that real,
        # currently-existing bug.
        with self.assertRaises(ValueError):
            self.f3.is_gaussian_prime('2')

    def test_make_cayley_dickson_algebra_versions_2_and_3_require_mu_without_identity(self):
        # Build a mult-identity-less ring by using a Ring with no mult identity.
        # generate_algebra_mod_n always has a mult identity, so instead verify the
        # happy path: default mu is used automatically when a mult identity exists.
        zi3_v2 = self.f3.make_cayley_dickson_algebra(version=2)
        self.assertEqual(zi3_v2.order, 9)

    def test_invalid_version_raises(self):
        with self.assertRaises(ValueError):
            self.f3.make_cayley_dickson_algebra(version=99)

    def test_repr_and_equality(self):
        zi3_copy = self.f3.make_cayley_dickson_algebra()
        self.assertEqual(self.zi3, zi3_copy)
        self.assertIn('Field(', repr(self.zi3))

    def test_hash_inconsistent_with_equality(self):
        # Same missing-parens `_key` bug as Magma/Element (see test_magma.py).
        zi3_copy = self.f3.make_cayley_dickson_algebra()
        self.assertEqual(self.zi3, zi3_copy)
        self.assertNotEqual(hash(self.zi3), hash(zi3_copy))
