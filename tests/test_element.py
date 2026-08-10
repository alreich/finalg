"""
Unit tests for finalg.element.Element.
"""

from unittest import TestCase

from finalg import generate_cyclic_group, generate_symmetric_group, generate_algebra_mod_n, make_finite_algebra
from finalg.element import Element


class TestElementConstruction(TestCase):

    def setUp(self):
        self.z4 = generate_cyclic_group(4)

    def test_valid_element(self):
        e = Element('1', self.z4)
        self.assertEqual(e.name, '1')
        self.assertIs(e.algebra, self.z4)

    def test_non_string_name_raises(self):
        with self.assertRaises(ValueError):
            Element(1, self.z4)

    def test_name_not_in_algebra_raises(self):
        with self.assertRaises(ValueError):
            Element('99', self.z4)

    def test_str_and_repr(self):
        e = Element('2', self.z4)
        self.assertEqual(str(e), '2')
        self.assertEqual(repr(e), "'2'")


class TestElementArithmeticOnGroup(TestCase):
    """Z4 is a Group: op is addition mod 4; has sub, inv, but no mult/div."""

    def setUp(self):
        self.z4 = generate_cyclic_group(4)

    def elem(self, name):
        return Element(name, self.z4)

    def test_add(self):
        self.assertEqual((self.elem('1') + self.elem('2')).name, '3')

    def test_sub_available_on_group(self):
        # Group defines .sub, so Element subtraction works even without a Ring.
        self.assertEqual((self.elem('1') - self.elem('2')).name, '3')  # 1 - 2 = -1 = 3 mod 4

    def test_neg(self):
        self.assertEqual((-self.elem('1')).name, '3')

    def test_mult_raises_when_algebra_has_no_mult(self):
        with self.assertRaises(ValueError):
            self.elem('1') * self.elem('2')

    def test_div_raises_when_algebra_has_no_div(self):
        with self.assertRaises(ValueError):
            self.elem('1') / self.elem('2')

    def test_pow_positive(self):
        self.assertEqual((self.elem('1') ** 3).name, '3')

    def test_pow_zero_returns_identity(self):
        self.assertEqual((self.elem('2') ** 0).name, '0')

    def test_pow_negative_uses_inverse(self):
        self.assertEqual((self.elem('1') ** -1).name, '3')

    def test_equality(self):
        self.assertEqual(self.elem('1'), Element('1', self.z4))

    def test_inequality_different_names(self):
        self.assertNotEqual(self.elem('1'), self.elem('2'))

    def test_equality_with_non_element_not_implemented(self):
        self.assertFalse(self.elem('1') == '1')

    def test_hash_consistent_with_equality(self):
        # Element.__hash__ was fixed to call __key() (previously hashed the bound
        # method itself). Equal Elements now hash equal.
        a = self.elem('1')
        b = Element('1', self.z4)
        self.assertEqual(a, b)
        self.assertEqual(hash(a), hash(b))


class TestElementArithmeticOnField(TestCase):
    """F5 is a Field: has add, mult, sub, div, inv, mult_inv."""

    def setUp(self):
        self.f5 = generate_algebra_mod_n(5)

    def elem(self, name):
        return Element(name, self.f5)

    def test_mult(self):
        self.assertEqual((self.elem('2') * self.elem('3')).name, '1')  # 6 mod 5 = 1

    def test_div(self):
        self.assertEqual((self.elem('2') / self.elem('3')).name, '4')  # 2 * inv(3) mod 5

    def test_add_and_sub(self):
        self.assertEqual((self.elem('2') + self.elem('4')).name, '1')
        self.assertEqual((self.elem('2') - self.elem('4')).name, '3')


class TestElementConjugation(TestCase):
    """S3 has a .conjugate method inherited from Group, so | is supported."""

    def setUp(self):
        self.s3 = generate_symmetric_group(3)

    def test_conjugate_matches_group_conjugate_method(self):
        a_name, b_name = self.s3.elements[1], self.s3.elements[2]
        a = Element(a_name, self.s3)
        b = Element(b_name, self.s3)
        result = a | b
        expected = self.s3.conjugate(a_name, b_name)
        self.assertEqual(result.name, expected)

    def test_conjugate_raises_when_unsupported(self):
        rps = make_finite_algebra('RPS', 'Rock Paper Scissors',
                                   ['r', 'p', 's'], [[0, 1, 0], [1, 1, 2], [0, 2, 2]])
        a = Element('r', rps)
        b = Element('p', rps)
        with self.assertRaises(ValueError):
            a | b


class TestElementOnNonAssociativeMagma(TestCase):

    def setUp(self):
        self.rps = make_finite_algebra('RPS', 'Rock Paper Scissors',
                                        ['r', 'p', 's'], [[0, 1, 0], [1, 1, 2], [0, 2, 2]])

    def test_add_uses_magma_op(self):
        a = Element('r', self.rps)
        b = Element('p', self.rps)
        self.assertEqual((a + b).name, self.rps.op('r', 'p'))

    def test_neg_without_inverses_raises(self):
        # FiniteAlgebra.inv() itself returns None when the algebra lacks inverses,
        # and Element.__neg__ then tries to build Element(None, algebra), which
        # raises ValueError (None is not a string) rather than yielding None.
        a = Element('r', self.rps)
        with self.assertRaises(ValueError):
            -a

    def test_pow_zero_without_identity_raises(self):
        a = Element('r', self.rps)
        with self.assertRaises(ValueError):
            a ** 0

    def test_pow_positive_ok(self):
        a = Element('r', self.rps)
        self.assertEqual((a ** 2).name, self.rps.op('r', 'r'))
