"""
Unit tests for finalg.abstract_matrix.AbstractMatrix.

@author: Alfred J. Reich
"""

from unittest import TestCase
import numpy as np

import finalg as alg
import finalg.abstract_matrix as mat
from finalg.ring import Ring


class TestAbstractMatrix(TestCase):

    def setUp(self) -> None:

        # The following duplication is for equality testing
        self.psr3 = alg.generate_powerset_ring(3)
        self.psr3x = alg.generate_powerset_ring(3)
        self.arr1 = [['{0, 1, 2}', '{2}', '{2}'],
                     ['{0, 1, 2}', '{2}', '{1, 2}'],
                     ['{0}', '{0, 1, 2}', '{0, 1, 2}']]
        self.arr1x = [['{0, 1, 2}', '{2}', '{2}'],
                      ['{0, 1, 2}', '{2}', '{1, 2}'],
                      ['{0}', '{0, 1, 2}', '{0, 1, 2}']]
        self.arr1t = [['{0, 1, 2}', '{0, 1, 2}', '{0}'],
                      ['{2}', '{2}', '{0, 1, 2}'],
                      ['{2}', '{1, 2}', '{0, 1, 2}']]
        self.arr2 = [['{0, 2}', '{2}', '{2}'],
                     ['{0, 2}', '{2}', '{2}'],
                     ['{0}', '{0, 2}', '{0, 2}']]
        self.arr3 = [['{0, 1, 2}', '{2}'],
                     ['{0}', '{0, 1, 2}']]
        self.arr4 = [['{1}', '{0, 1, 2}', '{0, 1, 2}'],
                     ['{}', '{0, 1, 2}', '{0, 1, 2}'],
                     ['{}', '{1}', '{}']]
        self.arr5 = [['{0, 1, 2}', '{0, 1, 2}', '{0}'],
                     ['{0, 1, 2}', '{0, 2}', '{0, 1}'],
                     ['{0}', '{0, 1}', '{0}']]
        self.arr6 = [['{}', '{0, 1}', '{0, 2}'],
                     ['{0, 1}', '{}', '{0}'],
                     ['{0, 2}', '{0}', '{}']]
        self.arr7 = [['{}', '{0, 1}', '{0, 2}'],
                     ['{0, 1}', '{}', '{0}'],
                     ['{0, 2}', '{0}', '{}']]
        self.mat1 = mat.AbstractMatrix(self.arr1, self.psr3)
        self.mat1x = mat.AbstractMatrix(self.arr1x, self.psr3x)
        self.mat1t = mat.AbstractMatrix(self.arr1t, self.psr3)
        self.mat2 = mat.AbstractMatrix(self.arr2, self.psr3)
        self.mat1minor = mat.AbstractMatrix(self.arr3, self.psr3)
        self.mat1cof = mat.AbstractMatrix(self.arr4, self.psr3)
        self.mat1xmat1t = mat.AbstractMatrix(self.arr5, self.psr3)
        self.mat1pmat1t = mat.AbstractMatrix(self.arr6, self.psr3)
        self.mat1mmat1t = mat.AbstractMatrix(self.arr7, self.psr3)

    def test_equality_of_matrices(self):
        self.assertTrue(self.mat1 == self.mat1x)

    def test_inequality_different_ring_instance_still_equal_by_value(self):
        # Ring equality is value-based (via __eq__), so a matrix built with an
        # independently-constructed but value-equal Ring still compares equal.
        self.assertEqual(self.mat1, self.mat1x)

    def test_inequality_different_content(self):
        self.assertNotEqual(self.mat1, self.mat2)

    def test_equality_with_non_matrix_not_implemented(self):
        self.assertFalse(self.mat1 == "not a matrix")

    def test_zeros(self):
        self.assertEqual(mat.AbstractMatrix.zeros((3, 4), self.psr3).array.tolist(),
                         [['{}', '{}', '{}', '{}'],
                          ['{}', '{}', '{}', '{}'],
                          ['{}', '{}', '{}', '{}']])

    def test_identity(self):
        self.assertEqual(mat.AbstractMatrix.identity(3, self.psr3).array.tolist(),
                         [['{0, 1, 2}', '{}', '{}'],
                          ['{}', '{0, 1, 2}', '{}'],
                          ['{}', '{}', '{0, 1, 2}']])

    def test_identity_returns_none_without_mult_identity(self):
        add_tbl = [[0, 1, 2], [1, 2, 0], [2, 0, 1]]
        mult_tbl = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]  # trivial mult: no identity
        zero_mult_ring = Ring('ZeroMultRing', 'desc', ['0', '1', '2'], add_tbl, mult_tbl)
        self.assertFalse(zero_mult_ring.has_mult_identity())
        self.assertIsNone(mat.AbstractMatrix.identity(2, zero_mult_ring))

    def test_random_shape_and_membership(self):
        np.random.seed(1)
        rand_mat = mat.AbstractMatrix.random((3, 4), self.psr3)
        self.assertEqual(rand_mat.shape, (3, 4))
        for row in rand_mat.array.tolist():
            for elem in row:
                self.assertIn(elem, self.psr3.elements)

    def test_array(self):
        self.assertEqual(self.mat1.array.tolist(), self.arr1)

    def test_shape(self):
        self.assertEqual(self.mat1.shape, (3, 3))

    def test_nrows(self):
        self.assertEqual(self.mat1.nrows, 3)

    def test_ncols(self):
        self.assertEqual(self.mat1.ncols, 3)

    def test_algebra(self):
        self.assertTrue(self.mat1.algebra == self.psr3x)

    def test_ring_property_same_as_algebra(self):
        self.assertIs(self.mat1.ring, self.mat1.algebra)

    def test_copy(self):
        copy = self.mat1.copy()
        self.assertEqual(self.mat1, copy)
        self.assertIsNot(self.mat1.array, copy.array)  # genuinely a new array

    def test_transpose(self):
        self.assertEqual(self.mat1.transpose(), self.mat1t)

    def test_scalar_mult_left(self):
        self.assertEqual(self.mat1.scalar_mult('{0, 2}'), self.mat2)

    def test_scalar_mult_right_differs_when_noncommutative(self):
        # Powerset-ring multiplication (intersection) is commutative, so left and
        # right scalar multiplication should agree here.
        left = self.mat1.scalar_mult('{0, 2}', left=True)
        right = self.mat1.scalar_mult('{0, 2}', left=False)
        self.assertEqual(left, right)

    def test_scalar_mult_invalid_scalar_raises(self):
        with self.assertRaises(ValueError):
            self.mat1.scalar_mult('not_an_element')

    def test_minor(self):
        self.assertEqual(self.mat1.minor(1, 1), self.mat1minor)

    def test_minor_invalid_row_raises(self):
        with self.assertRaises(ValueError):
            self.mat1.minor(99, 0)

    def test_minor_invalid_col_raises(self):
        with self.assertRaises(ValueError):
            self.mat1.minor(0, 99)

    def test_determinant(self):
        self.assertEqual(self.mat1.determinant(), '{1}')

    def test_determinant_1x1(self):
        m = mat.AbstractMatrix([['{0, 1}']], self.psr3)
        self.assertEqual(m.determinant(), '{0, 1}')

    def test_determinant_non_square_raises(self):
        m = mat.AbstractMatrix(self.arr3, self.psr3)  # 2x2, but let's make it non-square
        non_square = mat.AbstractMatrix([['{}', '{}', '{}'], ['{}', '{}', '{}']], self.psr3)
        with self.assertRaises(ValueError):
            non_square.determinant()

    def test_cofactor_matrix(self):
        self.assertEqual(self.mat1.cofactor_matrix(), self.mat1cof)

    def test_multiplication(self):
        self.assertEqual(self.mat1 * self.mat1.transpose(), self.mat1xmat1t)

    def test_multiplication_incompatible_shapes_raises(self):
        m1 = mat.AbstractMatrix([['{}', '{}', '{}'], ['{}', '{}', '{}']], self.psr3)  # 2x3
        m2 = mat.AbstractMatrix([['{}', '{}'], ['{}', '{}']], self.psr3)  # 2x2
        with self.assertRaises(ValueError):
            m1 * m2

    def test_multiplication_different_rings_raises(self):
        different_ring = alg.generate_algebra_mod_n(4)  # genuinely different table/elements
        m2 = mat.AbstractMatrix(self.arr1, different_ring)
        with self.assertRaises(ValueError):
            self.mat1 * m2

    def test_addition(self):
        self.assertEqual(self.mat1 + self.mat1.transpose(), self.mat1pmat1t)

    def test_addition_incompatible_shapes_raises(self):
        m2x3 = mat.AbstractMatrix([['{}', '{}', '{}'], ['{}', '{}', '{}']], self.psr3)
        with self.assertRaises(ValueError):
            self.mat1 + m2x3

    def test_addition_different_rings_raises(self):
        different_ring = alg.generate_algebra_mod_n(4)  # genuinely different table/elements
        m2 = mat.AbstractMatrix(self.arr1, different_ring)
        with self.assertRaises(ValueError):
            self.mat1 + m2

    def test_subtraction(self):
        self.assertEqual(self.mat1 - self.mat1.transpose(), self.mat1mmat1t)

    def test_subtraction_incompatible_shapes_raises(self):
        m2x3 = mat.AbstractMatrix([['{}', '{}', '{}'], ['{}', '{}', '{}']], self.psr3)
        with self.assertRaises(ValueError):
            self.mat1 - m2x3

    def test_getitem(self):
        self.assertEqual(self.mat1[1, 2], '{1, 2}')

    def test_setitem(self):
        self.mat1[1, 2] = '{}'
        self.assertEqual(self.mat1[1, 2], '{}')
        self.mat1[1, 2] = '{1, 2}'  # Put original value back

    def test_setitem_invalid_value_raises(self):
        with self.assertRaises(ValueError):
            self.mat1[0, 0] = 'not_an_element'

    def test_str(self):
        self.assertEqual(str(self.mat1), str(self.mat1.array))

    def test_repr_is_array2string(self):
        self.assertEqual(repr(self.mat1), np.array2string(self.mat1.array, separator=', '))

    def test_hash_consistent_with_equality(self):
        # AbstractMatrix.__hash__ was fixed to call __key() (previously hashed the
        # bound method itself), and __key() now returns nested tuples so it's
        # actually hashable. Equal matrices now hash equal.
        self.assertEqual(self.mat1, self.mat1x)
        self.assertEqual(hash(self.mat1), hash(self.mat1x))

    def test_to_tuple(self):
        result = self.mat1.to_tuple()
        self.assertEqual(result, tuple(tuple(row) for row in self.arr1))
        # Verify it's usable as a dict key.
        d = {result: 'ok'}
        self.assertEqual(d[result], 'ok')


class TestAbstractMatrixInverseOverField(TestCase):
    """Inverse only makes full sense over a Field, where every nonzero element
    (in particular, the determinant) has a multiplicative inverse."""

    def setUp(self):
        self.f5 = alg.generate_algebra_mod_n(5)
        self.m = mat.AbstractMatrix([['1', '2'], ['3', '4']], self.f5)

    def test_determinant(self):
        # det = 1*4 - 2*3 = 4 - 6 = -2 = 3 mod 5
        self.assertEqual(self.m.determinant(), '3')

    def test_inverse_times_original_is_identity(self):
        inv = self.m.inverse()
        product = self.m * inv
        expected_identity = mat.AbstractMatrix.identity(2, self.f5)
        self.assertEqual(product, expected_identity)
