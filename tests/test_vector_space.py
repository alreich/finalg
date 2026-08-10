"""
Unit tests for finalg.vector_space: VectorSpace and NDimensionalVectorSpace.
"""

from unittest import TestCase

from finalg import generate_algebra_mod_n
from finalg.vector_space import VectorSpace, NDimensionalVectorSpace
from finalg.module import module_sv_mult, Module


class TestNDimensionalVectorSpace(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.f5 = generate_algebra_mod_n(5)
        cls.vs3 = NDimensionalVectorSpace(cls.f5, 3)

    def test_type_and_dimensions(self):
        self.assertEqual(type(self.vs3).__name__, 'NDimensionalVectorSpace')
        self.assertIsInstance(self.vs3, VectorSpace)
        self.assertIsInstance(self.vs3, Module)
        self.assertEqual(self.vs3.dimensions, 3)

    def test_origin(self):
        self.assertEqual(self.vs3.origin, '0:0:0')

    def test_scalar_is_a_field(self):
        self.assertIs(self.vs3.scalar, self.f5)

    def test_vector_order(self):
        self.assertEqual(self.vs3.vector.order, 5 ** 3)

    def test_dot_product(self):
        self.assertEqual(self.vs3.dot_product('0:0:1', '0:0:2'), '2')

    def test_sv_mult(self):
        self.assertEqual(self.vs3.sv_mult('2', '1:1:1'), '2:2:2')


class TestVectorSpaceConstructionErrors(TestCase):

    def test_ring_that_is_not_a_field_raises(self):
        # R4 (mod 4) is a Ring but not a Field, since 4 is composite.
        r4 = generate_algebra_mod_n(4)
        with self.assertRaises(ValueError):
            NDimensionalVectorSpace(r4, 2)

    def test_direct_construction_with_non_field_raises(self):
        r4 = generate_algebra_mod_n(4)
        with self.assertRaises(ValueError):
            VectorSpace('bad', 'desc', r4, r4 ** 2, module_sv_mult(r4))

    def test_direct_construction_with_field_succeeds(self):
        f5 = generate_algebra_mod_n(5)
        group = f5 ** 2
        vs = VectorSpace('VS', 'desc', f5, group, module_sv_mult(f5))
        self.assertIs(vs.scalar, f5)
        self.assertIs(vs.vector, group)
