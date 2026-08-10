"""
Unit tests for finalg.module: Module, NDimensionalModule, and the module
condition-checking functions.
"""

import io
import contextlib
from unittest import TestCase

from finalg import generate_algebra_mod_n, generate_cyclic_group
from finalg.module import (
    Module,
    NDimensionalModule,
    module_sv_mult,
    module_dot_product,
    check_module_conditions,
    check_scaling_by_one,
    check_dist_of_scalars_over_vec_add,
    check_dist_of_vec_over_scalar_add,
    check_associativity,
)


class TestNDimensionalModule(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.r4 = generate_algebra_mod_n(4)
        cls.m3 = NDimensionalModule(cls.r4, 3)

    def test_type_and_dimensions(self):
        self.assertEqual(type(self.m3).__name__, 'NDimensionalModule')
        self.assertEqual(self.m3.dimensions, 3)

    def test_origin_is_all_zero_vector(self):
        self.assertEqual(self.m3.origin, '0:0:0')

    def test_scalar_and_vector_algebras(self):
        self.assertIs(self.m3.scalar, self.r4)
        self.assertEqual(self.m3.vector.order, self.r4.order ** 3)

    def test_vector_add(self):
        v1, v2 = self.m3.vector.elements[0], self.m3.vector.elements[1]
        self.assertEqual(self.m3.vector_add(v1, v2), self.m3.vector.op(v1, v2))

    def test_sv_mult(self):
        result = self.m3.sv_mult('2', '0:1:1')
        self.assertEqual(result, '0:2:2')

    def test_dot_product(self):
        # (0,0,1) . (0,0,2) = 0+0+2 = 2 mod 4
        self.assertEqual(self.m3.dot_product('0:0:1', '0:0:2'), '2')

    def test_repr(self):
        r = repr(self.m3)
        self.assertIn('NDimensionalModule', r)
        self.assertIn('Scalars:R4', r)

    def test_about_smoke(self):
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            result = self.m3.about()
        self.assertIsNone(result)
        out = buf.getvalue()
        self.assertIn('SCALARS', out)
        self.assertIn('VECTORS', out)


class TestModuleConstructionErrors(TestCase):

    def test_non_ring_scalar_raises(self):
        z4 = generate_cyclic_group(4)
        with self.assertRaises(ValueError):
            Module('bad', 'desc', z4, z4, lambda s, v: v)

    def test_invalid_sv_mult_raises(self):
        r4 = generate_algebra_mod_n(4)
        with self.assertRaises(ValueError):
            Module('bad', 'desc', r4, r4 ** 2, lambda s, v: v)  # doesn't satisfy module laws


class TestModuleSvMultAndDotProduct(TestCase):

    def test_module_sv_mult_scales_each_component(self):
        r4 = generate_algebra_mod_n(4)
        sv_mult = module_sv_mult(r4)
        self.assertEqual(sv_mult('2', '1:2:3'), '2:0:2')

    def test_module_dot_product(self):
        r4 = generate_algebra_mod_n(4)
        m3 = NDimensionalModule(r4, 3)
        result = module_dot_product(m3, '1:1:1', '1:1:1')
        self.assertEqual(result, '3')  # 1+1+1 mod 4


class TestModuleConditionCheckers(TestCase):

    def setUp(self):
        self.r4 = generate_algebra_mod_n(4)
        self.group = self.r4 ** 2
        self.sv_mult = module_sv_mult(self.r4)

    def test_check_scaling_by_one_true(self):
        self.assertTrue(check_scaling_by_one(self.r4, self.group, self.sv_mult))

    def test_check_dist_of_scalars_over_vec_add_true(self):
        self.assertTrue(check_dist_of_scalars_over_vec_add(self.r4, self.group, self.sv_mult))

    def test_check_dist_of_vec_over_scalar_add_true(self):
        self.assertTrue(check_dist_of_vec_over_scalar_add(self.r4, self.group, self.sv_mult))

    def test_check_associativity_true(self):
        self.assertTrue(check_associativity(self.r4, self.group, self.sv_mult))

    def test_check_module_conditions_true_for_valid_module(self):
        self.assertTrue(check_module_conditions(self.r4, self.group, self.sv_mult))

    def test_check_module_conditions_false_for_broken_sv_mult(self):
        broken_sv_mult = lambda s, v: v  # ignores the scalar entirely
        self.assertFalse(check_module_conditions(self.r4, self.group, broken_sv_mult))

    def test_verbose_smoke(self):
        buf = io.StringIO()
        broken_sv_mult = lambda s, v: v
        with contextlib.redirect_stdout(buf):
            result = check_module_conditions(self.r4, self.group, broken_sv_mult, verbose=True)
        self.assertFalse(result)
        self.assertIn('OK?', buf.getvalue())
