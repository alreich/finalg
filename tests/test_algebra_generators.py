"""
Unit tests for finalg.algebra_generators: functions that construct common
finite algebras (cyclic groups, symmetric groups, rings mod n, etc.).
"""

from unittest import TestCase

from finalg.algebra_generators import (
    generate_cyclic_group,
    generate_symmetric_group,
    generate_powerset_group,
    generate_commutative_monoid,
    generate_relative_primes_group,
    generate_powerset_ring,
    generate_algebra_mod_n,
    generate_nxn_matrix_algebra,
    generate_dihedral_group,
    generate_algebra_from_element_dict,
)
from finalg.group import Group
from finalg.ring import Ring
from finalg.field import Field
from finalg.monoid import Monoid
from finalg.abstract_matrix import AbstractMatrix


class TestGenerateCyclicGroup(TestCase):

    def test_default_naming(self):
        z5 = generate_cyclic_group(5)
        self.assertEqual(z5.name, 'Z5')
        self.assertEqual(z5.elements, ('0', '1', '2', '3', '4'))
        self.assertIsInstance(z5, Group)

    def test_custom_name_and_description(self):
        z5 = generate_cyclic_group(5, name='MyZ5', description='custom')
        self.assertEqual(z5.name, 'MyZ5')
        self.assertEqual(z5.description, 'custom')

    def test_elem_name_prefix(self):
        z5 = generate_cyclic_group(5, elem_name='g')
        self.assertEqual(z5.elements, ('g0', 'g1', 'g2', 'g3', 'g4'))

    def test_zfill(self):
        z11 = generate_cyclic_group(11, zfill=True)
        self.assertEqual(z11.elements[:3], ('00', '01', '02'))

    def test_operation_is_addition_mod_n(self):
        z5 = generate_cyclic_group(5)
        self.assertEqual(z5.op('3', '4'), '2')


class TestGenerateSymmetricGroup(TestCase):

    def test_order_and_name(self):
        s3 = generate_symmetric_group(3)
        self.assertEqual(s3.order, 6)
        self.assertEqual(s3.name, 'S3')
        self.assertIsInstance(s3, Group)

    def test_alternating_group(self):
        a3 = generate_symmetric_group(3, alternating=True)
        self.assertEqual(a3.order, 3)
        self.assertEqual(a3.name, 'A3')

    def test_cyclic_form_false_uses_tuple_names(self):
        s3 = generate_symmetric_group(3, cyclic_form=False)
        self.assertIn('(0, 1, 2)', s3.elements)

    def test_custom_name(self):
        s3 = generate_symmetric_group(3, name='Sym3')
        self.assertEqual(s3.name, 'Sym3')


class TestGeneratePowersetGroup(TestCase):

    def test_order_is_power_of_2(self):
        ps3 = generate_powerset_group(3)
        self.assertEqual(ps3.order, 8)

    def test_empty_set_element_name(self):
        ps3 = generate_powerset_group(3)
        self.assertIn('{}', ps3.elements)

    def test_is_a_group(self):
        ps3 = generate_powerset_group(3)
        self.assertIsInstance(ps3, Group)
        self.assertTrue(ps3.is_abelian())


class TestGenerateCommutativeMonoid(TestCase):

    def test_basic(self):
        m4 = generate_commutative_monoid(4)
        self.assertIsInstance(m4, Monoid)
        self.assertEqual(m4.elements, ('a0', 'a1', 'a2', 'a3'))
        self.assertEqual(m4.name, 'M4')

    def test_operation_is_mult_mod_n(self):
        m4 = generate_commutative_monoid(4)
        self.assertEqual(m4.op('a2', 'a3'), 'a2')  # 2*3 mod 4 = 2

    def test_custom_elem_name(self):
        m4 = generate_commutative_monoid(4, elem_name='x')
        self.assertEqual(m4.elements, ('x0', 'x1', 'x2', 'x3'))


class TestGenerateRelativePrimesGroup(TestCase):

    def test_relative_primes_of_8(self):
        rp = generate_relative_primes_group(8)
        self.assertEqual(set(rp.elements), {'1', '3', '5', '7'})
        self.assertEqual(rp.order, 4)
        self.assertIsInstance(rp, Group)

    def test_operation_is_mult_mod_n(self):
        rp = generate_relative_primes_group(8)
        self.assertEqual(rp.op('3', '5'), '7')  # 15 mod 8 = 7


class TestGeneratePowersetRing(TestCase):

    def test_order_and_type(self):
        psr = generate_powerset_ring(2)
        self.assertEqual(psr.order, 4)
        self.assertIsInstance(psr, Ring)

    def test_addition_is_symmetric_difference_multiplication_is_intersection(self):
        psr = generate_powerset_ring(2)
        # {0} + {0,1} = {1} (symmetric difference)
        self.assertEqual(psr.add('{0}', '{0, 1}'), '{1}')
        # {0} * {0,1} = {0} (intersection)
        self.assertEqual(psr.mult('{0}', '{0, 1}'), '{0}')

    def test_invalid_n_raises_on_default_description(self):
        with self.assertRaises(ValueError):
            generate_powerset_ring(0)


class TestGenerateAlgebraModN(TestCase):

    def test_prime_yields_field(self):
        f7 = generate_algebra_mod_n(7)
        self.assertIsInstance(f7, Field)
        self.assertEqual(f7.name, 'F7')

    def test_composite_yields_ring(self):
        r8 = generate_algebra_mod_n(8)
        self.assertIsInstance(r8, Ring)
        self.assertNotIsInstance(r8, Field)
        self.assertEqual(r8.name, 'R8')

    def test_operations(self):
        r8 = generate_algebra_mod_n(8)
        self.assertEqual(r8.add('5', '6'), '3')
        self.assertEqual(r8.mult('5', '6'), '6')


class TestGenerateNxNMatrixAlgebra(TestCase):

    def test_order_and_type(self):
        f3 = generate_algebra_mod_n(3)
        alg, elem_dict, rev_dict = generate_nxn_matrix_algebra(f3)
        self.assertIsInstance(alg, Ring)
        self.assertEqual(alg.order, 3 ** 4)  # 2x2 matrices over a 3-element ring

    def test_element_dict_values_are_abstract_matrices(self):
        f3 = generate_algebra_mod_n(3)
        _, elem_dict, _ = generate_nxn_matrix_algebra(f3)
        for name, mat in elem_dict.items():
            self.assertIsInstance(mat, AbstractMatrix)

    def test_rev_dict_round_trips(self):
        f3 = generate_algebra_mod_n(3)
        _, elem_dict, rev_dict = generate_nxn_matrix_algebra(f3)
        for name, mat in elem_dict.items():
            self.assertEqual(rev_dict[mat.to_tuple()], name)

    def test_custom_prefix(self):
        f3 = generate_algebra_mod_n(3)
        _, elem_dict, _ = generate_nxn_matrix_algebra(f3, element_name_prefix='m')
        self.assertTrue(all(name.startswith('m') for name in elem_dict))


class TestGenerateDihedralGroup(TestCase):

    def test_order_and_name(self):
        d4, elem_dict, iterations = generate_dihedral_group(4)
        self.assertEqual(d4.order, 8)
        self.assertEqual(d4.name, 'D8')
        self.assertIsInstance(d4, Group)

    def test_generator_names_present(self):
        d4, elem_dict, iterations = generate_dihedral_group(4)
        self.assertIn('e', elem_dict)
        self.assertIn('r', elem_dict)
        self.assertIn('f', elem_dict)

    def test_iterations_is_positive(self):
        d4, _, iterations = generate_dihedral_group(4)
        self.assertGreaterEqual(iterations, 1)

    def test_not_abelian_for_n_ge_3(self):
        d4, _, _ = generate_dihedral_group(4)
        self.assertFalse(d4.is_abelian())


class TestGenerateAlgebraFromElementDict(TestCase):

    def test_closure_of_sign_group(self):
        # {1, -1} under multiplication is already closed.
        alg, elem_dict, counter = generate_algebra_from_element_dict(
            {'e': 1, 'a': -1}, bin_op=lambda x, y: x * y)
        self.assertEqual(alg.order, 2)
        self.assertEqual(counter, 1)

    def test_expands_generators_to_closure(self):
        # Generate Z3 (as integers mod 3, represented via a custom bin_op) starting
        # from a single generator '1' under addition mod 3.
        alg, elem_dict, counter = generate_algebra_from_element_dict(
            {'a': 1},
            bin_op=lambda x, y: (x + y) % 3,
            max_iter=10,
        )
        self.assertEqual(alg.order, 3)
        self.assertGreaterEqual(counter, 1)

    def test_custom_elem_key_function(self):
        alg, elem_dict, counter = generate_algebra_from_element_dict(
            {'e': 1, 'a': -1},
            bin_op=lambda x, y: x * y,
            make_elem_key=lambda k1, k2, prod: f"{k1}*{k2}",
        )
        self.assertEqual(alg.order, 2)
