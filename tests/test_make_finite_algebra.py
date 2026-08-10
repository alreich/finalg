"""
Unit tests for finalg.make_finite_algebra.make_finite_algebra, which dispatches
to the correct FiniteAlgebra or FiniteCompositeAlgebra subclass based on its
inputs' mathematical properties.
"""

import json
import os
import tempfile
from unittest import TestCase

from finalg import make_finite_algebra, generate_algebra_mod_n, generate_cyclic_group
from finalg.magma import Magma
from finalg.quasigroup_and_loop import Quasigroup, Loop
from finalg.semigroup import Semigroup
from finalg.monoid import Monoid
from finalg.group import Group
from finalg.ring import Ring
from finalg.field import Field
from finalg.module import Module, module_sv_mult
from finalg.vector_space import VectorSpace


class TestFourArgDispatch(TestCase):
    """4 args (name, description, elements, table) selects among
    Magma/Quasigroup/Loop/Semigroup/Monoid/Group based on table properties."""

    def test_dispatches_to_magma(self):
        # RPS: not associative, no cancellation -> plain Magma
        alg = make_finite_algebra('RPS', 'rock paper scissors', ['r', 'p', 's'],
                                  [[0, 1, 0], [1, 1, 2], [0, 2, 2]])
        self.assertIsInstance(alg, Magma)
        self.assertEqual(type(alg).__name__, 'Magma')

    def test_dispatches_to_quasigroup(self):
        tbl9 = [[0, 4, 8, 2, 3, 9, 6, 7, 1, 5],
                [3, 6, 2, 8, 7, 1, 9, 5, 0, 4],
                [8, 9, 3, 1, 0, 6, 4, 2, 5, 7],
                [1, 7, 6, 5, 4, 8, 0, 3, 2, 9],
                [2, 1, 9, 0, 6, 7, 5, 8, 4, 3],
                [5, 2, 7, 4, 9, 3, 1, 0, 8, 6],
                [4, 3, 0, 6, 1, 5, 2, 9, 7, 8],
                [9, 8, 5, 7, 2, 0, 3, 4, 6, 1],
                [7, 0, 1, 9, 5, 4, 8, 6, 3, 2],
                [6, 5, 4, 3, 8, 2, 7, 1, 9, 0]]
        alg = make_finite_algebra('QG9', 'a quasigroup', [str(i) for i in range(10)], tbl9)
        self.assertIsInstance(alg, Quasigroup)
        self.assertNotIsInstance(alg, Loop)

    def test_dispatches_to_loop(self):
        loop_tbl = [[0, 1, 2, 3, 4, 5, 6],
                   [1, 2, 0, 5, 6, 4, 3],
                   [2, 0, 1, 6, 5, 3, 4],
                   [3, 6, 5, 4, 0, 1, 2],
                   [4, 5, 6, 0, 3, 2, 1],
                   [5, 3, 4, 2, 1, 6, 0],
                   [6, 4, 3, 1, 2, 0, 5]]
        alg = make_finite_algebra('L7', 'a loop', [str(i) for i in range(7)], loop_tbl)
        self.assertIsInstance(alg, Loop)

    def test_dispatches_to_semigroup(self):
        tbl = [[0, 3, 0, 3, 0, 3],
              [1, 4, 1, 4, 1, 4],
              [2, 5, 2, 5, 2, 5],
              [3, 0, 3, 0, 3, 0],
              [4, 1, 4, 1, 4, 1],
              [5, 2, 5, 2, 5, 2]]
        alg = make_finite_algebra('SG', 'a semigroup', ['a', 'b', 'c', 'd', 'e', 'f'], tbl)
        self.assertIsInstance(alg, Semigroup)
        self.assertNotIsInstance(alg, Monoid)

    def test_dispatches_to_monoid(self):
        tbl = [[(a * b) % 4 for b in range(4)] for a in range(4)]
        alg = make_finite_algebra('M4', 'mult mod 4', ['a0', 'a1', 'a2', 'a3'], tbl)
        self.assertIsInstance(alg, Monoid)
        self.assertNotIsInstance(alg, Group)

    def test_dispatches_to_group(self):
        tbl = [[0, 1, 2, 3], [1, 2, 3, 0], [2, 3, 0, 1], [3, 0, 1, 2]]
        alg = make_finite_algebra('Z4', 'add mod 4', ['0', '1', '2', '3'], tbl)
        self.assertIsInstance(alg, Group)


class TestFiveArgDispatch(TestCase):
    """5 args are either (Ring/Field tables) or (scalars, vectors, sv_mult)."""

    def test_dispatches_to_ring_when_not_a_field(self):
        elements = ['0', '1', '2', '3']
        add_tbl = [[0, 1, 2, 3], [1, 0, 3, 2], [2, 3, 0, 1], [3, 2, 1, 0]]
        mult_tbl = [[0, 0, 0, 0], [0, 1, 2, 3], [0, 2, 0, 2], [0, 3, 2, 1]]
        alg = make_finite_algebra('R4', 'ring mod 4', elements, add_tbl, mult_tbl)
        self.assertIsInstance(alg, Ring)
        self.assertNotIsInstance(alg, Field)

    def test_dispatches_to_field_for_prime_modulus(self):
        elements = ['0', '1', '2', '3', '4']
        add_tbl = [[(a + b) % 5 for b in range(5)] for a in range(5)]
        mult_tbl = [[(a * b) % 5 for b in range(5)] for a in range(5)]
        alg = make_finite_algebra('F5', 'field mod 5', elements, add_tbl, mult_tbl)
        self.assertIsInstance(alg, Field)

    def test_dispatches_to_vector_space_for_field_and_group(self):
        f5 = generate_algebra_mod_n(5)
        group = f5 ** 2
        vs = make_finite_algebra('VS', 'a vector space', f5, group, module_sv_mult(f5))
        self.assertIsInstance(vs, VectorSpace)

    def test_dispatches_to_module_for_ring_and_group(self):
        r4 = generate_algebra_mod_n(4)
        group = r4 ** 2
        mod = make_finite_algebra('MOD', 'a module', r4, group, module_sv_mult(r4))
        self.assertIsInstance(mod, Module)
        self.assertNotIsInstance(mod, VectorSpace)


class TestSixArgDispatch(TestCase):
    """6 args build a Ring/Field with a conjugate mapping, as with Cayley-Dickson algebras."""

    def test_rebuild_cayley_dickson_ring_with_conjugates(self):
        f3 = generate_algebra_mod_n(3)
        zi3 = f3.make_cayley_dickson_algebra()
        rebuilt = make_finite_algebra(
            zi3.name, zi3.description, list(zi3.elements),
            zi3.add_table.tolist(), zi3.mult_table.tolist(), zi3.conjugates())
        self.assertEqual(rebuilt, zi3)
        self.assertIsNotNone(rebuilt.conjugates())


class TestSingleArgDispatch(TestCase):
    """1 arg is either a JSON filepath or a dict, and must round-trip an algebra."""

    def test_from_dict(self):
        r4 = generate_algebra_mod_n(4)
        rebuilt = make_finite_algebra(r4.to_dict())
        self.assertEqual(rebuilt, r4)

    def test_from_json_file(self):
        r4 = generate_algebra_mod_n(4)
        with tempfile.TemporaryDirectory() as tmpdir:
            path = os.path.join(tmpdir, 'r4.json')
            r4.dump(path)
            rebuilt = make_finite_algebra(path)
        self.assertEqual(rebuilt, r4)

    def test_from_invalid_single_arg_type_raises(self):
        with self.assertRaises(ValueError):
            make_finite_algebra(12345)


class TestErrorPaths(TestCase):

    def test_wrong_number_of_args_raises(self):
        with self.assertRaises(ValueError):
            make_finite_algebra('a', 'b', 'c')

    def test_too_many_args_raises(self):
        with self.assertRaises(ValueError):
            make_finite_algebra(1, 2, 3, 4, 5, 6, 7)

    def test_duplicate_element_names_raise(self):
        with self.assertRaises(ValueError):
            make_finite_algebra('n', 'd', ['a', 'a'], [[0, 0], [0, 0]])

    def test_non_string_element_names_raise(self):
        with self.assertRaises(ValueError):
            make_finite_algebra('n', 'd', [1, 2], [[0, 0], [0, 0]])
