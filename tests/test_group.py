"""
Unit tests for finalg.group.Group.
"""

import io
import contextlib
from unittest import TestCase

from finalg import generate_symmetric_group, generate_cyclic_group, generate_dihedral_group, make_finite_algebra
from finalg.group import Group


def s3():
    return generate_symmetric_group(3)


class TestGroupConstruction(TestCase):

    def test_generate_symmetric_group_is_group(self):
        g = s3()
        self.assertEqual(type(g).__name__, 'Group')
        self.assertEqual(g.order, 6)

    def test_direct_construction_raises_without_inverses(self):
        # A monoid-like table (mult mod 4) is associative & has identity, but
        # 0 and 2 lack inverses, so it can't be a Group.
        tbl = [[(a * b) % 4 for b in range(4)] for a in range(4)]
        with self.assertRaises(ValueError):
            Group('bad', 'no inverses', ['a0', 'a1', 'a2', 'a3'], tbl)

    def test_direct_construction_raises_when_not_associative(self):
        with self.assertRaises(ValueError):
            Group('bad', 'not associative', ['r', 'p', 's'],
                  [[0, 1, 0], [1, 1, 2], [0, 2, 2]])


class TestBasicGroupOps(TestCase):

    def setUp(self):
        self.s3 = s3()
        self.a = self.s3.elements[1]
        self.b = self.s3.elements[2]

    def test_inv(self):
        for elem in self.s3.elements:
            self.assertEqual(self.s3.op(elem, self.s3.inv(elem)), self.s3.identity)

    def test_sub(self):
        # sub(x, y) = x * inv(y)
        result = self.s3.sub(self.a, self.b)
        self.assertEqual(result, self.s3.op(self.a, self.s3.inv(self.b)))

    def test_conjugate(self):
        result = self.s3.conjugate(self.a, self.b)
        expected = self.s3.op(self.b, self.s3.op(self.a, self.s3.inv(self.b)))
        self.assertEqual(result, expected)

    def test_commutator(self):
        result = self.s3.commutator(self.a, self.b)
        expected = self.s3.op(self.s3.op(self.s3.inv(self.a), self.s3.inv(self.b)),
                               self.s3.op(self.a, self.b))
        self.assertEqual(result, expected)

    def test_inverse_mapping(self):
        mapping = self.s3.inverse_mapping()
        self.assertEqual(set(mapping.keys()), set(self.s3.elements))
        for elem, inv in mapping.items():
            self.assertEqual(self.s3.op(elem, inv), self.s3.identity)


class TestCommutatorsAndSolvability(TestCase):

    def setUp(self):
        self.s3 = s3()
        self.z4 = generate_cyclic_group(4)

    def test_commutators_of_abelian_group_is_trivial(self):
        self.assertEqual(self.z4.commutators(), {self.z4.identity})

    def test_commutators_of_s3(self):
        self.assertEqual(self.s3.commutators(),
                          {'(0)(2)', '(0 1 2)', '(0 2 1)'})

    def test_commutator_subalgebra_is_a3(self):
        sub = self.s3.commutator_subalgebra()
        self.assertEqual(set(sub.elements), {'(0)(2)', '(0 1 2)', '(0 2 1)'})
        self.assertEqual(type(sub).__name__, 'Group')

    def test_s3_is_solvable(self):
        self.assertTrue(self.s3.is_solvable())

    def test_abelian_group_is_solvable(self):
        self.assertTrue(self.z4.is_solvable())

    def test_trivial_group_is_solvable(self):
        trivial = make_finite_algebra('Trivial', 'Trivial group', ['e'], [[0]])
        self.assertTrue(trivial.is_solvable())

    def test_a4_is_solvable(self):
        a4 = generate_symmetric_group(4, alternating=True)
        self.assertTrue(a4.is_solvable())

    def test_s5_is_not_solvable(self):
        s5 = generate_symmetric_group(5)
        self.assertFalse(s5.is_solvable())

    def test_a5_is_not_solvable(self):
        # A5 is simple and non-abelian, so its derived series stalls at A5 itself.
        a5 = generate_symmetric_group(5, alternating=True)
        self.assertFalse(a5.is_solvable())


class TestNormalSubgroupsAndQuotients(TestCase):

    def setUp(self):
        self.s3 = s3()

    def test_trivial_subgroup_and_full_group_are_normal(self):
        for sub in self.s3.trivial_subgroups():
            self.assertTrue(self.s3.is_normal(sub))

    def test_a3_is_normal_in_s3(self):
        a3 = self.s3.commutator_subalgebra()
        self.assertTrue(self.s3.is_normal(a3))

    def test_order_2_subgroups_not_normal_in_s3(self):
        order2_subs = [sub for sub in self.s3.proper_subalgebras() if sub.order == 2]
        self.assertEqual(len(order2_subs), 3)
        for sub in order2_subs:
            self.assertFalse(self.s3.is_normal(sub))

    def test_trivial_subgroups(self):
        subs = self.s3.trivial_subgroups()
        elemsets = sorted([sorted(s.elements) for s in subs], key=len)
        expected = sorted([[self.s3.identity], sorted(self.s3.elements)], key=len)
        self.assertEqual(elemsets, expected)

    def test_subgroups_includes_trivial_and_all_proper(self):
        subs = self.s3.subgroups()
        elemsets = sorted([sorted(s.elements) for s in subs])
        proper = sorted([sorted(s.elements) for s in self.s3.proper_subalgebras()])
        trivial = sorted([sorted(s.elements) for s in self.s3.trivial_subgroups()])
        self.assertEqual(elemsets, sorted(proper + trivial))

    def test_unique_proper_subgroups_deduplicates_conjugates(self):
        unique_subs = self.s3.unique_proper_subgroups()
        orders = sorted(s.order for s in unique_subs)
        # 3 conjugate order-2 subgroups collapse to 1 representative; the order-3
        # (normal) subgroup stands alone.
        self.assertEqual(orders, [2, 3])

    def test_quotient_group_by_normal_subgroup(self):
        a3 = self.s3.commutator_subalgebra()
        q = self.s3.quotient_group(a3)
        self.assertEqual(q.order, 2)
        self.assertEqual(type(q).__name__, 'Group')

    def test_truediv_operator_matches_quotient_group(self):
        a3 = self.s3.commutator_subalgebra()
        q1 = self.s3.quotient_group(a3)
        q2 = self.s3 / a3
        self.assertEqual(q1.elements, q2.elements)

    def test_quotient_group_by_non_normal_subgroup_raises(self):
        non_normal = [sub for sub in self.s3.proper_subalgebras() if sub.order == 2][0]
        with self.assertRaises(ValueError):
            self.s3.quotient_group(non_normal)

    def test_quotient_group_of_z4(self):
        z4 = generate_cyclic_group(4)
        sub = z4.subalgebra_from_elements(['0', '2'])
        q = z4.quotient_group(sub)
        self.assertEqual(q.order, 2)
        self.assertTrue(q.is_abelian())


class TestGroupAbout(TestCase):

    def test_about_smoke_and_content(self):
        g = s3()
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            result = g.about()
        # Unlike Magma.about (which returns None implicitly), Group.about
        # explicitly returns str(self).
        self.assertEqual(result, str(g))
        out = buf.getvalue()
        self.assertIn('** Group **', out)
        self.assertIn('Name: S3', out)
        self.assertIn('Inverse', out)

    def test_about_hides_table_when_too_large(self):
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            s3().about(max_size=1)
        self.assertIn('no table is printed', buf.getvalue())


def mult_action(modulus, k):
    """Action of the cyclic group Z_m (string-named integer elements) on the cyclic group Z_modulus,
    where h acts on x by multiplication by k^h, mod modulus."""
    return lambda h, x: str((int(x) * pow(k, int(h), modulus)) % modulus)


class TestSemidirectProduct(TestCase):

    def setUp(self):
        self.z2 = generate_cyclic_group(2)
        self.z3 = generate_cyclic_group(3)
        self.z4 = generate_cyclic_group(4)
        self.z7 = generate_cyclic_group(7)

    def test_z3_by_z2_with_inversion_is_s3(self):
        g = self.z3.semidirect_product(self.z2, mult_action(3, 2))
        self.assertIsInstance(g, Group)
        self.assertEqual(g.order, 6)
        self.assertEqual(g.name, 'Z3_sdp_Z2')
        self.assertEqual(g.elements, ('0:0', '0:1', '1:0', '1:1', '2:0', '2:1'))
        self.assertFalse(g.is_commutative())
        self.assertTrue(g.fast_isomorphic(s3()))

    def test_z4_by_z2_with_inversion_is_d8(self):
        g = self.z4.semidirect_product(self.z2, mult_action(4, 3))
        d8, _, _ = generate_dihedral_group(4)
        self.assertEqual(g.order, 8)
        self.assertTrue(g.fast_isomorphic(d8))

    def test_frobenius_group_of_order_21(self):
        g = self.z7.semidirect_product(self.z3, mult_action(7, 2))
        self.assertEqual(g.order, 21)
        self.assertFalse(g.is_commutative())
        orders = sorted(g.element_order(e) for e in g.elements)
        self.assertEqual(orders, [1] + [3] * 14 + [7] * 6)

    def test_product_rule(self):
        # (n1, h1) * (n2, h2) = (n1 + phi_h1(n2), h1 + h2), here with Z7 and phi_h(x) = 2^h x
        g = self.z7.semidirect_product(self.z3, mult_action(7, 2))
        self.assertEqual(g.op('3:1', '5:2'), '6:0')  # (3 + 2*5 mod 7, 1 + 2 mod 3) = (6, 0)
        self.assertEqual(g.op('0:1', '1:0'), '2:1')

    def test_trivial_action_gives_direct_product(self):
        g = self.z4.semidirect_product(self.z3, lambda h, x: x)
        self.assertEqual(g, self.z4 * self.z3)

    def test_n_is_normal_subgroup_with_h_as_complement(self):
        g = self.z7.semidirect_product(self.z3, mult_action(7, 2))
        n_sub = g.subalgebra_from_elements([f"{n}:0" for n in self.z7.elements], name='N')
        h_sub = g.subalgebra_from_elements([f"0:{h}" for h in self.z3.elements], name='H')
        self.assertTrue(g.is_normal(n_sub))
        self.assertFalse(g.is_normal(h_sub))
        self.assertEqual(n_sub.order * h_sub.order, g.order)
        self.assertEqual(set(n_sub.elements) & set(h_sub.elements), {'0:0'})

    def test_action_as_dict_of_dicts_or_dict_of_functions(self):
        expected = self.z3.semidirect_product(self.z2, mult_action(3, 2))
        by_dicts = self.z3.semidirect_product(
            self.z2, {'0': {'0': '0', '1': '1', '2': '2'},
                      '1': {'0': '0', '1': '2', '2': '1'}})
        by_funcs = self.z3.semidirect_product(
            self.z2, {'0': lambda x: x, '1': lambda x: str((-int(x)) % 3)})
        self.assertEqual(by_dicts, expected)
        self.assertEqual(by_funcs, expected)

    def test_custom_name_description_and_delimiter(self):
        z3 = generate_cyclic_group(3)
        z3.direct_product_delimiter('|')
        g = z3.semidirect_product(self.z2, mult_action(3, 2), name='MyS3', description='custom')
        self.assertEqual(g.name, 'MyS3')
        self.assertEqual(g.description, 'custom')
        self.assertEqual(g.elements[:3], ('0|0', '0|1', '1|0'))

    def test_action_that_is_not_one_to_one_raises(self):
        with self.assertRaises(ValueError):
            self.z4.semidirect_product(
                self.z2, lambda h, x: str((2 * int(x)) % 4) if h == '1' else x)

    def test_action_that_is_not_a_homomorphism_of_n_raises(self):
        swap_2_and_3 = {'0': '0', '1': '1', '2': '3', '3': '2'}
        with self.assertRaises(ValueError):
            self.z4.semidirect_product(
                self.z2, lambda h, x: swap_2_and_3[x] if h == '1' else x)

    def test_action_that_is_not_a_homomorphism_from_h_raises(self):
        # Inversion for both 1 and 2 in Z3 would need phi_2 = phi_1 o phi_1 = identity.
        with self.assertRaises(ValueError):
            self.z3.semidirect_product(
                self.z3, lambda h, x: x if h == '0' else str((-int(x)) % 3))

    def test_action_missing_an_element_of_h_raises(self):
        with self.assertRaises(ValueError):
            self.z3.semidirect_product(self.z2, {'0': lambda x: x})

    def test_action_mapping_outside_n_raises(self):
        with self.assertRaises(ValueError):
            self.z3.semidirect_product(self.z2, lambda h, x: '9')

    def test_other_must_be_a_group(self):
        with self.assertRaises(TypeError):
            self.z3.semidirect_product('not a group', lambda h, x: x)
