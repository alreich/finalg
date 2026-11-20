"""
Unit tests for finalg.group_construction: automorphisms, cyclic extensions, and groups_of_order.
"""

from unittest import TestCase

import numpy as np

from finalg import generate_cyclic_group
from finalg.group import Group
from finalg.group_construction import (automorphisms, inner_automorphisms, automorphism_group,
                                       valid_extension_data, reduce_by_conjugacy, build_cyclic_extension,
                                       cyclic_extensions, relabel_group, groups_of_order)
from finalg.group_catalog import small_groups

# The number of groups of order 1, 2, ..., 31 (OEIS A000001)
GROUP_COUNTS = (1, 1, 1, 2, 1, 2, 1, 5, 2, 2, 1, 5, 1, 2, 1, 14, 1, 5, 1, 5,
                2, 2, 1, 15, 2, 2, 5, 4, 1, 4, 1)


def small(order, name):
    """A group from the catalogue of small groups, by name."""
    return next(g for g in small_groups(order) if g.name == name)


def catalogue_name(group):
    """The name of the catalogue group that the given group is isomorphic to, found with fast_isomorphic."""
    matches = [s.name for s in small_groups(group.order) if group.fast_isomorphic(s)]
    assert len(matches) == 1, f"{group.name} matches {matches}"
    return matches[0]


class TestAutomorphisms(TestCase):

    def test_automorphism_counts(self):
        # (group, |Aut|, |Inn|)
        cases = [(generate_cyclic_group(1), 1, 1),
                 (small(5, "Z5"), 4, 1), (small(8, "Z8"), 4, 1), (small(4, "V4"), 6, 1),
                 (small(6, "D3"), 6, 6), (small(8, "D4"), 8, 4), (small(8, "Q8"), 24, 4),
                 (small(8, "Z2_x_Z2_x_Z2"), 168, 1)]
        for group, num_auts, num_inner in cases:
            with self.subTest(group=group.name):
                self.assertEqual(len(automorphisms(group)), num_auts)
                self.assertEqual(len(inner_automorphisms(group)), num_inner)

    def test_automorphisms_are_automorphisms(self):
        for group in [small(8, "Q8"), small(6, "D3"), small(8, "Z4_x_Z2"), small(12, "A4")]:
            T = np.array(group.table.table)
            n = group.order
            auts = automorphisms(group)
            with self.subTest(group=group.name):
                self.assertEqual(len(set(auts)), len(auts))
                for f in auts:
                    F = np.array(f)
                    self.assertEqual(sorted(f), list(range(n)))
                    self.assertTrue(np.array_equal(F[T], T[np.ix_(F, F)]))

    def test_inner_automorphisms_are_automorphisms(self):
        group = small(8, "D4")
        self.assertTrue(set(inner_automorphisms(group)) <= set(automorphisms(group)))

    def test_automorphism_group_of_cyclic_group_is_unit_group(self):
        aut = automorphism_group(small(8, "Z8"))
        self.assertIsInstance(aut, Group)
        self.assertEqual(aut.order, 4)
        self.assertTrue(aut.is_commutative())
        self.assertEqual(catalogue_name(aut), "V4")      # (Z/8)* is Z2 x Z2

    def test_automorphism_group_of_q8_is_s4(self):
        self.assertEqual(catalogue_name(automorphism_group(small(8, "Q8"))), "S4")

    def test_automorphism_group_names_and_maps(self):
        group = small(6, "D3")
        aut, maps = automorphism_group(group, return_maps=True)
        self.assertEqual(aut.identity, "f0")
        self.assertEqual(maps[0], tuple(range(group.order)))
        self.assertEqual(aut.elements, tuple(f"f{i}" for i in range(len(maps))))
        # f_i * f_j is the composition f_i o f_j
        for i in range(len(maps)):
            for j in range(len(maps)):
                composed = tuple(maps[i][maps[j][a]] for a in range(group.order))
                self.assertEqual(maps[int(aut.op(f"f{i}", f"f{j}")[1:])], composed)

    def test_automorphism_group_max_order(self):
        with self.assertRaises(ValueError):
            automorphism_group(small(16, "Z2_x_Z2_x_Z2_x_Z2"), max_order=100)


class TestCyclicExtensions(TestCase):

    def test_extensions_of_z4_by_z2(self):
        z4 = small(4, "Z4")
        names = sorted(catalogue_name(g) for g in cyclic_extensions(z4, 2))
        self.assertEqual(names, sorted(["Z8", "Z4_x_Z2", "D4", "Q8"]))

    def test_extension_data_for_nonsplit_cases(self):
        z4 = small(4, "Z4")
        by_name = {catalogue_name(g): (alpha, c) for g, alpha, c in cyclic_extensions(z4, 2, return_data=True)}
        identity_index = list(z4.elements).index(z4.identity)
        self.assertEqual(by_name["Z4_x_Z2"][1], identity_index)    # direct product: split, c = e
        self.assertEqual(by_name["D4"][1], identity_index)         # semidirect product: split, c = e
        self.assertNotEqual(by_name["Z8"][1], identity_index)      # non-split
        self.assertNotEqual(by_name["Q8"][1], identity_index)      # non-split

    def test_other_extensions(self):
        v4, z3, z5 = small(4, "V4"), small(3, "Z3"), small(5, "Z5")
        self.assertEqual(sorted(catalogue_name(g) for g in cyclic_extensions(v4, 2)),
                         sorted(["Z2_x_Z2_x_Z2", "Z4_x_Z2", "D4"]))
        self.assertEqual(sorted(catalogue_name(g) for g in cyclic_extensions(z3, 2)), sorted(["Z6", "D3"]))
        # m need not be prime
        self.assertEqual(sorted(catalogue_name(g) for g in cyclic_extensions(z3, 4)), sorted(["Z12", "Dic3"]))
        self.assertEqual(sorted(catalogue_name(g) for g in cyclic_extensions(z5, 4)),
                         sorted(["Z20", "F20", "Dic5"]))

    def test_extension_of_trivial_group_is_cyclic(self):
        for m in (2, 3, 7, 10):
            with self.subTest(m=m):
                exts = cyclic_extensions(generate_cyclic_group(1), m)
                self.assertEqual(len(exts), 1)
                self.assertTrue(exts[0].fast_isomorphic(generate_cyclic_group(m)))

    def test_valid_extension_data_satisfies_the_conditions(self):
        n_grp = small(8, "Q8")
        auts = automorphisms(n_grp)
        for m in (2, 3):
            data = valid_extension_data(n_grp, m, auts)
            self.assertTrue(data)
            for alpha, c in data:
                self.assertEqual(alpha[c], c)
                # alpha^m must be conjugation by c
                T = np.array(n_grp.table.table)
                inv = {int(r): int(col) for r, col in zip(*np.where(T == list(n_grp.elements).index(n_grp.identity)))}
                power = np.arange(n_grp.order)
                for _ in range(m):
                    power = np.array(alpha)[power]
                conj = np.array([T[T[c, a], inv[c]] for a in range(n_grp.order)])
                self.assertTrue(np.array_equal(power, conj))

    def test_reduce_by_conjugacy_keeps_at_least_one_and_no_more_than_all(self):
        n_grp = small(8, "D4")
        auts = automorphisms(n_grp)
        data = valid_extension_data(n_grp, 2, auts)
        reps = reduce_by_conjugacy(data, auts)
        self.assertTrue(0 < len(reps) <= len(data))
        self.assertTrue(set(reps) <= set(data))

    def test_build_cyclic_extension(self):
        z4 = small(4, "Z4")
        for g, alpha, c in cyclic_extensions(z4, 2, return_data=True):
            rebuilt = build_cyclic_extension(z4, 2, alpha, c)
            self.assertIsInstance(rebuilt, Group)
            self.assertEqual(rebuilt.order, 8)
            self.assertEqual(rebuilt.elements, g.elements)
            self.assertEqual(rebuilt.table, g.table)

    def test_extension_contains_n_as_normal_subgroup_with_cyclic_quotient(self):
        cases = [(small(4, "Z4"), 2), (small(4, "V4"), 2), (small(3, "Z3"), 4), (small(8, "Q8"), 2),
                 (small(6, "D3"), 3)]
        for n_grp, m in cases:
            delimiter = n_grp.direct_product_delimiter()
            for g in cyclic_extensions(n_grp, m):
                with self.subTest(N=n_grp.name, m=m, G=g.name):
                    sub = g.subalgebra_from_elements([f"{x}{delimiter}0" for x in n_grp.elements],
                                                     name="N_in_G", desc="copy of N")
                    self.assertTrue(g.is_normal(sub))
                    quotient = g / sub
                    self.assertEqual(quotient.order, m)
                    self.assertTrue(quotient.is_cyclic())
                    self.assertTrue(sub.fast_isomorphic(n_grp))


class TestGroupsOfOrder(TestCase):

    def test_counts_match_oeis_a000001(self):
        for n, count in enumerate(GROUP_COUNTS, start=1):
            with self.subTest(n=n):
                self.assertEqual(len(groups_of_order(n)), count)

    def test_groups_match_the_catalogue_of_small_groups(self):
        for n in range(1, 25):
            with self.subTest(n=n):
                names = [catalogue_name(g) for g in groups_of_order(n)]
                self.assertEqual(len(set(names)), len(names))
                self.assertEqual(set(names), {g.name for g in small_groups(n)})

    def test_groups_are_groups_with_standard_names(self):
        for n in (1, 6, 12):
            for k, g in enumerate(groups_of_order(n), start=1):
                self.assertIsInstance(g, Group)
                self.assertEqual(g.order, n)
                self.assertEqual(g.name, f"G{n}_{k}")
                self.assertEqual(g.elements, tuple(str(i) for i in range(n)))
                self.assertEqual(g.identity, "0")

    def test_the_groups_are_not_isomorphic_to_each_other(self):
        gs = groups_of_order(16)
        for i, a in enumerate(gs):
            for b in gs[i + 1:]:
                self.assertFalse(a.fast_isomorphic(b))

    def test_custom_cache_is_used_and_filled(self):
        cache = {}
        groups_of_order(12, cache)
        self.assertTrue({1, 2, 3, 4, 6, 12} <= set(cache))
        self.assertEqual(len(cache[12]), 5)

    def test_returned_list_is_a_copy(self):
        gs = groups_of_order(6)
        gs.clear()
        self.assertEqual(len(groups_of_order(6)), 2)

    def test_warns_for_orders_where_nonsolvable_groups_may_exist(self):
        with self.assertWarns(UserWarning):
            gs = groups_of_order(61)
        self.assertEqual(len(gs), 1)

    def test_no_warning_below_60(self):
        import warnings
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            groups_of_order(9)

    def test_invalid_orders(self):
        for bad in (0, -3, 2.5, "8", None, True):
            with self.subTest(n=bad):
                with self.assertRaises(ValueError):
                    groups_of_order(bad)


class TestRelabelGroup(TestCase):

    def test_relabel_group(self):
        q8 = small(8, "Q8")
        new = relabel_group(q8)
        self.assertEqual(new.name, "Q8")
        self.assertEqual(new.elements, tuple(str(i) for i in range(8)))
        self.assertEqual(new.identity, "0")
        self.assertTrue(new.fast_isomorphic(q8))
        self.assertEqual(relabel_group(q8, name="Quat").name, "Quat")
