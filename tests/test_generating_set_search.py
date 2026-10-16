"""Tests for the search for a small generating set (Magma._smallest_generating_set), which
fast_isomorphic depends on."""

import time
from unittest import TestCase

from finalg import (generate_cyclic_group, generate_algebra_mod_n, make_finite_algebra,
                    small_groups, identify_small_group)


def spans(alg, gens):
    """Brute-force closure of `gens` under op (plus nothing else)."""
    found = set(gens)
    changed = True
    while changed:
        changed = False
        for a in list(found):
            for b in list(found):
                c = alg.op(a, b)
                if c not in found:
                    found.add(c)
                    changed = True
    return found


class TestGeneratingSetSearch(TestCase):

    def test_every_small_group_is_generated_by_the_set(self):
        for grps in small_groups().values():
            for g in grps:
                gens = g._smallest_generating_set()
                closure = spans(g, gens) | {g.identity}
                self.assertEqual(closure, set(g.elements), g.name)

    def test_cyclic_groups_need_one_generator(self):
        for n in (1, 2, 7, 12, 32, 40):
            self.assertEqual(len(generate_cyclic_group(n)._smallest_generating_set()), 1, n)

    def test_group_of_order_one(self):
        g = generate_cyclic_group(1)
        self.assertEqual(g._smallest_generating_set(), [g.elements[0]])

    def test_klein_four_needs_two(self):
        self.assertEqual(len(small_groups(4)[1]._smallest_generating_set()), 2)

    def test_elementary_abelian_group_of_order_32_needs_five(self):
        z2 = generate_cyclic_group(2)
        g = z2 * z2 * z2 * z2 * z2
        self.assertEqual(len(g._smallest_generating_set()), 5)

    def test_ring_is_generated_additively(self):
        r = generate_algebra_mod_n(6)
        gens = r._smallest_generating_set()
        self.assertEqual(len(gens), 1)
        self.assertEqual(spans(r, gens) | {r.identity}, set(r.elements))

    def test_returned_list_is_a_copy(self):
        g = small_groups(8)[4]
        first = g._smallest_generating_set()
        first.append('junk')
        self.assertNotIn('junk', g._smallest_generating_set())

    def test_nonassociative_magma(self):
        # x*y = (x - y) mod 3 is not associative; the search must still produce a generating set
        elements = ['0', '1', '2']
        table = [[str((i - j) % 3) for j in range(3)] for i in range(3)]
        m = make_finite_algebra('Sub3', 'Subtraction mod 3', elements, [[int(c) for c in row] for row in table])
        gens = m._smallest_generating_set()
        self.assertEqual(spans(m, gens), set(m.elements))

    def test_order_32_isomorphism_is_fast(self):
        # Used to take about 40 seconds for some pairs
        for g in small_groups(32):
            relabeled = g.reorder_elements([g.elements[0]] + list(reversed(g.elements[1:])))
            start = time.time()
            self.assertIs(identify_small_group(relabeled), g, g.name)
            self.assertLess(time.time() - start, 30, g.name)
