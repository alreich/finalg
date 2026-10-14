"""
Unit tests for finalg.about_subalgebras: functions for analyzing and printing
partitions of a finite algebra's proper subalgebras, grouped by isomorphism.
"""

import io
import contextlib
from unittest import TestCase

from finalg import (generate_symmetric_group, generate_cyclic_group, generate_dicyclic_group,
                    generate_alternating_group, generate_dihedral_group, generate_algebra_mod_n)
from finalg.about_subalgebras import (
    partition_into_isomorphic_lists,
    about_isomorphic_partition,
    about_isomorphic_partitions,
    are_n,
    add_s,
    about_subalgebras,
    find_isomorphic_subalgebra,
)


class TestAreN(TestCase):

    def test_zero(self):
        self.assertEqual(are_n(0), 'are no')

    def test_one(self):
        self.assertEqual(are_n(1), 'is 1')

    def test_two_or_more(self):
        self.assertEqual(are_n(2), 'are 2')
        self.assertEqual(are_n(5), 'are 5')


class TestAddS(TestCase):

    def test_singular_unchanged(self):
        self.assertEqual(add_s('cat', 1), 'cat')

    def test_plural_gets_s(self):
        self.assertEqual(add_s('cat', 2), 'cats')
        self.assertEqual(add_s('cat', 0), 'cats')


class TestPartitionIntoIsomorphicLists(TestCase):

    def setUp(self):
        self.s3 = generate_symmetric_group(3)
        self.subs = self.s3.proper_subalgebras()

    def test_s3_subalgebras_partition_into_two_classes(self):
        partitions = partition_into_isomorphic_lists(self.subs)
        sizes = sorted(len(p) for p in partitions)
        self.assertEqual(sizes, [1, 3])

    def test_all_subs_in_a_partition_are_mutually_isomorphic(self):
        partitions = partition_into_isomorphic_lists(self.subs)
        for part in partitions:
            for sub in part:
                self.assertTrue(part[0].isomorphic(sub))

    def test_empty_input_returns_empty_list(self):
        self.assertEqual(partition_into_isomorphic_lists([]), [])

    def test_single_subalgebra(self):
        result = partition_into_isomorphic_lists([self.subs[0]])
        self.assertEqual(result, [[self.subs[0]]])


class TestAboutIsomorphicPartition(TestCase):

    def setUp(self):
        self.s3 = generate_symmetric_group(3)
        self.subs = self.s3.proper_subalgebras()
        self.partitions = partition_into_isomorphic_lists(self.subs)

    def test_prints_multi_member_partition(self):
        multi = [p for p in self.partitions if len(p) > 1][0]
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            about_isomorphic_partition(self.s3, multi)
        out = buf.getvalue()
        self.assertIn('Isomorphic', out)
        self.assertIn('order 2', out)

    def test_prints_single_member_partition(self):
        single = [p for p in self.partitions if len(p) == 1][0]
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            about_isomorphic_partition(self.s3, single)
        out = buf.getvalue()
        self.assertIn('order 3', out)
        self.assertIn('Normal', out)  # A3 is normal in S3

    def test_empty_partition_raises(self):
        # about_isomorphic_partition now checks for an empty partition before
        # touching part[0], so it raises the documented ValueError instead of
        # an incidental IndexError.
        with self.assertRaises(ValueError):
            about_isomorphic_partition(self.s3, [])


class TestAboutIsomorphicPartitions(TestCase):

    def test_prints_summary_for_nonempty_partitions(self):
        s3 = generate_symmetric_group(3)
        subs = s3.proper_subalgebras()
        partitions = partition_into_isomorphic_lists(subs)
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            about_isomorphic_partitions(s3, partitions)
        out = buf.getvalue()
        self.assertIn('unique proper subalgebras', out)
        self.assertIn('out of 4 total subalgebras', out)

    def test_prints_message_for_no_subalgebras(self):
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            about_isomorphic_partitions(None, [])
        self.assertIn('no proper subalgebras', buf.getvalue())


class TestAboutSubalgebras(TestCase):

    def test_returns_partitions_and_prints_summary(self):
        s3 = generate_symmetric_group(3)
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            result = about_subalgebras(s3)
        sizes = sorted(len(p) for p in result)
        self.assertEqual(sizes, [1, 3])
        self.assertIn('Subalgebras of S3', buf.getvalue())


class TestFindIsomorphicSubalgebra(TestCase):

    def setUp(self):
        self.s3 = generate_symmetric_group(3)
        self.subs = self.s3.proper_subalgebras()
        self.partitions = partition_into_isomorphic_lists(self.subs)

    def test_finds_matching_order2_subalgebra(self):
        z4 = generate_cyclic_group(4)
        order2_sub = z4.subalgebra_from_elements(['0', '2'])
        result = find_isomorphic_subalgebra(order2_sub, self.partitions)
        self.assertNotEqual(result, False)
        iso_mapping, iso_group = result
        self.assertEqual(iso_group.order, 2)
        self.assertIsInstance(iso_mapping, dict)

    def test_returns_false_when_no_matching_order(self):
        z5 = generate_cyclic_group(5)  # order 5; no subalgebra of S3 has order 5
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            result = find_isomorphic_subalgebra(z5, self.partitions)
        self.assertFalse(result)
        self.assertIn('Not found', buf.getvalue())

    def test_verbose_smoke(self):
        z4 = generate_cyclic_group(4)
        order2_sub = z4.subalgebra_from_elements(['0', '2'])
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            find_isomorphic_subalgebra(order2_sub, self.partitions, verbose=True)
        self.assertIn('Checking:', buf.getvalue())


def summary(alg, **kwargs):
    """Return what about_subalgebras prints for alg."""
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        about_subalgebras(alg, **kwargs)
    return buf.getvalue()


class TestSmallGroupNamesInSummary(TestCase):
    """The summary names the small group that each set of isomorphic subgroups is isomorphic to."""

    def test_s3_subgroups_are_z2_and_z3(self):
        out = summary(generate_symmetric_group(3))
        self.assertIn('(each isomorphic to the small group Z2)', out)
        self.assertIn('(isomorphic to the small group Z3)', out)
        self.assertEqual(out.count('small group'), 2)

    def test_q8_subgroups_are_z2_and_z4(self):
        q8, _, _ = generate_dicyclic_group(2)
        out = summary(q8)
        self.assertIn('(isomorphic to the small group Z2)', out)
        self.assertIn('(each isomorphic to the small group Z4)', out)
        self.assertNotIn('V4', out)

    def test_a4_subgroups_are_z2_z3_and_v4(self):
        out = summary(generate_alternating_group(4))
        self.assertIn('(each isomorphic to the small group Z2)', out)
        self.assertIn('(each isomorphic to the small group Z3)', out)
        self.assertIn('(isomorphic to the small group V4)', out)

    def test_note_comes_directly_after_the_heading_line(self):
        lines = summary(generate_symmetric_group(3)).splitlines()
        for i, line in enumerate(lines):
            if 'small group' in line:
                self.assertIn(' of order ', lines[i - 1])
                self.assertTrue(lines[i + 1].lstrip().startswith('Group:'))

    def test_s4_subgroups_include_d3_d4_and_a4(self):
        out = summary(generate_symmetric_group(4))
        for name in ('Z2', 'Z3', 'Z4', 'V4', 'D3', 'D4', 'A4'):
            self.assertIn(f'the small group {name})', out)
        self.assertEqual(out.count('small group'), 7)

    def test_dihedral_group_of_order_20_subgroups(self):
        out = summary(generate_dihedral_group(10)[0])
        for name in ('Z2', 'V4', 'Z5', 'D5', 'Z10'):
            self.assertIn(f'the small group {name})', out)
        self.assertNotIn('Dic', out)

    def test_identify_small_false_leaves_the_notes_out(self):
        out = summary(generate_symmetric_group(3), identify_small=False)
        self.assertNotIn('small group', out)
        self.assertIn('Subalgebras of S3', out)

    def test_about_isomorphic_partition_respects_identify_small(self):
        s3 = generate_symmetric_group(3)
        part = [p for p in partition_into_isomorphic_lists(s3.proper_subalgebras()) if len(p) == 3][0]
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            about_isomorphic_partition(s3, part, identify_small=False)
        self.assertNotIn('small group', buf.getvalue())
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            about_isomorphic_partition(s3, part)
        self.assertIn('small group Z2', buf.getvalue())

    def test_about_isomorphic_partitions_respects_identify_small(self):
        s3 = generate_symmetric_group(3)
        partitions = partition_into_isomorphic_lists(s3.proper_subalgebras())
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            about_isomorphic_partitions(s3, partitions, identify_small=False)
        self.assertNotIn('small group', buf.getvalue())

    def test_subgroups_of_order_greater_than_20_get_no_note(self):
        # The proper subgroups of Z42 have orders 2, 3, 6, 7, 14, and 21. Order 21 is too big for the
        # catalogue of small groups, so only the others are named.
        out = summary(generate_cyclic_group(42))
        self.assertIn('of order 21', out)
        self.assertIn('of order 14', out)
        self.assertIn('the small group Z14)', out)
        self.assertNotIn('the small group Z21)', out)

    def test_subalgebras_that_are_not_groups_get_no_note(self):
        out = summary(generate_algebra_mod_n(4))
        self.assertNotIn('small group', out)

    def test_return_value_is_unchanged_by_the_note(self):
        s3 = generate_symmetric_group(3)
        with contextlib.redirect_stdout(io.StringIO()):
            with_notes = about_subalgebras(s3)
            without_notes = about_subalgebras(s3, identify_small=False)
        self.assertEqual([[g.name for g in p] for p in with_notes],
                         [[g.name for g in p] for p in without_notes])
