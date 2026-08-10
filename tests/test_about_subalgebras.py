"""
Unit tests for finalg.about_subalgebras: functions for analyzing and printing
partitions of a finite algebra's proper subalgebras, grouped by isomorphism.
"""

import io
import contextlib
from unittest import TestCase

from finalg import generate_symmetric_group, generate_cyclic_group
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
        # NOTE: about_isomorphic_partition accesses part[0] unconditionally near the
        # top of the function (before the `size == 0` branch that would raise
        # ValueError is ever reached), so an empty partition actually raises
        # IndexError rather than the ValueError the final `else` clause suggests.
        # This test documents that real, currently-existing behavior.
        with self.assertRaises(IndexError):
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
