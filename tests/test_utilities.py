"""
Unit tests for finalg.utilities: assorted helper functions.
"""

from unittest import TestCase
import numpy as np
import networkx as nx

from finalg.utilities import (
    np_arr_to_tuple,
    delete_row_col,
    get_duplicates,
    all_strings,
    replace_item,
    yes_or_no,
    symm_diff_of_two_lists_of_lists,
    same_lists_of_lists,
    powerset,
    cayley_graph_to_json,
    make_table_from_xml,
    compress_runs,
)


class TestNpArrToTuple(TestCase):

    def test_conversion(self):
        arr = np.array([[1, 2], [3, 4]])
        self.assertEqual(np_arr_to_tuple(arr), ((1, 2), (3, 4)))

    def test_result_is_hashable(self):
        arr = np.array([[1, 2], [3, 4]])
        d = {np_arr_to_tuple(arr): 'ok'}
        self.assertEqual(d[((1, 2), (3, 4))], 'ok')


class TestDeleteRowCol(TestCase):

    def test_deletes_specified_row_and_col(self):
        arr = np.array([[1, 2], [3, 4]])
        result = delete_row_col(arr, 0, 1)
        self.assertEqual(result.tolist(), [[3]])

    def test_does_not_modify_original(self):
        arr = np.array([[1, 2], [3, 4]])
        delete_row_col(arr, 0, 1)
        self.assertEqual(arr.tolist(), [[1, 2], [3, 4]])

    def test_larger_array(self):
        arr = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        result = delete_row_col(arr, 1, 2)
        self.assertEqual(result.tolist(), [[1, 2], [7, 8]])


class TestGetDuplicates(TestCase):

    def test_finds_duplicates(self):
        self.assertEqual(sorted(get_duplicates([1, 2, 2, 3, 3, 3])), [2, 3])

    def test_no_duplicates_returns_empty(self):
        self.assertEqual(get_duplicates([1, 2, 3]), [])

    def test_works_with_strings(self):
        self.assertEqual(get_duplicates(['a', 'b', 'a']), ['a'])


class TestAllStrings(TestCase):

    def test_all_strings_true(self):
        self.assertTrue(all_strings(['a', 'b', 'c']))

    def test_mixed_types_false(self):
        self.assertFalse(all_strings(['a', 1]))

    def test_empty_iterable_true(self):
        self.assertTrue(all_strings([]))


class TestReplaceItem(TestCase):

    def test_replace_in_list(self):
        self.assertEqual(replace_item([1, 2, 3, 2], 2, 99), [1, 99, 3, 99])

    def test_replace_in_tuple_preserves_type(self):
        result = replace_item((1, 2, 3, 2), 2, 99)
        self.assertEqual(result, (1, 99, 3, 99))
        self.assertIsInstance(result, tuple)

    def test_no_match_leaves_unchanged(self):
        self.assertEqual(replace_item([1, 2, 3], 99, 0), [1, 2, 3])


class TestYesOrNo(TestCase):

    def test_true_is_yes(self):
        self.assertEqual(yes_or_no(True), "Yes")

    def test_false_is_no(self):
        self.assertEqual(yes_or_no(False), "No")


class TestSymmDiffOfTwoListsOfLists(TestCase):

    def test_symmetric_difference(self):
        result = symm_diff_of_two_lists_of_lists([[1, 2], [3, 4]], [[3, 4], [5, 6]])
        self.assertEqual(result, {(1, 2), (5, 6)})

    def test_identical_lists_gives_empty_set(self):
        result = symm_diff_of_two_lists_of_lists([[1, 2], [3, 4]], [[1, 2], [3, 4]])
        self.assertEqual(result, set())


class TestSameListsOfLists(TestCase):

    def test_same_regardless_of_order(self):
        self.assertTrue(same_lists_of_lists([[1, 2], [3, 4]], [[3, 4], [1, 2]]))

    def test_different_lists_are_not_same(self):
        self.assertFalse(same_lists_of_lists([[1, 2]], [[3, 4]]))


class TestPowerset(TestCase):

    def test_powerset_of_three_elements(self):
        result = list(powerset([1, 2, 3]))
        expected = [(), (1,), (2,), (3,), (1, 2), (1, 3), (2, 3), (1, 2, 3)]
        self.assertEqual(result, expected)

    def test_powerset_of_empty_iterable(self):
        self.assertEqual(list(powerset([])), [()])

    def test_powerset_size(self):
        self.assertEqual(len(list(powerset(range(4)))), 16)


class TestCayleyGraphToJson(TestCase):

    def test_returns_node_link_dict(self):
        g = nx.MultiDiGraph()
        g.add_node('a')
        g.add_edge('a', 'a', label='x')
        result = cayley_graph_to_json(g)
        self.assertIsInstance(result, dict)
        self.assertIn('nodes', result)
        self.assertIn('edges', result)
        self.assertTrue(result['directed'])
        self.assertTrue(result['multigraph'])


class TestMakeTableFromXml(TestCase):

    def test_basic_conversion(self):
        result = make_table_from_xml('0 1\n1 0')
        self.assertEqual(result, [[0, 1], [1, 0]])

    def test_larger_table(self):
        result = make_table_from_xml('0 1 2\n1 2 0\n2 0 1')
        self.assertEqual(result, [[0, 1, 2], [1, 2, 0], [2, 0, 1]])


class TestCompressRuns(TestCase):

    def test_repeated_chars_compressed(self):
        self.assertEqual(compress_runs('ffrrrf'), 'f^2r^3f')

    def test_all_repeats(self):
        self.assertEqual(compress_runs('aaabbc'), 'a^3b^2c')

    def test_no_repeats_unchanged(self):
        self.assertEqual(compress_runs('abcd'), 'abcd')

    def test_all_same_char(self):
        self.assertEqual(compress_runs('aaaa'), 'a^4')

    def test_single_char(self):
        self.assertEqual(compress_runs('a'), 'a')

    def test_empty_string(self):
        self.assertEqual(compress_runs(''), '')

    def test_mixed_runs(self):
        self.assertEqual(compress_runs('ffrrrfrrf'), 'f^2r^3fr^2f')
