"""
Unit tests for finalg.utilities: assorted helper functions.
"""

import json
import os
import tempfile
from unittest import TestCase
import numpy as np
import networkx as nx

import finalg.small_algebras as sgt
from finalg import make_finite_algebra, generate_cyclic_group
from finalg.group import Group

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
    mult_action,
    action_from_generators,
    pretty_json,
    dump_pretty,
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
        # noinspection SpellCheckingInspection
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


def as_json_data(d):
    """The dictionary as it looks after a trip through JSON, where tuples have become lists."""
    return json.loads(json.dumps(d))


def table_and_elements(alg):
    """What defines an algebra, apart from its name and description."""
    d = alg.to_dict()
    return d['elements'], d['table']


class TestPrettyJson(TestCase):

    def test_the_text_is_valid_json_that_gives_back_the_dictionary(self):
        for name in ('Z1', 'Q8', 'SL2_3', 'SG32_8'):
            d = sgt.load_group(name).to_dict()
            self.assertEqual(json.loads(pretty_json(d)), as_json_data(d))

    def test_each_row_of_the_cayley_table_is_on_its_own_line(self):
        g = sgt.load_group('D4')
        lines = pretty_json(g.to_dict()).splitlines()
        rows = [ln for ln in lines if ln.startswith('    [')]
        self.assertEqual(len(rows), g.order)
        self.assertEqual(lines[0], '{')
        self.assertEqual(lines[-1], '}')

    def test_each_key_is_on_its_own_line(self):
        lines = pretty_json(sgt.load_group('Q8').to_dict()).splitlines()
        for key in ('name', 'description', 'elements', 'table'):
            self.assertEqual(sum(ln.startswith(f'  "{key}":') for ln in lines), 1, key)

    def test_long_element_lists_are_wrapped(self):
        text = pretty_json(sgt.load_group('SG32_8').to_dict(), width=100)
        lines = text.splitlines()
        start = next(i for i, ln in enumerate(lines) if ln.startswith('  "elements":'))
        end = next(i for i in range(start, len(lines)) if lines[i].strip() in (']', '],'))
        self.assertGreater(end - start, 2)
        for ln in lines[start:end + 1]:
            self.assertLessEqual(len(ln), 100)

    def test_a_ring_keeps_both_tables(self):
        # Pretty-printing is not just for groups: any finalg algebra with a second table works
        ring = make_finite_algebra("R", "Z3 as a ring", ["0", "1", "2"],
                                   [[(a + b) % 3 for b in range(3)] for a in range(3)],
                                   [[(a * b) % 3 for b in range(3)] for a in range(3)])
        d = ring.to_dict()
        self.assertIn('table2', d)
        text = pretty_json(d)
        self.assertEqual(json.loads(text), as_json_data(d))
        self.assertEqual(sum(ln.startswith('    [') for ln in text.splitlines()), 6)

    def test_the_files_made_by_finalg_code_are_in_the_pretty_format(self):
        for r in sgt.groups_index():
            if r['source'] != 'finalg':
                continue
            text = (sgt.data_dir() / 'json' / r['file']).read_text(encoding='utf-8')
            self.assertEqual(text, pretty_json(sgt.load_group((r['order'], r['gap_id'])).to_dict()), r['file'])

    def test_dump_pretty_writes_a_file_that_make_finite_algebra_reads_back(self):
        g = sgt.load_group('D4') * sgt.load_group('Z3')
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, 'D4_x_Z3.json')
            dump_pretty(g, path)
            with open(path, encoding='utf-8') as fin:
                self.assertEqual(fin.read(), pretty_json(g.to_dict()))
            h = make_finite_algebra(path)
        self.assertIsInstance(h, Group)
        self.assertEqual(as_json_data(h.to_dict()), as_json_data(g.to_dict()))


class TestConstructionHelpers(TestCase):

    def test_mult_action_multiplies_by_a_power_of_k(self):
        act = mult_action(8, 5)
        self.assertEqual(act('0', '3'), '3')
        self.assertEqual(act('1', '3'), '7')       # 3 * 5 mod 8
        self.assertEqual(act('2', '3'), '3')       # 3 * 25 mod 8
        self.assertEqual(act('1', '0'), '0')

    def test_mult_action_builds_the_modular_group_of_order_16(self):
        z8, z2 = generate_cyclic_group(8), generate_cyclic_group(2)
        m16 = z8.semidirect_product(z2, mult_action(8, 5), name='M16')
        self.assertEqual(m16.order, 16)
        self.assertFalse(m16.is_commutative())
        self.assertEqual(table_and_elements(m16), table_and_elements(sgt.load_group('M16')))

    def test_mult_action_with_an_invalid_multiplier_is_rejected_by_semidirect_product(self):
        z8, z2 = generate_cyclic_group(8), generate_cyclic_group(2)
        with self.assertRaises(ValueError):
            z8.semidirect_product(z2, mult_action(8, 3 + 1))    # x -> 4x is not an automorphism of Z8

    def test_action_from_generators_builds_sl2_3_from_q8_and_z3(self):
        q8, z3 = sgt.load_group('Q8'), generate_cyclic_group(3)
        action = action_from_generators(q8, z3, {'1': {'a': 'ax', 'x': 'a'}})
        sl23 = q8.semidirect_product(z3, action, name='SL2_3')
        self.assertEqual(sl23.order, 24)
        self.assertEqual(table_and_elements(sl23), table_and_elements(sgt.load_group('SL2_3')))

    def test_action_from_generators_covers_every_element_of_both_groups(self):
        q8, z3 = sgt.load_group('Q8'), generate_cyclic_group(3)
        action = action_from_generators(q8, z3, {'1': {'a': 'ax', 'x': 'a'}})
        self.assertEqual(set(action), set(z3.elements))
        for h in z3.elements:
            self.assertEqual(set(action[h]), set(q8.elements))
            self.assertEqual(set(action[h].values()), set(q8.elements))   # each phi_h is a bijection

    def test_a_missing_image_means_the_identity_automorphism(self):
        z4, z2 = generate_cyclic_group(4), generate_cyclic_group(2)
        action = action_from_generators(z4, z2, {'1': {}})
        for h in z2.elements:
            self.assertEqual(action[h], {n: n for n in z4.elements})
        self.assertTrue(z4.semidirect_product(z2, action).is_commutative())     # the direct product

    def test_generators_that_do_not_generate_the_acting_group_are_rejected(self):
        z4, z2 = generate_cyclic_group(4), generate_cyclic_group(2)
        with self.assertRaises(ValueError):
            action_from_generators(z2, z4, {'2': {'1': '1'}})     # '2' generates only {0, 2} in Z4
