"""
Unit tests for the loading functions of finalg.small_algebras, and for the serialized small groups
(all groups of order 1 through 40) that they read, which are stored as JSON files in the package
finalg.data.small_groups.

The tests need nothing but finalg: they do not use GAP. (The GAP SmallGroup ids that are stored in
index.json were checked when the library was built.)
"""

import itertools as it
import json
import os
import shutil
import tempfile
from pathlib import Path
from unittest import TestCase

# If the loader is ever moved or renamed, this import is the only line in this file that has to change.
import finalg.small_algebras as sgt

from finalg import (generate_cyclic_group, generate_dihedral_group,
                    generate_dicyclic_group, generate_symmetric_group, generate_alternating_group)
from finalg.group import Group

# Numbers of groups of order 1, 2, ..., 40 (OEIS A000001) ...
GROUP_COUNTS = [1, 1, 1, 2, 1, 2, 1, 5, 2, 2, 1, 5, 1, 2, 1, 14, 1, 5, 1, 5,
                2, 2, 1, 15, 2, 2, 5, 4, 1, 4, 1, 51, 1, 2, 1, 14, 1, 2, 2, 14]
# ... and numbers of abelian groups of those orders (OEIS A000688)
ABELIAN_COUNTS = [1, 1, 1, 2, 1, 1, 1, 3, 2, 1, 1, 2, 1, 1, 1, 5, 1, 2, 1, 2,
                  1, 1, 1, 3, 2, 1, 3, 2, 1, 1, 1, 7, 1, 1, 1, 4, 1, 1, 1, 3]
TOTAL_GROUPS = sum(GROUP_COUNTS)   # 181

# The groups of order 32 that are not constructed by finalg code, but were converted from GAP
GAP_GROUPS = [(32, 8), (32, 15), (32, 32)]


def as_json_data(d):
    """The dictionary as it looks after a trip through JSON, where tuples have become lists."""
    return json.loads(json.dumps(d))


def table_and_elements(alg):
    """What defines an algebra, apart from its name and description."""
    d = alg.to_dict()
    return d['elements'], d['table']


class TestDataPackage(TestCase):
    """The files that are shipped inside the package."""

    def test_the_data_folder_has_an_index_and_a_json_folder(self):
        folder = sgt.data_dir()
        self.assertTrue((folder / 'index.json').is_file())
        self.assertTrue((folder / 'json').is_dir())

    def test_every_file_in_the_index_exists_and_every_json_file_is_in_the_index(self):
        listed = {r['file'] for r in sgt.groups_index()}
        on_disk = {p.name for p in (sgt.data_dir() / 'json').iterdir() if p.name.endswith('.json')}
        self.assertEqual(listed, on_disk)
        self.assertEqual(len(on_disk), TOTAL_GROUPS)

    def test_pyproject_ships_the_data(self):
        # This guards against forgetting the package-data entry, without having to build a wheel.
        pyproject = Path(__file__).resolve().parent.parent / 'pyproject.toml'
        if not pyproject.is_file():
            self.skipTest("pyproject.toml is not available")
        text = pyproject.read_text(encoding='utf-8')
        self.assertIn('finalg.data.small_groups', text)
        self.assertIn('json/*.json', text)

    def test_the_loader_does_not_depend_on_the_working_directory(self):
        here = os.getcwd()
        with tempfile.TemporaryDirectory() as elsewhere:
            os.chdir(elsewhere)
            try:
                self.assertEqual(sgt.load_group('Q8').order, 8)
                self.assertEqual(len(sgt.groups_index()), TOTAL_GROUPS)
            finally:
                os.chdir(here)


class TestGroupsIndex(TestCase):

    def test_there_is_one_record_for_every_group_of_order_1_through_40(self):
        index = sgt.groups_index()
        self.assertEqual(len(index), TOTAL_GROUPS)
        self.assertEqual(TOTAL_GROUPS, 181)
        for n in range(1, 41):
            self.assertEqual(sum(1 for r in index if r['order'] == n), GROUP_COUNTS[n - 1], f"order {n}")

    def test_records_have_the_documented_keys(self):
        for r in sgt.groups_index():
            self.assertEqual(set(r), {'order', 'gap_id', 'name', 'file', 'structure', 'abelian', 'source'})
            self.assertIn(r['source'], ('finalg', 'gap'))
            self.assertIsInstance(r['abelian'], bool)

    def test_records_are_sorted_by_order_and_smallgroup_index_without_gaps(self):
        index = sgt.groups_index()
        self.assertEqual([(r['order'], r['gap_id']) for r in index],
                         [(n, i) for n in range(1, 41) for i in range(1, GROUP_COUNTS[n - 1] + 1)])

    def test_names_ids_and_files_are_unique(self):
        index = sgt.groups_index()
        for key in ('name', 'file'):
            values = [r[key] for r in index]
            self.assertEqual(len(values), len(set(values)), key)

    def test_file_names_follow_the_naming_convention(self):
        for r in sgt.groups_index():
            self.assertEqual(r['file'], f"G{r['order']:02d}_{r['gap_id']:02d}_{r['name']}.json")

    def test_number_of_abelian_groups_of_each_order(self):
        index = sgt.groups_index()
        for n in range(1, 41):
            self.assertEqual(sum(r['abelian'] for r in index if r['order'] == n), ABELIAN_COUNTS[n - 1],
                             f"order {n}")

    def test_exactly_three_groups_come_from_gap(self):
        self.assertEqual([(r['order'], r['gap_id']) for r in sgt.groups_index() if r['source'] == 'gap'],
                         GAP_GROUPS)

    def test_changing_the_returned_list_does_not_change_the_index(self):
        index = sgt.groups_index()
        index.pop()
        index[0]['name'] = 'changed'
        fresh = sgt.groups_index()
        self.assertEqual(len(fresh), TOTAL_GROUPS)
        self.assertEqual(fresh[0]['name'], 'Z1')


class TestLoadGroup(TestCase):

    def test_by_name(self):
        for name, order in (('Z1', 1), ('Q8', 8), ('F21', 21), ('SL2_3', 24), ('S4', 24), ('Heis3', 27),
                            ('Q32', 32), ('D20', 40)):
            g = sgt.load_group(name)
            self.assertIsInstance(g, Group)
            self.assertEqual((g.name, g.order), (name, order))

    def test_by_smallgroup_id_as_a_tuple_or_a_list(self):
        self.assertEqual(sgt.load_group((8, 4)).name, 'Q8')
        self.assertEqual(sgt.load_group([8, 4]).name, 'Q8')
        self.assertEqual(sgt.load_group((24, 12)).name, 'S4')
        self.assertEqual(sgt.load_group((32, 8)).name, 'SG32_8')

    def test_by_file_stem_or_file_name(self):
        self.assertEqual(sgt.load_group('G08_04').name, 'Q8')
        self.assertEqual(sgt.load_group('G08_04_Q8').name, 'Q8')
        self.assertEqual(sgt.load_group('G08_04_Q8.json').name, 'Q8')

    def test_by_path_of_a_json_file(self):
        path = sgt.data_dir() / 'json' / 'G08_04_Q8.json'
        self.assertEqual(sgt.load_group(Path(str(path))).name, 'Q8')
        self.assertEqual(sgt.load_group(str(path)).name, 'Q8')

    def test_every_way_of_naming_a_group_gives_the_same_group(self):
        reference = table_and_elements(sgt.load_group('Z3_sdp_Z8'))
        for key in ((24, 1), 'G24_01', 'G24_01_Z3_sdp_Z8', 'G24_01_Z3_sdp_Z8.json'):
            self.assertEqual(table_and_elements(sgt.load_group(key)), reference, key)

    def test_each_call_returns_a_new_object(self):
        self.assertIsNot(sgt.load_group('Q8'), sgt.load_group('Q8'))

    def test_a_loaded_group_works_with_the_rest_of_finalg(self):
        q8 = sgt.load_group('Q8')
        self.assertFalse(q8.is_commutative())
        self.assertEqual(sorted(q8.element_order(e) for e in q8.elements), [1, 2, 4, 4, 4, 4, 4, 4])
        self.assertEqual(len(q8.center()), 2)

    def test_unknown_names_and_ids_raise_key_error(self):
        for bad in ('Nope', 'G99_01', (99, 1), (8, 6), (0, 1)):
            with self.assertRaises(KeyError, msg=repr(bad)):
                sgt.load_group(bad)

    def test_the_error_for_an_unknown_name_suggests_similar_names(self):
        with self.assertRaises(KeyError) as ctx:
            sgt.load_group('Dic')
        self.assertIn('Dic3', str(ctx.exception))

    def test_keys_of_the_wrong_type_raise_type_error(self):
        for bad in (8, 3.5, None, (1, 2, 3), {'name': 'Q8'}):
            with self.assertRaises(TypeError, msg=repr(bad)):
                sgt.load_group(bad)


class TestLoadGroups(TestCase):

    def test_all_groups_by_default(self):
        groups = sgt.load_groups()
        self.assertEqual(len(groups), TOTAL_GROUPS)
        self.assertEqual([(g.order, g.name) for g in groups],
                         [(r['order'], r['name']) for r in sgt.groups_index()])

    def test_by_order(self):
        self.assertEqual([g.name for g in sgt.load_groups(orders=4)], ['Z4', 'V4'])
        self.assertEqual(len(sgt.load_groups(orders=16)), 14)
        self.assertEqual(len(sgt.load_groups(orders=32)), 51)
        self.assertEqual({g.order for g in sgt.load_groups(orders=[6, 8])}, {6, 8})
        self.assertEqual(len(sgt.load_groups(orders=range(1, 21))), 54)
        self.assertEqual(len(sgt.load_groups(orders=range(21, 41))), 127)

    def test_orders_without_groups_give_an_empty_result(self):
        self.assertEqual(sgt.load_groups(orders=41), [])
        self.assertEqual(sgt.load_groups(orders=[]), [])

    def test_by_name(self):
        groups = sgt.load_groups(names=['Q8', 'Pauli', 'F21'])
        self.assertEqual([g.name for g in groups], ['Q8', 'Pauli', 'F21'])   # in library order, not request order
        self.assertEqual([g.name for g in sgt.load_groups(names='Q8')], ['Q8'])

    def test_unknown_names_raise_key_error(self):
        with self.assertRaises(KeyError):
            sgt.load_groups(names=['Q8', 'Nope'])

    def test_orders_and_names_must_both_match(self):
        self.assertEqual([g.name for g in sgt.load_groups(orders=8, names=['Q8', 'F21'])], ['Q8'])

    def test_by_keys_in_the_order_given(self):
        groups = sgt.load_groups(keys=[(32, 8), 'D4', 'G24_03', 'Z2'])
        self.assertEqual([g.name for g in groups], ['SG32_8', 'D4', 'SL2_3', 'Z2'])

    def test_as_a_dictionary_keyed_by_smallgroup_id(self):
        d = sgt.load_groups(orders=[6, 8], as_dict=True)
        self.assertEqual(sorted(d), [(6, 1), (6, 2), (8, 1), (8, 2), (8, 3), (8, 4), (8, 5)])
        self.assertEqual(d[(8, 4)].name, 'Q8')
        everything = sgt.load_groups(as_dict=True)
        self.assertEqual(len(everything), TOTAL_GROUPS)
        self.assertEqual(everything[(32, 20)].name, 'Q32')


class TestEveryGroupInTheLibrary(TestCase):

    @classmethod
    def setUpClass(cls):
        cls.index = sgt.groups_index()
        cls.groups = sgt.load_groups()

    def test_every_algebra_is_a_group_with_the_right_order_and_name(self):
        for r, g in zip(self.index, self.groups):
            self.assertIsInstance(g, Group, r['file'])
            self.assertEqual(g.order, r['order'], r['file'])
            self.assertEqual(g.name, r['name'], r['file'])

    def test_the_abelian_flag_in_the_index_is_correct(self):
        for r, g in zip(self.index, self.groups):
            self.assertEqual(g.is_commutative(), r['abelian'], r['file'])

    def test_every_group_has_an_identity_and_inverses_and_is_associative(self):
        for r, g in zip(self.index, self.groups):
            self.assertTrue(g.has_identity(), r['file'])
            self.assertTrue(g.has_inverses(), r['file'])
            self.assertTrue(g.is_associative(), r['file'])

    def test_the_identity_is_the_first_element(self):
        for r, g in zip(self.index, self.groups):
            self.assertEqual(g.identity, g.elements[0], r['file'])

    def test_element_names_are_unique_strings(self):
        for r, g in zip(self.index, self.groups):
            self.assertTrue(all(isinstance(e, str) for e in g.elements), r['file'])
            self.assertEqual(len(set(g.elements)), g.order, r['file'])

    def test_groups_of_order_20_or_less_are_pairwise_non_isomorphic(self):
        for n in range(1, 21):
            for g1, g2 in it.combinations([g for g in self.groups if g.order == n], 2):
                self.assertFalse(g1.fast_isomorphic(g2), f"{g1.name} and {g2.name} are isomorphic")

    def test_groups_of_the_same_order_differ_in_simple_invariants_or_are_few(self):
        # A cheap sanity check for the large orders, where an isomorphism test is slow: no two groups
        # of one order are *identical* as Cayley tables (the tables themselves would then be equal).
        for n in range(21, 41):
            tables = [json.dumps(table_and_elements(g)) for g in self.groups if g.order == n]
            self.assertEqual(len(tables), len(set(tables)), f"order {n}")

    def test_the_groups_that_came_from_gap_have_order_32_and_are_not_abelian(self):
        for gid in GAP_GROUPS:
            g = sgt.load_group(gid)
            self.assertEqual(g.order, 32)
            self.assertFalse(g.is_commutative())
            self.assertEqual(g.name, f"SG{gid[0]}_{gid[1]}")


class TestAgreementWithFinalgGenerators(TestCase):
    """Groups in the library that finalg can also generate by itself must be the very same Cayley tables."""

    def assert_same_group(self, generated, name=None):
        loaded = sgt.load_group(name or generated.name)
        self.assertEqual(table_and_elements(generated), table_and_elements(loaded), generated.name)

    def test_cyclic_groups(self):
        for n in range(1, 41):
            self.assert_same_group(generate_cyclic_group(n))

    def test_dihedral_groups_of_order_6_through_40(self):
        for n in range(3, 21):
            self.assert_same_group(generate_dihedral_group(n)[0])

    def test_dicyclic_groups_of_order_8_through_40(self):
        for n in range(2, 11):
            self.assert_same_group(generate_dicyclic_group(n)[0])

    def test_symmetric_and_alternating_groups(self):
        self.assert_same_group(generate_symmetric_group(4))
        self.assert_same_group(generate_alternating_group(4))

    def test_s3_is_isomorphic_to_d3(self):
        # The elements of S3 are named differently from those of D3, so only isomorphism can be tested
        self.assertTrue(generate_symmetric_group(3).fast_isomorphic(sgt.load_group('D3')))


class TestAnotherDirectory(TestCase):
    """The `directory` argument, and verify_library."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.copy = Path(self.tmp.name) / 'library'
        shutil.copytree(str(sgt.data_dir()), str(self.copy), ignore=shutil.ignore_patterns('__pycache__'))

    def edit_index(self, change):
        path = self.copy / 'index.json'
        index = json.loads(path.read_text(encoding='utf-8'))
        change(index)
        path.write_text(json.dumps(index), encoding='utf-8')

    def test_a_copy_of_the_library_loads_like_the_original(self):
        self.assertEqual(len(sgt.groups_index(self.copy)), TOTAL_GROUPS)
        self.assertEqual(table_and_elements(sgt.load_group('Q8', directory=self.copy)),
                         table_and_elements(sgt.load_group('Q8')))
        self.assertEqual(len(sgt.load_groups(orders=[4, 6], directory=str(self.copy))), 4)

    def test_verify_library_accepts_the_real_library(self):
        self.assertEqual(sgt.verify_library(verbose=False), TOTAL_GROUPS)
        self.assertEqual(sgt.verify_library(self.copy, verbose=False), TOTAL_GROUPS)

    def test_verify_library_notices_a_wrong_name(self):
        self.edit_index(lambda index: index[3].update(name='Wrong'))
        with self.assertRaises(AssertionError):
            sgt.verify_library(self.copy, verbose=False)

    def test_verify_library_notices_a_wrong_order(self):
        self.edit_index(lambda index: index[10].update(order=12))
        with self.assertRaises(AssertionError):
            sgt.verify_library(self.copy, verbose=False)

    def test_verify_library_notices_a_missing_group(self):
        self.edit_index(lambda index: index.pop(4))     # Z5 is the only group of order 5
        with self.assertRaises(AssertionError):
            sgt.verify_library(self.copy, verbose=False)

    def test_verify_library_notices_a_missing_file(self):
        (self.copy / 'json' / 'G08_04_Q8.json').unlink()
        with self.assertRaises(FileNotFoundError):
            sgt.verify_library(self.copy, verbose=False)

    def test_changes_to_a_copy_are_seen_at_once(self):
        # The index of a user-supplied folder is not cached, but read again at every call
        self.assertEqual(sgt.groups_index(self.copy)[3]['name'], 'Z4')
        self.edit_index(lambda index: index[3].update(name='Renamed'))
        self.assertEqual(sgt.groups_index(self.copy)[3]['name'], 'Renamed')
        self.assertEqual(sgt.groups_index()[3]['name'], 'Z4')
