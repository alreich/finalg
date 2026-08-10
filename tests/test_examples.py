"""
Unit tests for finalg.examples.Examples, the convenience class that loads
example algebras listed in a JSON manifest (examples.json).
"""

import io
import os
import json
import shutil
import tempfile
import contextlib
from pathlib import Path
from unittest import TestCase

import finalg
from finalg.examples import Examples
from finalg.magma import Magma
from finalg.group import Group


ALGEBRAS_DIR = Path(__file__).resolve().parents[1] / 'src' / 'finalg' / 'data' / 'algebras'


class TestPrebuiltExamplesSingleton(TestCase):
    """finalg.examples is constructed automatically by finalg/__init__.py."""

    def test_is_an_examples_instance(self):
        self.assertIsInstance(finalg.examples, Examples)

    def test_len_matches_examples_json(self):
        with (ALGEBRAS_DIR / 'examples.json').open('r') as fin:
            names = json.load(fin)
        self.assertEqual(len(finalg.examples), len(names))

    def test_getitem_returns_finite_algebra(self):
        alg = finalg.examples[0]
        self.assertIsInstance(alg, Magma)  # Every FiniteAlgebra is-a Magma

    def test_getitem_indices_are_stable_and_distinct(self):
        first = finalg.examples[0]
        second = finalg.examples[1]
        self.assertNotEqual(first.name, second.name)
        # Fetching the same index twice returns the same (already-built) object.
        self.assertIs(finalg.examples[0], first)

    def test_contains_rock_paper_scissors_example(self):
        names = [finalg.examples[i].name for i in range(len(finalg.examples))]
        self.assertIn('RPS', names)

    def test_about_prints_all_examples_with_indices(self):
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            result = finalg.examples.about()
        self.assertIsNone(result)
        out = buf.getvalue()
        self.assertIn('Example Algebras', out)
        self.assertIn(f'{len(finalg.examples)} example algebras are available', out)
        self.assertIn('0: ', out)


class TestExamplesConstructionWithCustomManifest(TestCase):
    """Build a fresh, small Examples instance from a temporary directory to
    exercise the constructor and __len__/__getitem__/about directly, rather
    than relying solely on the pre-built finalg.examples singleton."""

    def setUp(self):
        self.tmpdir = tempfile.TemporaryDirectory()
        tmp_path = Path(self.tmpdir.name)
        self.selected_files = ['rock_paper_scissors.json', 'v4_klein_4_group.json']
        for filename in self.selected_files:
            shutil.copy(ALGEBRAS_DIR / filename, tmp_path / filename)
        with (tmp_path / 'my_examples.json').open('w') as fout:
            json.dump(self.selected_files, fout)
        self.examples = Examples(tmp_path, 'my_examples.json')

    def tearDown(self):
        self.tmpdir.cleanup()

    def test_len(self):
        self.assertEqual(len(self.examples), 2)

    def test_filenames_list_preserved_in_order(self):
        self.assertEqual(self.examples.filenames_list, self.selected_files)

    def test_getitem_order_matches_manifest_order(self):
        self.assertEqual(self.examples[0].name, 'RPS')
        self.assertEqual(self.examples[1].name, 'V4')

    def test_algebras_are_correct_types(self):
        self.assertIsInstance(self.examples[0], Magma)
        self.assertNotIsInstance(self.examples[0], Group)  # RPS is not a group
        self.assertIsInstance(self.examples[1], Group)  # V4 is a group

    def test_default_manifest_filename(self):
        # Examples defaults filenames_json to 'examples.json'; verify that
        # passing the default name explicitly gives an equivalent instance
        # when a real examples.json is present.
        default_examples = Examples(ALGEBRAS_DIR)
        self.assertEqual(len(default_examples), len(finalg.examples))

    def test_about_reports_correct_count_and_names(self):
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            self.examples.about()
        out = buf.getvalue()
        self.assertIn('2 example algebras are available', out)
        self.assertIn('0: RPS', out)
        self.assertIn('1: V4', out)


class TestExamplesConstructionErrors(TestCase):

    def test_missing_manifest_file_raises(self):
        with tempfile.TemporaryDirectory() as tmp:
            with self.assertRaises(FileNotFoundError):
                Examples(Path(tmp), 'nonexistent_manifest.json')

    def test_manifest_referencing_missing_algebra_file_raises(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            with (tmp_path / 'bad_manifest.json').open('w') as fout:
                json.dump(['does_not_exist.json'], fout)
            with self.assertRaises(FileNotFoundError):
                Examples(tmp_path, 'bad_manifest.json')
