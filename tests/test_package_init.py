"""
Unit tests for finalg's package-level public API, as defined by __init__.py's
__all__ list. These are integration-style tests over the package surface
itself, rather than any single algebra-class module.
"""

from unittest import TestCase


class TestStarImportExposesPublicApi(TestCase):
    """`from finalg import *` should bring every name in __all__ into scope
    and each should be usable."""

    def test_star_import_and_use_each_generator(self):
        ns = {}
        exec("from finalg import *", ns)

        # Spot-check a representative subset of __all__ actually works end to end.
        rps = ns['make_finite_algebra']('RPS', 'Rock, Paper, Scissors',
                                        ['r', 'p', 's'], [[0, 1, 0], [1, 1, 2], [0, 2, 2]])
        self.assertEqual(rps.elements, ('r', 'p', 's'))

        z5 = ns['generate_cyclic_group'](5)
        self.assertEqual(z5.order, 5)

        s3 = ns['generate_symmetric_group'](3)
        self.assertEqual(s3.order, 6)

        ps3 = ns['generate_powerset_group'](3)
        self.assertEqual(ps3.order, 8)

        m4 = ns['generate_commutative_monoid'](4)
        self.assertEqual(m4.order, 4)

        rp = ns['generate_relative_primes_group'](8)
        self.assertEqual(rp.order, 4)

        psr = ns['generate_powerset_ring'](2)
        self.assertEqual(psr.order, 4)

        f7 = ns['generate_algebra_mod_n'](7)
        self.assertEqual(f7.order, 7)

        d8, _, _ = ns['generate_dihedral_group'](4)
        self.assertEqual(d8.order, 8)

        table = ns['make_cayley_table']([[0, 1], [1, 0]], ('a', 'b'))
        self.assertEqual(table.order, 2)

        with ns['InfixNotation'](z5) as f:
            self.assertEqual((f['1'] + f['2']).name, '3')

        p = ns['Perm']((1, 0, 2))
        self.assertEqual(p.values, (1, 0, 2))

        self.assertGreater(len(ns['examples']), 0)

    def test_all_names_are_importable(self):
        import finalg
        for name in finalg.__all__:
            self.assertTrue(hasattr(finalg, name), f"finalg.{name} is missing")
