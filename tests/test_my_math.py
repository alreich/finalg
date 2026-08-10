"""
Unit tests for finalg.my_math: number-theoretic utility functions.
"""

from unittest import TestCase

from finalg.my_math import is_relatively_prime, relative_primes, totient, divisors, xgcd


class TestIsRelativelyPrime(TestCase):

    def test_relatively_prime_pair(self):
        self.assertTrue(is_relatively_prime(8, 15))

    def test_non_relatively_prime_pair(self):
        self.assertFalse(is_relatively_prime(8, 4))

    def test_one_is_relatively_prime_to_everything(self):
        self.assertTrue(is_relatively_prime(1, 100))

    def test_equal_numbers_greater_than_1_not_relatively_prime(self):
        self.assertFalse(is_relatively_prime(6, 6))


class TestRelativePrimes(TestCase):

    def test_relative_primes_of_8(self):
        self.assertEqual(relative_primes(8), [1, 3, 5, 7])

    def test_relative_primes_of_prime_number(self):
        self.assertEqual(relative_primes(7), [1, 2, 3, 4, 5, 6])

    def test_relative_primes_of_1(self):
        # range(1) is just [0], and gcd(0, 1) == 1, so 0 counts as relatively prime to 1.
        self.assertEqual(relative_primes(1), [0])


class TestTotient(TestCase):

    def test_totient_of_8(self):
        self.assertEqual(totient(8), 4)

    def test_totient_of_prime(self):
        self.assertEqual(totient(7), 6)

    def test_totient_matches_len_of_relative_primes(self):
        for n in [5, 12, 20]:
            self.assertEqual(totient(n), len(relative_primes(n)))


class TestDivisors(TestCase):

    def test_non_trivial_divisors_of_12(self):
        self.assertEqual(divisors(12), [2, 3, 4, 6])

    def test_all_divisors_of_12(self):
        self.assertEqual(divisors(12, non_trivial=False), [1, 2, 3, 4, 6, 12])

    def test_prime_has_no_non_trivial_divisors(self):
        self.assertEqual(divisors(7), [])

    def test_prime_all_divisors(self):
        self.assertEqual(divisors(7, non_trivial=False), [1, 7])


class TestXgcd(TestCase):

    def test_gcd_value(self):
        g, x, y = xgcd(240, 46)
        self.assertEqual(g, 2)

    def test_bezout_identity_holds(self):
        a, b = 240, 46
        g, x, y = xgcd(a, b)
        self.assertEqual(g, a * x + b * y)

    def test_coprime_numbers_give_gcd_1(self):
        g, x, y = xgcd(17, 5)
        self.assertEqual(g, 1)
        self.assertEqual(17 * x + 5 * y, 1)

    def test_one_argument_zero(self):
        g, x, y = xgcd(5, 0)
        self.assertEqual(g, 5)
        self.assertEqual(5 * x + 0 * y, 5)
