# PERMUTATION CLASS

import random
from collections.abc import Sequence

class Perm(Sequence):

    def __init__(self, values):
        """A Perm represents a specific ordering of a finite set of n consecutive positive
        integers, where the minimum number ("base") is either 0 or 1. For example, a base-0
        Perm is an ordering (tuple) of the elements of the set {0,...,n-1}, and a base-1 Perm
        is an ordering of the elements of {1,...,n}. Same size Perms (n) with the same base
        can be composed ("multiplied") to obtain a new Perm. Composition is from right-to-left.
        Perms can also be inverted, using the method, inverse(), and they have a sign (+1 or -1),
        depending on whether their parity is even or odd, respectively. The Cauchy form is the
        default format for a Perm's internal representation, but they can be printed in cycle form,
        using the method, to_cycles(). Also, a Perm can be instantiated from a list of lists (cycles)
        using the static method, from_cycles.

        NOTE: If the base value (0 or 1) is a cycle, then it will be included in the cycle notation
        as (0,) or (1,), or in string form, '(0)' or '(1)' (i.e, without the comma). This also holds
        for cycles consisting of the maximum value, n or n-1.

        values: a list or tuple of integers representing the bottom line of Cauchy's 2-line notation.
        """
        self._values = tuple(values)
        self._base = min(values)

        if self._base not in (0, 1):
            raise ValueError("Minimum value must be 0 or 1.")

        # Normalize to 0-based indexing internally
        self._map = [v - self._base for v in self._values]

        # Validation check
        if sorted(self._map) != list(range(len(self._map))):
            raise ValueError(f"The input {values} is not a valid permutation.")

        self._is_even = None  # memoized at first call to 'is_even' method

    def __mul__(self, other):
        """Composes two permutations in right-to-left order. That is, when viewed
        as functions, (self * other)(x) = self(other(x))"""
        if len(self._map) != len(other.mapping):
            raise ValueError("Permutations must be of the same size to multiply.")
        if self._base != other.base:
            raise ValueError("Permutations must be of the same base to multiply.")
        product_values = [self._map[other.mapping[i]] + self._base for i in range(len(self._map))]
        return Perm(product_values)

    def __repr__(self):
        """Return a representation of the permutation that can be evaluated to create the same permutation."""
        return f"Perm({self._values})"

    def __str__(self):
        """Return a string representation of the permutation in cycle form."""
        s = str(self.to_cycles())[1:-1]
        return s.replace(',', '').replace(') (', ')(')

    def __len__(self):
        """Returns the size of the permutation."""
        return len(self._values)

    def __eq__(self, other):
        """Return True if this Perm equals the other Perm."""
        return self._values == other.values

    def __ne__(self, other):
        return not (self == other)

    def __hash__(self):
        """Allow the permutation to be hashable, e.g., used as dictionary key
        or member of a set."""
        return hash(self._values)

    def __getitem__(self, index):
        """Return the ith element of the permutation."""
        return self._values[index]

    def __call__(self, seq):
        """Apply the permutation to a sequence so that the sequence is rearranged according
        to the permutation's ordering. The length of the sequence has to equal the size of
        the permutation. If seq is a string, tuple, list, or even another Perm, then the
        return value will have the same type, resp. If seq is a range, then a list will be
        returned.
        Examples:
            p = Perm((4, 2, 1, 5, 3))
            p("ABCDE") ==> 'DBAEC'
            p((1, 2, 3, 4, 5)) ==> (4, 2, 1, 5, 3)
            p([1, 2, 3, 4, 5]) ==> [4, 2, 1, 5, 3]
            p(range(1, 6)) ==> [4, 2, 1, 5, 3]
            p(Perm((1, 2, 3, 4, 5))) ==> Perm((4, 2, 1, 5, 3))
        """
        return self._apply(seq)

    def _apply(self, seq: Sequence):
        if len(self) != len(seq):
            raise ValueError("The sequence must have the same length as the permutation.")
        q = [seq[i] for i in self.mapping]
        if isinstance(seq, Perm):
            return Perm(q)
        elif isinstance(seq, str):
            return "".join(q)
        elif isinstance(seq, tuple):
            return tuple(q)
        else:
            return q

    @property
    def values(self):
        """Return the 0-bases or 1-based values the permutation was constructed with."""
        return self._values

    @property
    def base(self):
        """Return the minimum of the permutation's values, either 0 or 1."""
        return self._base

    @property
    def mapping(self):
        """Internally, regardless of the input base, the permutation uses a 0-based range of
        numbers for some operations, such as multiplication. That is, mapping = (0,...,n-1),
        where n is the size of the permutation."""
        return self._map

    @property
    def size(self):
        """Returns same value as len(self). Redundant, but convenient."""
        return len(self)

    def inverse(self):
        """Returns the inverse of permutation."""
        inv_values = [0] * len(self)
        for i, val in enumerate(self._map):
            inv_values[val] = i + self._base
        return Perm(inv_values)

    def id(self):
        """Return the identity permutation for this permutation (i.e., same size and base)."""
        return Perm([i + self._base for i in range(len(self))])

    @staticmethod
    def identity(n, base=0):
        """Returns the identity permutation for n elements with the given base (default 0)."""
        if n > 0:
            return Perm([i + base for i in range(n)])
        else:
            raise ValueError("The size (n) must be > 0.")

    @staticmethod
    def random(n, base=0):
        """Generate a random permutation of size n with base = 0 (default) or 1."""
        if base not in (0, 1):
            raise ValueError("base must be 0 or 1")
        mapping = random.sample(range(n), k=n)
        values = [i + base for i in mapping]
        return Perm(values)

    @staticmethod
    def from_cycles(cycles, base=0):
        """Returns a Perm, given a list of lists (cycles) that represent a permutation
        in cycle format. This works for either 0-base or 1-base permutation cycles.
        Examples:
            Perm.from_cycles([[1, 3, 2], [4, 5]])  # ==> Perm((3, 1, 2, 5, 4))
            Perm.from_cycles([[0, 2, 1], [3, 4]])  # ==> Perm((2, 0, 1, 4, 3))
        """
        n = max(map(max, cycles))
        p = list(range(base, n + 1))
        for cycle in cycles:
            for i in range(len(cycle)):
                current_val = cycle[i]
                next_val = cycle[(i + 1) % len(cycle)]
                p[current_val - base] = next_val
        return Perm(p)

    def to_cycles(self):
        """Returns the cycle form of the permutation as a list of lists (cycles).
        This works for either 0-base or 1-base permutations in the default Cauchy form.
        Examples:
            Perm((3, 1, 2, 5, 4)).to_cycles()  # ==> [[1, 3, 2], [4, 5]]
            Perm((2, 0, 1, 4, 3)).to_cycles()  # ==> [[0, 2, 1], [3, 4]]
        """
        n = len(self)
        visited = [False] * n
        cycles = []
        for i in range(n):
            if not visited[i]:
                curr = i
                cycle = ()
                while not visited[curr]:
                    visited[curr] = True
                    # cycle.append(curr + self._base)
                    cycle = cycle + (curr + self._base,)
                    curr = self._map[curr]
                if len(cycle) > 1:
                    cycles.append(cycle)
                # Only keep a singleton if it contains the max value or min value
                # elif len(cycle) == 1 and cycle[0] == n - 1 + self._base:
                elif len(cycle) == 1 and (cycle[0] == n - 1 + self._base
                                          or cycle[0] == self._base):
                    cycles.append(cycle)
                else:
                    pass
        return cycles

    @property
    def is_even(self):
        """Returns True if the permutation is even or False if it is odd."""
        if self._is_even is None:
            self._is_even = self._is_even_fnc()
            return self._is_even
        else:
            return self._is_even

    def _is_even_fnc(self):
        """An internal method that checks if the permutation is even or odd."""
        inversions = 0
        # Iterate through all possible pairs of elements
        size = len(self)
        for i in range(size):
            for j in range(i + 1, size):
                # If a pair is out of natural order, it's an inversion
                if self._map[i] > self._map[j]:
                    inversions += 1
        # If the total inversions are even, the permutation is even
        return inversions % 2 == 0

    @property
    def sign(self):
        """Return +1 or -1 depending on whether the permutation is even or odd, respectively."""
        if self.is_even:
            return 1
        else:
            return -1

    @property
    def parity(self):
        """Returns the string "even" or "odd" depending on whether the permutation is even or odd."""
        if self.is_even:
            return "even"
        else:
            return "odd"
